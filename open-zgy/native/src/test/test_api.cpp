// Copyright 2017-2021, Schlumberger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "test_all.h"
#include "test_utils.h"
#include "test_sdutils.h"
#include "../api.h"
#include "../writerargs.h"
#include "../iocontext.h"
#include "../exception.h"
#include "../impl/environment.h"
#include "../impl/mtguard.h"
#include "../impl/timer.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <memory>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <limits>
#include <thread>
#include <chrono>
#include <algorithm>
#include <functional>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
#include <unistd.h>
#endif

using namespace OpenZGY;
using namespace OpenZGY::Formatters;
using Test_Utils::LocalFileAutoDelete;
using Test_Utils::CloudFileAutoDelete;
using Test_Utils::must_throw;

namespace {
  template<typename T, std::size_t  N>
  std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a)
  {
    os << "[";
    for (std::size_t ii=0; ii<N; ++ii)
      os << a[ii] << (ii == N-1 ? "" : ", ");
    os << "]";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const SampleStatistics& in)
  {
    os << "cnt: " << in.cnt
       << " sum: " << in.sum
       << " ssq: " << in.ssq
       << " min: " << in.min
       << " max: " << in.max;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const SampleHistogram& in)
  {
    os << "cnt: " << in.samplecount
       << " min: " << in.minvalue
       << " max: " << in.maxvalue
       << " bincount: " << in.bins.size();
    return os;
  }
}

namespace Test_API {
#if 0
}
#endif

/**
 * This is a messy hack to see if we are running under valgring and must
 * expect tests to take a lot longer than usual. Possibly so long that
 * the tests aren't feasible to run.
 *
 * The official way of doing this is RUNNING_ON_VALGRIND macro,
 * but I don't want a dependency on <valgrind.h>, and making a local
 * copy of that file is just as messy.
 */
static bool
is_running_on_valgrind()
{
  std::string p = InternalZGY::Environment::getStringEnv("LD_LIBRARY_PATH", "");
  return (strstr (p.c_str(), "/valgrind/") != nullptr ||
          strstr (p.c_str(), "/vgpreload") != nullptr);
}

static std::string
get_testdata(const std::string& name)
{
  using InternalZGY::Environment;
#ifdef _WIN32
  std::string result = Environment::getStringEnv("OPENZGY_TESTDATA", "..\\..\\build\\testdata");
  if (result.back() != '\\')
    result += "\\";
  result += name;
  return result;
#else
  return
    Environment::getStringEnv("OPENZGY_TESTDATA", "../../build/testdata")
    + "/" + name;
#endif
}

static bool
similar(double a, double b, double eps)
{
  return std::abs(a - b)  <= eps * 0.5 * (std::abs(a) + std::abs(b));
}

static void
dump_api(std::shared_ptr<OpenZGY::IZgyTools> rr, std::ostream& out)
{
  const OpenZGY::IZgyTools& r = *rr.get();
  std::streamsize oldprec = std::cout.precision();
  std::ios_base::fmtflags oldflags = std::cout.flags();
  out << "File format and version        = " << r.datatype()
      << " ZGY version " << r.filestats()->fileVersion() << "\n";
  out << "Size I,J,K                     = " << r.size() << "\n";
  out << "Brick size I,J,K               = " << r.bricksize() << "\n";
  out << "Number of bricks I,J,K         = " << r.brickcount()[0] << "\n";
  out << "Number of LODs                 = " << r.nlods() << "\n";
  out << "Coding range min/max           = " << r.datarange() << "\n";
  out << "Statistical min/max/count      = " << r.statistics() << "\n";
  out << "Histogram range min/max/count  = " << r.histogram() << "\n";
  out << "Inline start/increment/count   = "
      << r.annotstart()[0] << " "
      << r.annotinc()[0] << " "
      << r.size()[0] << "\n";
  out << "Xline  start/increment/count   = "
      << r.annotstart()[1] << " "
      << r.annotinc()[1] << " "
      << r.size()[1] << "\n";
  out << "Sample start/increment/count   = "
      << r.zstart() << " "
      << r.zinc() << " "
      << r.size()[2] << "\n";
  out << "Horizontal dim/factor/name     = "
      << r.hunitdim() << " "
      << r.hunitfactor() << " '"
      << r.hunitname() << "'\n";
  out << "Vertical dim/factor/name       = "
      << r.zunitdim() << " "
      << r.zunitfactor() << " '"
      << r.zunitname() << "'\n";
  out << "Ordered Corner Points Legend   = [  <i>,   <j>] { <inline>,   <xline>} (  <easting>,  <northing>)" << "\n";
  r.filestats()->dump(out);
  for (int ii=0; ii<4; ++ii)
    out << "Ordered Corner Point " << ii << "         = ["
        << std::fixed << std::setprecision(0)
        << std::setw(5) << r.indexcorners()[ii][0] << ", "
        << std::setw(5) << r.indexcorners()[ii][1] << "] {"
        << std::scientific << std::setprecision(6)
        << std::setw(9) << r.annotcorners()[ii][0] << ", "
        << std::setw(9) << r.annotcorners()[ii][1] << "} ("
        << std::fixed << std::setprecision(2)
        << std::setw(11) << r.corners()[ii][0] << ", "
        << std::setw(11) << r.corners()[ii][1] << ")\n"
        << std::scientific << std::setprecision(oldprec);
  // Note, I'd like defaultfloat instead of scientific
  // but the former is not supported on all compilers.
  // This is also why I need to reset the flags.
  std::cout.flags(oldflags);
}

/**
 * Read metadata from the well-known "Empty-v3.zgy" file.
 * Either directly from the reader, or from a writer that
 * is about to copy the file.
 */
static void
do_test_readmeta(const OpenZGY::IZgyTools& r, bool readonly, int version)
{
  TEST_CHECK(r.datatype() == SampleDataType::int8);
  TEST_EQUAL(r.size()[0], 181);
  TEST_EQUAL(r.size()[1], 241);
  TEST_EQUAL(r.size()[2], 169);
  TEST_EQUAL(r.bricksize()[0], 64);
  TEST_EQUAL(r.bricksize()[1], 64);
  TEST_EQUAL(r.bricksize()[2], 64);
  TEST_EQUAL(r.brickcount()[0][0], 3);
  TEST_EQUAL(r.brickcount()[0][1], 4);
  TEST_EQUAL(r.brickcount()[0][2], 3);
  TEST_EQUAL(r.nlods(), (readonly ? 3 : 1));
  TEST_EQUAL_FLOAT(r.datarange()[0], -10038.5, 0.5);
  TEST_EQUAL_FLOAT(r.datarange()[1], +9761.62, 0.5);
  if (readonly) {
    // Skip these tests if still writing the file.
    TEST_EQUAL_FLOAT(r.statistics().min, -10038.5, 0.5);
    TEST_EQUAL_FLOAT(r.statistics().max, +9761.62, 0.5);
    TEST_EQUAL_FLOAT(r.histogram().minvalue, -10038.5, 0.5);
    TEST_EQUAL_FLOAT(r.histogram().maxvalue, +9761.62, 0.5);
    TEST_EQUAL(r.statistics().cnt, (version==1 ? 0 : 7371949));
    TEST_EQUAL(r.histogram().samplecount, (version==1 ? 0 : 7371949));
  }
  else {
    TEST_EQUAL(r.statistics().min, 0);
    TEST_EQUAL(r.statistics().max, 0);
    TEST_EQUAL(r.histogram().minvalue, 0);
    TEST_EQUAL(r.histogram().maxvalue, 0);
    TEST_EQUAL(r.statistics().cnt, 0);
    TEST_EQUAL(r.histogram().samplecount, 0);
  }
  TEST_EQUAL(r.annotstart()[0], 640);
  TEST_EQUAL(r.annotstart()[1], 920);
  TEST_EQUAL(r.annotinc()[0], 2);
  TEST_EQUAL(r.annotinc()[1], 2);
  TEST_EQUAL(r.zstart(), 648);
  TEST_EQUAL(r.zinc(), 6);
  TEST_CHECK(r.hunitdim() == UnitDimension::unknown);
  TEST_EQUAL(r.hunitfactor(), 1.0);
  TEST_EQUAL(r.hunitname(), "");
  TEST_CHECK(r.zunitdim() == UnitDimension::unknown);
  TEST_EQUAL(r.zunitfactor(), 1.0);
  TEST_EQUAL(r.zunitname(), "");
  std::shared_ptr<const FileStatistics>filestats = r.filestats();
  TEST_EQUAL(filestats->fileVersion(), version);
  if (readonly) {
    // Skip these tests of wtill writing the file.
    TEST_EQUAL(filestats->fileSize(), (version==1 ? 11867826 : 12320768));
    TEST_EQUAL(filestats->alphaNormalCount(), (version==1 ? 17 : 5));
    TEST_EQUAL(filestats->alphaNormalSizePerEntry(), 64*64);
    TEST_EQUAL(filestats->brickNormalCount(), 45);
    TEST_EQUAL(filestats->brickNormalSizePerEntry(), 64*64*64);
  }

  // actual
  const IZgyReader::corners_t& index = r.indexcorners();
  const IZgyReader::corners_t& annot = r.annotcorners();
  const IZgyReader::corners_t& world = r.corners();

  // expect
  double ibeg[2] { 0, 0 };
  double iend[2] { (double)r.size()[0] - 1, (double)r.size()[1] - 1 };
  double abeg[2] { r.annotstart()[0], r.annotstart()[1] };
  double aend[2] { abeg[0] + r.annotinc()[0] * (r.size()[0] - 1),
                   abeg[1] + r.annotinc()[1] * (r.size()[1] - 1)};

  TEST_CHECK(index[0][0] == ibeg[0] && index[0][1] == ibeg[1]);
  TEST_CHECK(index[1][0] == iend[0] && index[1][1] == ibeg[1]);
  TEST_CHECK(index[2][0] == ibeg[0] && index[2][1] == iend[1]);
  TEST_CHECK(index[3][0] == iend[0] && index[3][1] == iend[1]);

  TEST_CHECK(annot[0][0] == abeg[0] && annot[0][1] == abeg[1]);
  TEST_CHECK(annot[1][0] == aend[0] && annot[1][1] == abeg[1]);
  TEST_CHECK(annot[2][0] == abeg[0] && annot[2][1] == aend[1]);
  TEST_CHECK(annot[3][0] == aend[0] && annot[3][1] == aend[1]);

  TEST_EQUAL_FLOAT(world[0][0],  564344.97, 0.1);
  TEST_EQUAL_FLOAT(world[0][1], 5917369.23, 0.1);
  TEST_EQUAL_FLOAT(world[1][0],  568209.86, 0.1);
  TEST_EQUAL_FLOAT(world[1][1], 5915035.21, 0.1);
  TEST_EQUAL_FLOAT(world[2][0],  567457.07, 0.1);
  TEST_EQUAL_FLOAT(world[2][1], 5922522.52, 0.1);
  TEST_EQUAL_FLOAT(world[3][0],  571321.95, 0.1);
  TEST_EQUAL_FLOAT(world[3][1], 5920188.49, 0.1);

  std::string nulluuid = "00000000-0000-0000-0000-000000000000";
  TEST_EQUAL(r.verid().size(), nulluuid.size());
  bool verid_missing = (r.verid() == nulluuid);
  TEST_CHECK(version == 1 ? verid_missing : !verid_missing);
}

static void
test_dump()
{
  std::string fname = get_testdata("Empty-v3.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  std::stringstream ss;
  reader->dump(ss);
  if (verbose()) {
    std::cout << "\n" << ss.str() << std::flush;
  }
  reader->close();
}

#if 0 // Locking is disabled, set to complain-only mode.
static void
test_locks()
{
  LocalFileAutoDelete lad("testlocks.zgy");
  ZgyWriterArgs args = ZgyWriterArgs().filename(lad.name()).size(64,64,128);
  ZgyWriterArgs rargs = ZgyWriterArgs().filename(lad.name());
  std::shared_ptr<IZgyWriter> writer;
  std::shared_ptr<IZgyReader> reader1, reader2;

  // Read / write / update disallowed on a file already open for write.
  writer = IZgyWriter::open(args);
  must_throw("Already opened for write", [&](){
    IZgyReader::open(lad.name());});
  must_throw("Already opened for write", [&](){
    IZgyWriter::open(args);});
  must_throw("Already opened for write", [&](){
    IZgyWriter::reopen(rargs);});
  writer->close();

  // Read / overwrite / update is good when nothing else is open.
  reader1 = IZgyReader::open(lad.name());
  reader1->close();
  writer = IZgyWriter::open(args);
  writer->close();
  writer = IZgyWriter::reopen(rargs);
  writer->close();

  // Multiple open for read is allowed.
  reader1 = IZgyReader::open(lad.name());
  reader2 = IZgyReader::open(lad.name());

  // Open for overwrite or update only allowed once all readers are closed.
  must_throw("Already opened for read", [&](){
    IZgyWriter::reopen(rargs);});
  reader1->close();
  must_throw("Already opened for read", [&](){
    IZgyWriter::reopen(rargs);});
  reader2->close();
  writer = IZgyWriter::reopen(rargs);
  writer->close();

  // Letting a reader go out of scope will unlock it, in spite of it
  // not being completely closed.
  reader1 = IZgyReader::open(lad.name());
  reader1.reset();
  writer = IZgyWriter::reopen(rargs);
  writer->close();

  // Letting a writer go out of scope will close and unlock it.
  writer = IZgyWriter::reopen(rargs);
  writer.reset();
  writer = IZgyWriter::reopen(rargs);
  writer->close();
}
#endif

static void
test_readmeta_r()
{
  std::string fname = get_testdata("Empty-v3.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  if (verbose()) {
    std::cout << "\n";
    dump_api(reader, std::cout);
    reader->filestats()->dump(std::cout, "filestats: ");
  }
  do_test_readmeta(*reader, /*readonly=*/true, 3);
  reader->close();
}

static void
test_readmeta_w()
{
  std::string fname = get_testdata("Empty-v3.zgy");
  LocalFileAutoDelete lad("testfile.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  ZgyWriterArgs args = ZgyWriterArgs()
    .metafrom(reader)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  if (verbose()) {
    std::cout << "\n";
    dump_api(writer, std::cout);
    writer->filestats()->dump(std::cout, "filestats: ");
  }
  // After July 2023, a completely empty file allows reading lod data.
  // After November 2023, writable files always have nlods()==1.
  do_test_readmeta(*writer, /*readonly=*/false, 3);
  writer->finalize();
  writer->close();
  reader->close();
}

static void
test_readmeta_v1_r()
{
  std::string fname = get_testdata("Empty-v1.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  if (verbose()) {
    std::cout << "\n";
    dump_api(reader, std::cout);
    reader->filestats()->dump(std::cout, "filestats: ");
  }
  do_test_readmeta(*reader, /*readonly=*/true, 1);
  reader->close();
}

#if 0 // The required data file is not checked in yet.
static void
test_readcmeta()
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Compressed.zgy"));
  if (verbose()) {
    std::cout << "\n";
    dump_api(reader, std::cout);
    reader->filestats()->dump(std::cout, "filestats: ");
  }
}
#endif

static void
test_readconst()
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v3.zgy"));
  std::pair<bool,double> c = reader->readconst(std::array<std::int64_t,3>{0,0,0}, reader->size(), 0, true);
  //std::cout << "constant? " << std::boolalpha << c.first << std::noboolalpha << " value " << c.second << "\n";
  TEST_CHECK(!c.first);
  reader->close();
}

static void
test_readbulk()
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v3.zgy"));
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size = reader->size();
  std::unique_ptr<float[]>buf(new float[size[0] * size[1] * size[2]]);
  memset(buf.get(), 0, size[0] * size[1] * size[2] * sizeof(float));
  reader->read(orig, size, buf.get(), 0);
  // storage zero maps to 99.6137
  TEST_CHECK(std::abs(buf[0]+99.6137) < 0.001);
  reader->close();
}

template<typename T>
static void
do_readwriteT(const std::string& filename, SampleDataType dtype)
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  {
    std::shared_ptr<IZgyWriter> writer =
      IZgyWriter::open(ZgyWriterArgs()
                       .filename(filename)
                       .bricksize(16, 8, 32)
                       .size(16, 8, 160)
                       .datatype(dtype)
                       .datarange(-1,1));

    const T     fortytwo{42};
    const T     fifteen{15};
    const float pointnine{0.9f};
    const float pointeight{0.8f};

    writer->writeconst(size3i_t{0,0,0}, writer->bricksize(), &fortytwo);
    writer->write     (size3i_t{0,0,32}, size3i_t{1,1,1}, &fifteen);
    writer->write     (size3i_t{0,0,64}, size3i_t{1,1,1}, &pointnine);
    writer->writeconst(size3i_t{0,0,128}, writer->bricksize(), &pointeight);

    // Reading back from the still open writer.
    // There was no test coverage of these for int8 and int16.
    std::vector<T> wdata(160, -1);
    writer->read(size3i_t{0,0,0}, size3i_t{1,1,160}, wdata.data());
    auto wconst1 = writer->readconst(size3i_t{0,0,0},   writer->bricksize(), false);
    auto wconst2 = writer->readconst(size3i_t{0,0,32},  writer->bricksize(), false);
    auto wconst3 = writer->readconst(size3i_t{0,0,64},  writer->bricksize(), false);
    auto wconst4 = writer->readconst(size3i_t{0,0,96},  writer->bricksize(), false);
    auto wconst5 = writer->readconst(size3i_t{0,0,128}, writer->bricksize(), true);

    TEST_EQUAL((int)wdata[0],  42);
    TEST_EQUAL((int)wdata[32], 15);
    TEST_EQUAL((int)wdata[64], dtype==SampleDataType::int8 ? 114 : 29490);
    TEST_EQUAL((int)wdata[96], 0);
    TEST_EQUAL((int)wdata[128], dtype==SampleDataType::int8 ? 102 : 26214);
    TEST_CHECK(wconst1.first);
    TEST_EQUAL(wconst1.second, 42);
    TEST_CHECK(!wconst2.first);
    TEST_CHECK(!wconst3.first);
    TEST_CHECK(wconst4.first);
    TEST_EQUAL(wconst4.second, 0);
    TEST_CHECK(wconst5.first);
    TEST_EQUAL_FLOAT(wconst5.second, 0.8, 0.01);

    writer->close();
  }

  {
    // Close, open for read, try again to read the file.
    // These were already covered, but might as well do them here as well.
    std::shared_ptr<IZgyReader> reader = IZgyReader::open(filename);
    std::vector<T> rdata(160, -1);
    reader->read(size3i_t{0,0,0}, size3i_t{1,1,160}, rdata.data(), /*lod=*/0);
    auto rconst1 = reader->readconst(size3i_t{0,0,0},   reader->bricksize(), 0, false);
    auto rconst2 = reader->readconst(size3i_t{0,0,32},  reader->bricksize(), 0, false);
    auto rconst3 = reader->readconst(size3i_t{0,0,64},  reader->bricksize(), 0, false);
    auto rconst4 = reader->readconst(size3i_t{0,0,96},  reader->bricksize(), 0, false);
    auto rconst5 = reader->readconst(size3i_t{0,0,128}, reader->bricksize(), 0, true);
    TEST_EQUAL((int)rdata[0],  42);
    TEST_EQUAL((int)rdata[32], 15);
    TEST_EQUAL((int)rdata[64], dtype==SampleDataType::int8 ? 114 : 29490);
    TEST_EQUAL((int)rdata[96], 0);
    TEST_EQUAL((int)rdata[128], dtype==SampleDataType::int8 ? 102 : 26214);
    TEST_CHECK(rconst1.first);
    TEST_EQUAL(rconst1.second, 42);
    TEST_CHECK(!rconst2.first);
    TEST_CHECK(!rconst3.first);
    TEST_CHECK(rconst4.first);
    TEST_EQUAL(rconst4.second, 0);
    TEST_CHECK(rconst5.first);
    TEST_EQUAL_FLOAT(rconst5.second, 0.8, 0.01);

    reader->close();
  }
}

static void
test_readwrite_int8()
{
  LocalFileAutoDelete lad("readwrite_int8.zgy");
  do_readwriteT<std::int8_t>(lad.name(), SampleDataType::int8);
}

static void
test_readwrite_int16()
{
  LocalFileAutoDelete lad("readwrite_int16.zgy");
  do_readwriteT<std::int16_t>(lad.name(), SampleDataType::int16);
}

/**
 * The first brick of the (now confusingly named) Empty-v3
 * has been patched by hand to have this test pattern:
 *
 *   offset             v3 sample  v1 sample
 *    0*64 ..  0*64+9 :  0 .. 9     0 .. 9
 *    8*64 ..  8*64+9 : 10 .. 19   20 .. 29
 *   64*64 .. 64*64+9 : 20 .. 29   10 .. 10
 *  all others        : 0          0
 *
 * Test reading back this brick from both files.
 * the result that the application sees should be the same.
 */
static void
test_readsubtiles()
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v3.zgy"));
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size{64,64,64};
  std::unique_ptr<std::int8_t[]>buf(new std::int8_t[64*64*64]);
  memset(buf.get(), 0, 64*64*64);
  reader->read(orig, size, buf.get(), 0);
  TEST_EQUAL((int)buf[0], 0);
  TEST_EQUAL((int)buf[1], 1);
  TEST_EQUAL((int)buf[2], 2);
  TEST_EQUAL((int)buf[8*64+0], 10);
  TEST_EQUAL((int)buf[8*64+1], 11);
  TEST_EQUAL((int)buf[8*64+2], 12);
  TEST_EQUAL((int)buf[64*64+0], 20);
  TEST_EQUAL((int)buf[64*64+1], 21);
  TEST_EQUAL((int)buf[64*64+2], 22);
  reader->close();

  reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v1.zgy"));
  memset(buf.get(), 0, 64*64*64);
  reader->read(orig, size, buf.get(), 0);
  TEST_EQUAL((int)buf[0], 0);
  TEST_EQUAL((int)buf[1], 1);
  TEST_EQUAL((int)buf[2], 2);
  TEST_EQUAL((int)buf[8*64+0], 10);
  TEST_EQUAL((int)buf[8*64+1], 11);
  TEST_EQUAL((int)buf[8*64+2], 12);
  TEST_EQUAL((int)buf[64*64+0], 20);
  TEST_EQUAL((int)buf[64*64+1], 21);
  TEST_EQUAL((int)buf[64*64+2], 22);
  reader->close();
}

static void
test_readbadvt()
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v3.zgy"));
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size{16,16,16};
  std::unique_ptr<std::int16_t[]>buf(new std::int16_t[size[0] * size[1] * size[2]]);
  must_throw("storage cannot be null", [&](){
    reader->read(orig, size, (std::int8_t*)nullptr, 0);});
  must_throw("data type not supported", [&](){
    reader->read(orig, size, buf.get(), 0);});
  reader->close();
}

static void
test_readbadpos()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v3.zgy"));
  std::unique_ptr<std::int8_t[]>buf(new std::int8_t[1024]);
  // Negative size is caught when ZgyReader::read() wraps our raw buffer into
  // a DataBuffer, so the check in ZgyInternalBulk::readToExistingBuffer()
  // is not used.
  must_throw("size must be positive", [&](){
    reader->read(size3i_t{0,0,0}, size3i_t{1,-1,1}, buf.get(), 0);});
  must_throw("size must be positive", [&](){
    reader->read(size3i_t{0,0,0}, size3i_t{0,0,0}, buf.get(), 0);});
  must_throw("outside the valid range", [&](){
    reader->read(size3i_t{0,0,0}, size3i_t{1,1,1000}, buf.get(), 0);});
  must_throw("lod -1 is outside the valid range", [&](){
    reader->read(size3i_t{0,0,0}, size3i_t{1,1,1}, buf.get(), -1);});
  must_throw("lod 3 is outside the valid range", [&](){
    reader->read(size3i_t{0,0,0}, size3i_t{1,1,1}, buf.get(), 3);});
  // Repeat the tests using readconst.
  must_throw("region is empty", [&](){
    reader->readconst(size3i_t{0,0,0}, size3i_t{1,-1,1}, 0, true);});
  must_throw("region is empty", [&](){
    reader->readconst(size3i_t{0,0,0}, size3i_t{0,0,0}, 0, true);});
  must_throw("outside the valid range", [&](){
    reader->readconst(size3i_t{0,0,0}, size3i_t{1,1,1000}, 0, true);});
  must_throw("lod -1 is outside the valid range", [&](){
    reader->readconst(size3i_t{0,0,0}, size3i_t{1,1,1}, -1, true);});
  must_throw("lod 3 is outside the valid range", [&](){
    reader->readconst(size3i_t{0,0,0}, size3i_t{1,1,1}, 3, true);});
  reader->close();
}

static void
test_readnotopen()
{
  std::string fname = get_testdata("Fancy-int8.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  reader->close();

  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size{1,1,1};

  must_throw("not open for read", [&](){
    reader->readconst(orig, size, 0, true);
  });

  must_throw("not open for read", [&](){
    std::int8_t data{0};
    reader->read(orig, size, &data);
  });

#if 0 // Need an int16_t file for this.
  must_throw("not open for read", [&](){
    std::int16_t data{0};
    reader->read(orig, size, &data);
  });
#endif

  must_throw("not open for read", [&](){
    float data{0};
    reader->read(orig, size, &data);
  });
}

/**
 * The file's type will match the template argumnent.
 * Access will be attempted both using T and float.
 */
template<typename T>
static void
test_writenotopen()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("testfile.zgy");
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> next{0,0,64};
  const std::array<std::int64_t,3> size{1,1,1};
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename(lad.name())
    .size(33, 28, 92)
    .datatype(typeid(T)==typeid(std::int8_t) ? SampleDataType::int8 :
              typeid(T)==typeid(std::int16_t) ? SampleDataType::int16 :
              SampleDataType::float32)
    .datarange(-100,+100);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  const T     tdata[6]{42,43,44,45,46,47};
  const float fdata[6]{42,43,44,45,46,47};
  writer->writeconst(size3i_t{0,0,0}, size3i_t{1,1,1}, &tdata[0]);
  writer->writeconst(size3i_t{0,0,1}, size3i_t{1,1,1}, &fdata[1]);
  writer->write     (size3i_t{0,0,2}, size3i_t{1,1,1}, &tdata[2]);
  writer->write     (size3i_t{0,0,3}, size3i_t{1,1,1}, &fdata[3]);
  writer->writeconst(size3i_t{0,0,4}, size3i_t{1,1,1}, &tdata[4]);
  writer->writeconst(size3i_t{0,0,5}, size3i_t{1,1,1}, &fdata[5]);
  writer->finalize();
  writer->close();

  // Overwriting
  must_throw("not open for write", [&](){
    const T data{0};
    writer->write(orig, size, &data);
  });
  must_throw("not open for write", [&](){
    const float data{0};
    writer->write(orig, size, &data);
  });
  must_throw("not open for write", [&](){
    const T data{0};
    writer->writeconst(orig, size, &data);
  });
  must_throw("not open for write", [&](){
    const float data{0};
    writer->writeconst(orig, size, &data);
  });

  // Writing new
  must_throw("not open for write", [&](){
    const T data{-1};
    writer->write(next, size, &data);
  });
  // Now that the error is reported at a high level
  // it will not be treated as an I/O error and will
  // not set the error flag. Which is correct.
  // Fake an I/O error.
  writer->set_errorflag(true);
  must_throw("previous errors", [&](){
    const float data{-1};
    writer->write(next, size, &data);
  });
  writer->set_errorflag(false);
  must_throw("not open for write", [&](){
    const float data{-1};
    writer->write(next, size, &data);
  });
  must_throw("not open for write", [&](){
    const T data{-1};
    writer->writeconst(next, size, &data);
  });
  must_throw("not open for write", [&](){
    const float data{-1};
    writer->writeconst(next, size, &data);
  });

  // Finish
  must_throw("not open for write", [&](){
    writer->finalize();
  });
  // A double close is acceptable.
  //must_throw("not open for write", [&](){
  writer->close();
  //});

  // This is a mostly unrelated test.
  // Check that writeconst() uses a read-modify-write when needed
  // and doesn't just clobber the entire brick when the caller
  // didn't provide a brick aligned region.
  std::shared_ptr<IZgyReader> reader = IZgyReader::open(lad.name());
  float tcheck[6]{0};
  float fcheck[6]{0};
  reader->read(size3i_t{0,0,0}, size3i_t{1,1,6}, &fcheck[0], 0);
  reader->read(size3i_t{0,0,0}, size3i_t{1,1,6}, &tcheck[0], 0);
  reader->close();
  // Verify using the same valuetype as written:
  // storage for even numbers, float for odd.
  if (verbose())
    std::cout << "Read back: "
              << tcheck[0] << " "<< fcheck[1] << " "
              << tcheck[2] << " "<< fcheck[3] << " "
              << tcheck[4]<< " " << fcheck[5]
              << std::endl;
}

static void
test_createargs()
{
  LocalFileAutoDelete lad("createargs.zgy");
  const std::string filename = lad.name();
  lad.disarm(); // Don't expect to actually create it.
  must_throw("cannot be empty or negative", [&lad](){
    IZgyWriter::open(ZgyWriterArgs().filename(lad.name()).size(0,1,1));
    lad.arm();
  });
  must_throw("too large", [&lad](){
    IZgyWriter::open(ZgyWriterArgs().filename(lad.name()).size(1,1,1LL<<31));
    lad.arm();
  });
}

static void
test_ioerror()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("ioerror.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(128, 128, 128)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  writer->set_errorflag(true);
  // TODO-Test: It would be better to provoke a real I/O error here,
  // but this might not be possible without breaking encapsulation.
  TEST_CHECK(writer->errorflag());
  must_throw("previous error", [writer](){
    const float fortytwo{42};
    writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  });
  TEST_CHECK(writer->errorflag());
  // finalize() and close() are now allowed even after en error,
  // as they will try to salvage what they can and clear lowres etc.
  //must_throw("previous error", [writer](){
    writer->finalize();
  //});
  //must_throw("previous error", [writer](){
    writer->close();
  //});
}

static void
do_test_finalize(int mode)
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("incomplete_" + std::to_string(mode) + ".zgy");
  ZgyWriterArgs args = OpenZGY::ZgyWriterArgs()
    .size(128, 128, 128)
    .filename(lad.name());
  if (mode == 7)
    args.size(7, 13, 19);
  const std::int64_t samplecount{(mode==7) ? 7*13*19 : 128*128*128};
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  const float fortyone{41};
  const float fortytwo{42};
  switch (mode) {
  case 1: // write, implicit finalize, close.
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->close();
    break;
  case 2: // write, finalize, close.
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->finalize();
    writer->close();
    break;
  case 3: // write, finalize, write, finalize again, close.
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->finalize();
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->finalize();
    writer->close();
    break;
  case 4: // write, finalize, close_incomplete should then be same as close.
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->finalize();
    writer->close_incomplete();
    break;
  case 5: // write, close without finalize
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->close_incomplete();
    break;
  case 6: // write, finalize, write, close without needed second finalize
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->finalize();
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->close_incomplete();
    break;
  case 7: // as 6 but with a tiny file that had no lowres anyway.
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortyone);
    writer->finalize();
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    writer->close_incomplete();
    break;
  default:
    TEST_CHECK_(false, "unrecognized testcase %d", mode);
    break;
  }
  writer.reset();
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(lad.name());
  SampleStatistics stat = reader->statistics();
  SampleHistogram hist = reader->histogram();
  float first{0};
  reader->read(size3i_t{0,0,0}, size3i_t{1,1,1}, &first, 0);
  // Testing the #7 corner case in particular.
  TEST_CHECK_(first == 42, "first expect 42 got %g", first);
  switch (mode) {
  case 1:
  case 2:
  case 3:
  case 4: // Complete file ought to have been created.
    {
      TEST_EQUAL(reader->nlods(), 2);
      std::pair<bool,double> constval =
        reader->readconst(size3i_t{0,0,0}, size3i_t{64,64,64}, 1, true);
      TEST_CHECK(constval.first);
      TEST_EQUAL(constval.second, 42);
      TEST_EQUAL(stat.cnt, samplecount);
      TEST_EQUAL(stat.sum, samplecount*42.0);
      TEST_EQUAL(stat.ssq, samplecount*42.0*42.0);
      TEST_EQUAL(stat.min, 42.0);
      TEST_EQUAL(stat.max, 42.0);
      TEST_EQUAL(hist.samplecount, samplecount);
      //FAILS! TEST_CHECK_(hist.minvalue == 42, "minvalue expect 42 got %lg", hist.minvalue);
      // SINGLEPASS: The range can end up wider due to dynamic histogram.
      TEST_CHECK_(hist.maxvalue >= 42, "maxvalue expect >= 42 got %lg", hist.maxvalue);
    }
    break;
  case 5:
  case 6:
  case 7: // Lowres data, histogram, statistics not available.
    {
      // SINGLEPASS: Lowres is now available, close_incomplete to be retired.
      // However, case 7 (single brick) has no lod 1 in the first place.
      if (mode == 7) {
        TEST_EQUAL(reader->nlods(), 1);
        must_throw("outside the valid range", [&](){
          reader->readconst(size3i_t{0,0,0}, size3i_t{2,2,2}, 1, true);
        });
      }
      else {
        TEST_EQUAL(reader->nlods(), 2/*was 1*/);
        reader->readconst(size3i_t{0,0,0}, size3i_t{2,2,2}, 1, true);
      }
      // SINGLEPASS: Now thay *are* available.
      TEST_EQUAL(stat.cnt, /*was 0*/ samplecount);
      TEST_EQUAL_FLOAT(stat.sum, /*was 0*/ samplecount*42.0, 0.01);
      TEST_EQUAL_FLOAT(stat.ssq, /*was 0*/ samplecount*42.0*42.0, 10.0);
      TEST_CHECK(stat.min <= stat.max);
      TEST_EQUAL(hist.samplecount, stat.cnt);
      // SINGLEPASS: The range can end up wider due to dynamic histogram.
      TEST_CHECK_(hist.minvalue <= 42, "minvalue expect <= 42 got %lg", hist.minvalue);
      TEST_CHECK_(hist.maxvalue >= 42, "maxvalue expect >= 42 got %lg", hist.maxvalue);
    }
    break;
  default:
    TEST_CHECK_(false, "unrecognized testcase %d", mode);
    break;
  }
  reader->close();
  reader.reset();
}

template<int MODE>
static void
test_finalize()
{
  do_test_finalize(MODE);
}

/*
 * Testing the behavior of constant-value bricks and low resolution computation.
 *
 * Consider this survey seen from above:
 *      0    32    64    96    128
 *  0   +--+--+--+--+--+--+--+--+-+
 *      |  |  |  |  |  |  |  |  |C|
 * 16   +--+--+--+--+--+--+--+--+-+
 *      |  |  |  |  |  |  |  |  | |
 * 32   +--+--+--+--+--+--+--+--+-+
 *      |  |  |AA|AA|BB|  |  |  |D|      A/B: (32-3,32-2)..(80+1,64+4)
 * 48   +--+--+--+--+--+--+--+--+-+      C:   (128-2,0)..(136,16+4)
 *      |  |  |AA|AA|BB|  |  |  |D|      D:   (128,32)..(136,64)
 * 64   +--+--+--+--+--+--+--+--+-+
 *      |  |  |  |  |  |  |  |  | |
 * 80   +--+--+--+--+--+--+--+--+-+
 *
 * bricksize  = (16,16,16) * char
 * surveysize = (136,80,28)
 * surveysize in bricks: 8.5,5,1.75
 *
 * Bricks A,B,C,D are all explicitly set to all const, with a value
 * that is not the same as the default value for missing bricks.
 * Actually the region of all-const will be slightly larger than
 * shown.
 *
 * Bricks not labeled are all set to some generated pattern.
 *
 * Case A sees 4 brick-columns with all dead traces. The lod1
 * brick-column will all have dead traces. The decimation algorithm
 * should not be called at all. For performance reasons it is not
 * acceptable to first create regular bricks filled with a single
 * value, and later have the code in the writer turn that into a
 * scalar brick.
 *
 * CAVEAT: This test cannot actually verify the behavior above,
 * because if the brick was inflated then it would get deflated again
 * by write(). After having caused a performance problem that can't be
 * seen in this tiny test. Check it manually by examining the log
 * output and/or temporarily disable calls to isAllSame() in
 * ZgyInternalBulk::readToNewBuffer[s] and GenLodImpl::_calculate().
 * I think that will prevent the automatic deflation.
 *
 * Case B sees only some of the 4 input brick-columns being dead. For
 * performance reasons it is not acceptable to run decimation on the
 * constant value bricks.
 *
 * Also in this case the automated test will not be able to spot
 * problems, because it cannot see whether half of the samples were
 * written by a flood fill instead of running the decimation.
 *
 * Case C sees 4 brick columns where two are outside the survey and
 * one is dead. There will be real lowres data from the remaining
 * column. No automatic verification.
 *
 * Case D sees 4 brick columns where two are outside the survey and
 * two are dead. The low resolution brick should be all-const with the
 * same scalar value as D. Take care that the "D" bricks and the
 * "outside" bricks don't end up with different scalars and that this
 * leads to the brick above becoming non-const. This is critical for
 * cloud access where we can only write each brick once.
 *
 * Finally, here is something the automated test can actually check.
 * If the lod1 brick for area D isn't constant then it probably
 * ended up as a mix of our novalue (42) and the system's (0).
 *
 * It is also possible to verify the number of normal vs. scalar
 * bricks on the file. Lod0 has 9 * 5 brick columns total, 9 of them
 * scalar. So, 36 normal. lod1 has 5 * 3 brick columns total, 2 scalar
 * (A and D) and 13 normal. lod2 has 3 * 2 brick columns, all normal.
 * lod3 has 2 * 1 normal brick columns and the final lod4 has 1
 * normal. Grand total 36+13+6+2+1=58 normal, 9+2=11 scalar. Each
 * brick column in lod0 has 2 vertical bricks. Bringing the brick count up to
 * 2*36+13+6+2+1=94 normal, 2*9+2=20 scalar.
 */
static void
test_genlod()
{
  LocalFileAutoDelete lad("testgenlod.zgy");
  std::vector<std::int8_t> dead(28, 42);
  std::vector<std::int8_t> live(28, 0);
  int value{0};
  for (auto& it : live)
    it = static_cast<std::int8_t>(++value);
  const std::array<std::int64_t,3> size{136,80,28};
  const std::array<std::int64_t,3> bricksize{16,16,16};
  const std::array<std::int64_t,3> zero3d{0,0,0};
  std::vector<std::int8_t> survey(size[0] * size[1] * size[2]);
  for (std::size_t ii = 0; ii < survey.size(); ii += size[2])
    std::copy(live.data(), live.data() + size[2], &survey.data()[ii]);
  // Region A and B
  for (std::int64_t ii=29; ii<81; ++ii)
    for (std::int64_t jj=30; jj<68; ++jj)
      std::copy(dead.data(), dead.data() + size[2],
                &survey.data()[ii*size[1]*size[2] + jj*size[2]]);
  // Region C
  for (std::int64_t ii=126; ii<136; ++ii)
    for (std::int64_t jj=0; jj<20; ++jj)
      std::copy(dead.data(), dead.data() + size[2],
                &survey.data()[ii*size[1]*size[2] + jj*size[2]]);
  // Region D
  for (std::int64_t ii=128; ii<136; ++ii)
    for (std::int64_t jj=32; jj<64; ++jj)
      std::copy(dead.data(), dead.data() + size[2],
                &survey.data()[ii*size[1]*size[2] + jj*size[2]]);
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename(lad.name())
    .datatype(SampleDataType::int8)
    .datarange(-128, +127)
    .bricksize(bricksize[0], bricksize[1], bricksize[2])
    .size(size[0], size[1], size[2]);
  std::shared_ptr<OpenZGY::IZgyWriter> writer =
    OpenZGY::IZgyWriter::open(args);
  const std::int8_t fortytwo{42};
  writer->writeconst(zero3d, size, &fortytwo);
  writer->write(zero3d, size, survey.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Average}, nullptr);
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(lad.name());
  for (std::int64_t ii=0; ii<(size[0]+1)/2; ii += bricksize[0]) {
    for (std::int64_t jj=0; jj<(size[1]+1)/2; jj += bricksize[1]) {
      std::pair<bool,double> c = reader->readconst
        (std::array<std::int64_t,3>{ii, jj, 0}, reader->bricksize(), 1, false);
      if (verbose()) {
        std::cout << std::boolalpha
                  << "(" << ii/bricksize[0] << ", " << jj/bricksize[1] << ")"
                  << " -> " << c.first << " " << c.second << std::endl
                  << std::noboolalpha;
      }
      if (ii/bricksize[0] == 1 && jj/bricksize[1] == 1) {
        TEST_CHECK(c.first);
      }
      else if (ii/bricksize[0] == 4 && jj/bricksize[1] == 1) {
        TEST_CHECK(c.first);
      }
      else  {
        TEST_CHECK(!c.first);
      }
    }
  }
  std::shared_ptr<const FileStatistics> stats = reader->filestats();
  TEST_EQUAL(stats->brickNormalCount(), 94);
  TEST_EQUAL(stats->brickCompressedCount(), 0);
  TEST_EQUAL(stats->brickMissingCount(), 0);
  TEST_EQUAL(stats->brickConstantCount(), 20);
}

/**
 * Write a huge survey consisting almost exclusively of empty bricks.
 * On finalize, a shortcut should ensure that the genlod algorithm
 * isn't run for all bricks. If it is, the test will take a very
 * long time.
 *
 * See also test_ambig2 which also creates a huge file.
 */
static void
test_genlod2()
{
  // The test depends on measuring elapsed time, so valgrind is out.
  if (is_running_on_valgrind()) {
    if (verbose())
      std::cout << "\nSkipping api.genlod2 when running under valgrind.\n";
    return;
  }
  InternalZGY::Timer timer(true, "test_genlod2");

  LocalFileAutoDelete lad("genlod2.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(131313,131,1313) // 2052 x 3 x 21 = 129276 bricks ~= 63 GB
    .datatype(SampleDataType::int16)
    .datarange(-32768, +32767)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  const std::array<std::int64_t,3> zero{0,0,0};
  const std::int16_t fillvalue{1962};
  const std::int16_t fortytwo{42};
  writer->writeconst(zero, writer->size(), &fillvalue);
  writer->writeconst(zero, IZgyWriter::size3i_t{1,1,32}, &fortytwo);
  writer->finalize();
  writer->close();
  timer.stop();
  // A single test run on a powerful machine took 0.7 seconds with the
  // old code, and 66 seconds with buggy code that did not have the
  // shortcut. 6 seconds should be a safe value to test for.
  // Not on Windows, though.
  // Not on GitLab builds either. Maybe just drop it.
  //#ifndef _WIN32
  //TEST_CHECK(timer.getTotal() <= 6);
  //#endif
}

#ifdef HAVE_SD

static bool
read_first_sample(const std::shared_ptr<OpenZGY::IZgyReader>& reader,
                  bool expect_ok,
                  const std::string& message)
{
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size{1,1,1};
  float bulk[1]{0};
  bool ok = false;

  if (verbose())
    std::cout << message << " is being tested." << std::endl;
  try {
    reader->read(orig, size, &bulk[0], 0);
    if (!expect_ok)
      std::cout << message << " ERROR, should not have worked."
                << std::endl;
    else if (std::abs(bulk[0]-0.0039) > 0.001)
      std::cout << message << " FAILED with wrong sample value "
                << bulk[0] << std::endl;
    else
      ok = true;
  }
  catch (const std::exception& ex) {
    if (expect_ok) {
      std::cout << message << " FAILED with: " << ex.what() << std::endl;
      ok = false;
    }
    else {
      ok = true;
      if (verbose())
        std::cout << "Got expected: " << ex.what() << std::endl;
    }
  }
  return ok;
}

static std::string
cloud_synt2_name()
{
  std::string testfolder = InternalZGY::Environment::getStringEnv("OPENZGY_SDTESTDATA", "sd://opendes/slb-testdata");
  if (testfolder.back() != '/')
    testfolder += "/";
  return testfolder + "Synt2.zgy";
}

/**
 * \brief Use the SDAPI token callback mechanism.
 *
 * Normal tokens expire in 60 minutes, and storage tokens can add
 * another 60 minutes to that number. If an application waits 60+
 * minutes with no activity then both tokens should have timed out.
 * If the application reads data ever so often then 120 minutes
 * might be needed.
 *
 * To properly run this test the app should be started with credentials
 * that allow refresh (e.g. a client credentials grant) and should be
 * allowed to run for 65+ minutes, setting $OPENZGY_TEST_TOKENCB_SLEEP
 * to the desired time. Since the test needs to be run manually anyway
 * the result should also be verified manually. With .verbose enabled
 * there should be exactly two redacted tokens logged for the first
 * reader, and they should be different. The second reader might see
 * more requests and all the results should be the same.
 *
 * If the credentials are not refreshable then the two readers should
 * behave the same. And if the sleep is less than 60 minutes the
 * credentials might or might not time out.
 *
 * Another thing to check is how the second reader behaves after the
 * credentials has timed out. Specifically that it doesn't keep
 * invoking the callback in an infinite loop. The system will probably
 * do a fixed number of retries with exponential backoff.
 */
static void
test_tokencb2()
{
  using InternalZGY::Environment;
  if (verbose())
    std::cout << std::endl;

  static auto redact = [](const std::string& s) {
                         return s.size() < 20 ? s :
                           s.substr(0,7) + "..." + s.substr(s.size() - 6);
                       };

  std::function<std::string()> cb = Test_Utils::get_token_callback();
  std::string initial = cb(); // for the non-refreshing token.
  std::function<std::string()> refreshed_cb =
    [cb]() {
      std::string token = cb();
      if (verbose())
          std::cout << "refreshed token: " << redact(token) << std::endl;
      return token;
    };
  std::function<std::string()> norefresh_cb =
    [initial]() {
      std::string token = initial;
      if (verbose())
          std::cout << "norefresh token: " << redact(token) << std::endl;
      return token;
    };
  auto refreshed_context = SeismicStoreIOContext()
    .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
    .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
    .sdtokencb(refreshed_cb);
  auto norefresh_context = SeismicStoreIOContext()
    .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
    .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
    .sdtokencb(norefresh_cb);

  std::shared_ptr<IZgyReader> refreshed_reader =
    IZgyReader::open(cloud_synt2_name(), &refreshed_context);
  std::shared_ptr<IZgyReader> norefresh_reader =
    IZgyReader::open(cloud_synt2_name(), &norefresh_context);

  TEST_CHECK(read_first_sample(refreshed_reader, true, "Normal read #1"));
  TEST_CHECK(read_first_sample(norefresh_reader, true, "Normal read #2"));

  int sleeptime = Environment::getNumericEnv("OPENZGY_TEST_TOKENCB_SLEEP", -1);
  if (sleeptime >= 0) {
    std::cout << "Sleeping " << sleeptime
              << " minutes waiting for tokens to expire."
              << std::endl;
    std::this_thread::sleep_for(std::chrono::minutes(sleeptime));
    std::cout << "Awake!" << std::endl;
  }

  // Note the ordering. Storage tokens are shared among managers,
  // so if the "refreshed" read is done first it will make a valid
  // storage token that the next read can borrow.
  TEST_CHECK(read_first_sample(norefresh_reader, true, "Second read norefresh"));
  TEST_CHECK(read_first_sample(refreshed_reader, true, "Second read refreshed"));

  try {
    refreshed_reader->close();
    norefresh_reader->close();
  }
  catch (const std::exception& ex)
  {
    std::cout << "Exception closing the datasets: " << ex.what() << std::endl;
    TEST_CHECK(false && "failed in close");
  }
}

#endif

static void
test_ZgyWriterArgs()
{
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename("testfile")
    .size(33, 28, 92)
    .datatype(SampleDataType::int16)
    .zunit(UnitDimension::time, "ms", 1000);
  std::stringstream ss;
  args.dump(ss);
  if (verbose())
    std::cout << ss.str() << std::flush;
  TEST_CHECK(ss.str().find("filename:    \"testfile\"") != std::string::npos);
  TEST_CHECK(ss.str().find("size:        (33,28,92)") != std::string::npos);
  TEST_CHECK(ss.str().find("bricksize:   (64,64,64)") != std::string::npos);
  TEST_CHECK(ss.str().find("\"ms\"") != std::string::npos);
}

static void
do_write_once(const std::string& filename, const IOContext *context = nullptr)
{
  ZgyWriterArgs args = ZgyWriterArgs()
    .iocontext(context)
    .filename(filename)
    .size(33, 28, 92)
    .datatype(SampleDataType::int16)
    .datarange(-32768,+32767)
    .ilstart(1).ilinc(2)
    .xlstart(500).xlinc(5)
    .zstart(100).zinc(4)
    .hunit(UnitDimension::length, "m", 1)
    .zunit(UnitDimension::time, "ms", 1000)
    .corners(ZgyWriterArgs::corners_t{{{5,7},{5,107},{205,7},{205,107}}});
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  std::vector<float> data(2*3*4, -1000);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t bsize{64,64,64};
  const OpenZGY::IZgyWriter::size3i_t count{2,3,4};
  float fortytwo{42};
  writer->writeconst(origin, bsize, &fortytwo);
  writer->write(origin, count, data.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{}, nullptr);
  writer->close();
}

static void
do_check_written(const std::string& filename, const IOContext* context = nullptr)
{
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(filename, context);
  if (verbose()) {
    std::cout << "\n";
    dump_api(reader, std::cout);
    //reader->dump(std::cout);
    reader->filestats()->dump(std::cout, "filestats: ");
  }
  std::unique_ptr<float[]> checkdata(new float[64*64*64]);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t bsize{64,64,64};
  reader->read(origin, bsize, checkdata.get(), 0);
  TEST_EQUAL(checkdata[0], -1000);
  TEST_EQUAL(checkdata[63], 42);
  const SampleHistogram h = reader->histogram();
  TEST_EQUAL_FLOAT(h.minvalue, -32768, 1e-5);
  TEST_EQUAL_FLOAT(h.maxvalue, +32767, 1e-5);
  TEST_EQUAL(h.samplecount, 33*28*92);
  TEST_EQUAL(h.bins[124], 2*3*4);
  TEST_EQUAL(h.bins[128], 33*28*92 - 2*3*4);
  const ZgyWriterArgs::corners_t corners = reader->corners();
  const double eps = 1.0e-10;
  TEST_EQUAL_FLOAT(corners[0][0],   5, eps);
  TEST_EQUAL_FLOAT(corners[0][1],   7, eps);
  TEST_EQUAL_FLOAT(corners[1][0],   5, eps);
  TEST_EQUAL_FLOAT(corners[1][1], 107, eps);
  TEST_EQUAL_FLOAT(corners[2][0], 205, eps);
  TEST_EQUAL_FLOAT(corners[2][1],   7, eps);
  TEST_EQUAL_FLOAT(corners[3][0], 205, eps);
  TEST_EQUAL_FLOAT(corners[3][1], 107, eps);
  TEST_CHECK(reader->size() == (OpenZGY::IZgyWriter::size3i_t{33, 28, 92}));
  TEST_CHECK(reader->datatype() == SampleDataType::int16);
  //TEST_CHECK(reader->hunitname() == "m");
  //TEST_CHECK(reader->zunitname() == "ms");
  TEST_EQUAL(reader->annotstart()[0], 1);
  TEST_EQUAL(reader->annotstart()[1], 500);
  TEST_EQUAL(reader->zstart(), 100);
  reader->close();
}

static void
test_write()
{
  LocalFileAutoDelete lad("testfile.zgy");
  do_write_once(lad.name());
  do_check_written(lad.name());
}

/**
 * Writing with a dummy compressor means all bricks end up uncompressed
 * but there are still a few minor differences when a compressor is given:
 *  \li Writing the same block twice, including Read/modify/write is not allowed.
 *  \li Version cannot be 3 for compressed and cannot be 4 for uncompressed.
 *  \li Compressed files might (but currently do not) skip the header padding.
 *      See comments in ZgyInternalMeta::flushMeta().
 */
static std::shared_ptr<const FileStatistics>
do_test_copy_slurp_8(const std::string& iname, const std::string& oname, const std::string& compressor)
{
  // TODO-Low: Consolidate with test_copy.

  std::vector<std::string> compress_args;
  if (compressor == "ZFP") {
    compress_args.push_back("snr");
    compress_args.push_back("30");
  }

  // read
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(iname);
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size = reader->size();
  std::unique_ptr<float[]>buf(new float[size[0] * size[1] * size[2]]);
  reader->read(orig, size, buf.get(), 0);

  // write
  ZgyWriterArgs args = ZgyWriterArgs()
    .metafrom(reader)
    .datatype(SampleDataType::float32) // Also for uncompressed in this case.
    .filename(oname);
  if (!compressor.empty()) {
    args
      .datatype(SampleDataType::float32)
      .compressor(compressor, compress_args);
    // Note that lodcompressor defaults to compressor if not set.
  }

  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  writer->write(origin, size, buf.get());
  if (compressor.empty()) {
    writer->write(origin, size, buf.get());
  }
  else {
    // Compressed files do not allow updates.
    must_throw("compressed data in update mode", [&](){
      writer->write(origin, size, buf.get());
    });
  }
  writer->finalize();
  reader->close();
  writer->close();

  // read back the file statistics only
  std::shared_ptr<OpenZGY::IZgyReader> copy = OpenZGY::IZgyReader::open(oname);
  std::shared_ptr<const FileStatistics> filestats = copy->filestats();
  copy->close();
  if (verbose())
    filestats->dump(std::cout);
  return filestats;
}

static void
test_compress_noop()
{
  const std::string filename = get_testdata("Empty-v3.zgy");
  LocalFileAutoDelete lad("tempcopy.zgy");
  std::shared_ptr<const FileStatistics> stats =
    do_test_copy_slurp_8(filename, lad.name(), "Null");

  TEST_EQUAL(stats->fileVersion(), 3); // v4 only set if actual compressed
  TEST_CHECK(stats->isCompressed() == false); // looked for compressed bricks.
  TEST_EQUAL(stats->fileSize(), 4*64*64*64*4); // header plus 3 normal bricks
  TEST_EQUAL(stats->brickNormalCount(), 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_EQUAL(stats->brickCompressedCount(), 0);
  TEST_EQUAL(stats->brickMissingCount(), 0);
  TEST_EQUAL(stats->brickConstantCount(), 42);
  TEST_CHECK(stats->compressionFactor() > 0.99);
}

static void
test_compress_zfp()
{
  const std::string filename = get_testdata("Empty-v3.zgy");
  LocalFileAutoDelete lad("tempcopy.zgy");
  std::shared_ptr<const FileStatistics> stats =
    do_test_copy_slurp_8(filename, lad.name(), "ZFP");

  TEST_CHECK(stats->fileVersion() != 3);
  TEST_CHECK(stats->isCompressed() == true);
  TEST_CHECK(stats->fileSize() < 3*64*64*64*4); // header plus 3 comp. bricks
  TEST_CHECK(stats->fileSize() > 1*64*64*64*4); // still a full header
  TEST_EQUAL(stats->brickNormalCount(), 0);
  TEST_EQUAL(stats->brickCompressedCount(), 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_EQUAL(stats->brickMissingCount(), 0);
  TEST_EQUAL(stats->brickConstantCount(), 42);
  TEST_CHECK(stats->compressionFactor() < 0.9);
  TEST_CHECK(stats->brickCompressedSize() < 2*64*64*64);
}

/**
 * This is here just for completeness, there is near 100% overlap with
 * other tests made in this file.
 */
static void
test_compress_off()
{
  const std::string filename = get_testdata("Empty-v3.zgy");
  LocalFileAutoDelete lad("tempcopy.zgy");
  std::shared_ptr<const FileStatistics> stats =
    do_test_copy_slurp_8(filename, lad.name(), "");

  TEST_CHECK(stats->fileVersion() != 4);
  TEST_CHECK(stats->isCompressed() == false); // looked for compressed bricks.
  TEST_EQUAL(stats->fileSize(), 4*64*64*64*4); // header plus 3 normal bricks
  TEST_EQUAL(stats->brickNormalCount(), 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_EQUAL(stats->brickCompressedCount(), 0);
  TEST_EQUAL(stats->brickMissingCount(), 0);
  TEST_EQUAL(stats->brickConstantCount(), 42);
  TEST_CHECK(stats->compressionFactor() > 0.99);
}

/**
 * If a file is written with real data but is then cleared before closing
 * it for the first time, then the histogram might end up permanently
 * useless. See ZgyInternalBulk::newTrackedBricksTryEnable().
 *
 * TODO-WIP-BrickedAPI: Document these "features" to and users?
 */
static void
test_inflated_constant()
{
  typedef std::array<std::int64_t,3> size3i_t;
  Test_Utils::LocalFileAutoDelete lad("api_inflated_constant.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .filename(lad.name())
     .size(64, 64, 1024));
  const float sixteen{16};
  const float fortytwo{42};
  const float ninetynine{99};
  writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
  writer->writeconst(size3i_t{0,0,0}, size3i_t{1,1,1}, &sixteen);
  writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
  writer->close();

  // The first write should ensure that the old (zero) defaultvalue
  // is not included in the statistics. The histogram might contain zero
  // anyway? The second write allocates one brick, and the third one
  // resets everything but might not deallocate / leak the brick.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .filename(lad.name()));
  writer->writeconst(size3i_t{0,0,64}, size3i_t{1,1,1}, &ninetynine);
  writer->close();
  writer.reset();

  std::shared_ptr<OpenZGY::IZgyReader> reader;
  reader = OpenZGY::IZgyReader::open(lad.name());

  SampleStatistics stat = reader->statistics();
  SampleHistogram hist = reader->histogram();
  std::shared_ptr<const FileStatistics> filestats = reader->filestats();
  reader->close();
  reader.reset();

  if (verbose())
    filestats->dump(std::cout, "   ");

  // Will "16" still be remembered as part of the statistics?
  // In the old code this depends on whether trimRange() was
  // called from HistogramBuilder::operator+= or not. (?)
  // In the new code, as long as we did the final clear-survey
  // inside the first session then it would have been recognized
  // as such and caused the statistics to be reset.
  TEST_EQUAL(stat.min, 42);
  TEST_EQUAL(stat.max, 99);
  TEST_EQUAL(stat.cnt, 64*64*1024);

  // Will zero or 16 be included in the histogram range?
  // In the old code, 16 gets included because the limits
  // are tracked during the first pass. And zero apparently
  // gets included as well. I already documented that
  // feature as unspecified.
  // In the new code they shouldn't, since the histogram
  // was completely cleared by the third write and should
  // have caused a ~42..42 range. Due to how the dynamic
  // histogram works, the range will actually be somewhat
  // wider. But not as wide as to include zero or 16.
  //
  // SINGLEPASS:
  // Will the new "99" value added in the second session
  // be included? In the old code, yes, because the histogram
  // is rebuilt. In the new code, no, because the file isn't
  // recognized as all-const (see newTrackedBricksTryEnable)
  // which means the histogram range cannot change.
  // Unless the leaking bricks feature is on. Probably.
  //
  TEST_CHECK(hist.minvalue >= 20);
  TEST_CHECK(hist.maxvalue >= 99);
  TEST_EQUAL(hist.samplecount, 64*64*1024);

  // There will be 4 real lowres bricks and two fullres.
  // Or just one fullres if the space for the first brick
  // got deliberately leaked in writeAlignedBrickList().

  // SINGLEPASS: Currently, leaking is ON. So, 5 real bricks.
  TEST_EQUAL(filestats->brickNormalCount(),     6-1);
  TEST_EQUAL(filestats->brickCompressedCount(), 0);
  TEST_EQUAL(filestats->brickMissingCount(),    0);
  TEST_EQUAL(filestats->brickConstantCount(),  25+1);
}

/**
 * In an int8 / int16 file, verify that statistics and histogram
 * both reflect the values after clipping to the coding range.
 */
static void
test_stats_after_clip()
{
  typedef std::array<std::int64_t,3> size3i_t;
  Test_Utils::LocalFileAutoDelete lad("api_stats_after_clip.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .filename(lad.name())
     .datatype(SampleDataType::int8)
     .datarange(-1280, +1270)
     .size(64, 64, 100));
  const float data[]{160, -50, 420, 9999};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,4}, &data[0]);
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader;
  reader = OpenZGY::IZgyReader::open(lad.name());

  SampleStatistics stat = reader->statistics();
  SampleHistogram hist = reader->histogram();
  reader->close();
  reader.reset();

  TEST_EQUAL(stat.min, -50);
  TEST_EQUAL(stat.max, 1270); // 9999 was clipped to this.
  TEST_EQUAL(stat.cnt, 64*64*100);

  TEST_EQUAL(hist.minvalue, -1280);
  TEST_EQUAL(hist.maxvalue, +1270);
  TEST_EQUAL(hist.samplecount, 64*64*100);
  TEST_EQUAL(hist.bins.back(), 1); // Bucket for +1270
}

#ifdef HAVE_SD
static void
test_write_cloud()
{
  Test_Utils::CloudFileAutoDelete cad("writecloud.zgy", Test_Utils::default_sd_context());
  do_write_once(cad.name(), Test_Utils::default_sd_context());
  do_check_written(cad.name(), Test_Utils::default_sd_context());
}

static void
test_write_cloud_mt()
{
  // Note: I cannot use do_write_once() and do_check_written() because
  // I need a file with a larger number of normal bricks.
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  SeismicStoreIOContext not_default_context(*Test_Utils::default_sd_context());
  not_default_context.segsize(7); // In MB, so 12 bricks fit in one segment
  not_default_context.segsplit(3); // Write segment 1&2&3 in parallel, then...
  Test_Utils::CloudFileAutoDelete cad("writecloud_mt.zgy", Test_Utils::default_sd_context());

  // write
  const std::string filename = cad.name();
  ZgyWriterArgs args = ZgyWriterArgs()
    .iocontext(&not_default_context)
    .filename(filename)
    .size(128, 128, 1024)
    .datatype(SampleDataType::float32);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  const size3i_t size = writer->size();
  const size3i_t bs = writer->bricksize();
  float data = 1;
  for (std::int64_t ii=0; ii<size[0]; ii+=bs[0]) {
    for (std::int64_t jj=0; jj<size[1]; jj+=bs[1]) {
      for (std::int64_t kk=0; kk<size[2]; kk+=bs[2]) {
        writer->write(size3i_t{ii,jj,kk}, size3i_t{1,1,1}, &data);
        data += 1;
      }
    }
  }
  writer->finalize();
  writer->close();
  writer.reset();

  // read back and check
  not_default_context.cputhreads(3); // Let OpenZGY do its own MT also on read
  not_default_context.iothreads(3);
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(filename, &not_default_context);
  data = 1;
  for (std::int64_t ii=0; ii<size[0]; ii+=bs[0]) {
    for (std::int64_t jj=0; jj<size[1]; jj+=bs[1]) {
      for (std::int64_t kk=0; kk<size[2]; kk+=bs[2]) {
        float actual{0};
        reader->read(size3i_t{ii,jj,kk}, size3i_t{1,1,1}, &actual, 0);
        TEST_CHECK(similar(data, actual, 1.0e-5));
        data += 1;
      }
    }
  }
  reader->close();
  reader.reset();
}

static void
test_alturl()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  const std::string filename = cloud_synt2_name();
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils("sd://", Test_Utils::default_sd_context());
  const std::string alturl = utils->alturl(filename);
  if (verbose()) {
    std::cout << "alturl: " << alturl << std::endl;
  }
  TEST_CHECK(alturl.size() >= 3*filename.size());
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(filename, Test_Utils::default_sd_context());
  size3i_t actual_size = reader->size();
  size3i_t expect_size{181,241,169};
  TEST_CHECK(actual_size == expect_size);
}

/**
 * This test should work also if the token is expired.
 * It might not work for a bogus string, see
 * CallbackAuthProvider::getServiceAuthTokenImpl() and
 * utils:getAuthTokenExpiration()
 */
static void
test_idtoken()
{
  using InternalZGY::Environment;
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils("sd://", Test_Utils::default_sd_context());
  std::string token = utils->idtoken();
  std::string token_in_context = Environment::getStringEnv("OPENZGY_TOKEN");
  // If the token contains client credentials then idtoken() will
  // create a real token from it, so they are not expected to match.
  // This applies to other SAuth specific credential types as well.
  // But those are not currently not expected in unit tests.
  if (token_in_context.substr(0, 14) == "ServiceSauthV2" ||
      token_in_context.substr(0, 21) == "SlbAuthServiceAccount") {
    TEST_CHECK(!token.empty());
  }
  else {
    TEST_CHECK(token == token_in_context);
  }
}

static void
test_basicinfo()
{
  const std::string filename = cloud_synt2_name();
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils("sd://", Test_Utils::default_sd_context());
  std::map<std::string,std::string> info = utils->basicinfo(filename);
  if (verbose()) {
    std::stringstream ss;
    ss << "\ninfo(\"" << filename << "\")\n";
    for (const auto& it : info)
      ss << "     \"" << it.first << "\": \"" << it.second << "\"\n";
    std::cout << ss.str() << std::flush;
  }
  const auto tag = info.find("legaltag");
  if (TEST_CHECK(tag != info.end()))
    TEST_EQUAL(tag->second, "opendes-slb-synthetic-seismic");
}

static void
test_sharecred()
{
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils("sd://", Test_Utils::default_sd_context());
  TEST_EQUAL(utils.use_count(), 1);
  {
    auto ctxt = SeismicStoreIOContext().credentialsFrom(utils);
    TEST_EQUAL(utils.use_count(), 2); // utils itself, and a lambda inside ctxt.
    Test_Utils::CloudFileAutoDelete cad("sharedcred.zgy", &ctxt);
    TEST_EQUAL(utils.use_count(), 3); // Added the ZgyUtils in "cad"
    do_write_once(cad.name(), &ctxt);
    do_check_written(cad.name(), &ctxt);
  }
  TEST_EQUAL(utils.use_count(), 1);
}

#endif

static void
test_noinfo()
{
  const std::string filename = get_testdata("Empty-v3.zgy");
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils(filename, nullptr);
  std::map<std::string,std::string> info = utils->basicinfo(filename);
  // No extra info has been implemented for plain files.
  TEST_EQUAL(info.size(), 0);
}

static void
test_historange()
{
  LocalFileAutoDelete lad("testhisto.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename(lad.name())
    .size(33, 28, 92)
    .datatype(SampleDataType::float32);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  std::vector<float> data(2*3*4, 0);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t bottom{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t bsize{64,64,64};
  const OpenZGY::IZgyWriter::size3i_t count{2,3,4};
  float fortytwo{42};
  writer->writeconst(origin, bsize, &fortytwo);
  data[0] = -1000;
  writer->write(bottom, count, data.data());
  data[0] = -500;
  writer->write(bottom, count, data.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate}, nullptr);
  writer->close();

  // The writer has seen "42" (as a constant brick), "0", "-500", and "-1000".
  // The latter was overwritten so the true value range is now -500..+42
  // but by design and for implementation reasons the histogram range will be
  // large enough to also hold the "-1000". A narrower range than -500..+42
  // (e.g. because I forgot to include the writeconst) would be a bug.
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(lad.name());
  const SampleHistogram h = reader->histogram();
  TEST_CHECK(h.minvalue <= -500);
  TEST_CHECK(h.maxvalue >= +42);
  // SINGLEPASS: The next two are implememtation details.
  //TEST_EQUAL_FLOAT(h.minvalue, -1000, 1e-5);
  //TEST_EQUAL_FLOAT(h.maxvalue, +42, 1e-5);
  TEST_EQUAL(h.samplecount, 33*28*92);
  reader->close();
}

static void
test_lod(OpenZGY::DecimationType decimation)
{
  LocalFileAutoDelete lad("testlods.zgy");
  ZgyWriterArgsV2 args = ZgyWriterArgsV2()
    .decimation(std::vector<OpenZGY::DecimationType>{decimation})
    .filename(lad.name())
    .size(3, 2, 101)
    .datatype(SampleDataType::float32);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  std::vector<float> data(3*2*101, 100);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t ssize{3,2,101};
  const OpenZGY::IZgyWriter::size3i_t lod1size{2,1,51};
  data[16] = 50; // cube[0,0,16]
  data[2*(101*2)+0*(101)+32] = 200; // cube[2,0,32]
  writer->write(origin, ssize, data.data());
  writer->finalize();
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(lad.name());
  std::unique_ptr<float[]> lowres(new float[2*1*51]);
  reader->read(origin, lod1size, lowres.get(), 1);

  // The only input samples that differ from 100 are [0,0,16] and [2,0,32]
  // corresponding to [0,0,8] and [1,0,16] in LOD 1. So the asserts
  // check those two and (for the first value) some samples above and below.

  // Expected result around [0,0,8].
  // The first "expect" is based on a previous run but you can see
  // that this looks like a lowpass filter. Ditto for the second which
  // should be like average, just one sample != 100, but lower because
  // the single "50" input carries more weight. The last one is computed
  // trivially by hand.
  static double expect_zero[16]{0};
  static double expect_lowpass[16]
    {100.0,  100.0, 100.0, 100.0, 100.0, 100.0, 103.62, 93.39,
     76.64, 105.51, 95.84, 100.0, 100.0, 100.0,  100.0, 100.0};
  static double expect_weighted[16]
    {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
     50.57, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
  static double expect_average[16]
    {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
     93.75, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};

  const double* expect = expect_zero;
  switch (decimation) {
  case OpenZGY::DecimationType::LowPass:         expect = expect_lowpass; break;
  case OpenZGY::DecimationType::WeightedAverage: expect = expect_weighted; break;
  case OpenZGY::DecimationType::Average:         expect = expect_average; break;
  default: break;
  }

  bool decimation_ok{true};
  for (int ii=0; ii<16; ++ii)
    decimation_ok = decimation_ok && similar(lowres[ii], expect[ii], 0.02);
  TEST_CHECK(decimation_ok);
  if (!decimation_ok || verbose()) {
    std::cout << "    | ix:   expect     actual |\n"
              << "    | --  --------   -------- |\n";
    for (int ii=0; ii<16; ++ii)
      std::cout << "    | "
                << std::setw(2) << ii << ": "
                << std::setw(8) << expect[ii] << " | "
                << std::setw(8) << lowres[ii] << " |\n";
    std::cout << std::flush;
  }

  // Expected result around [1,0,16].
  // Source is [2,0,32]=200, [2,1,32]=100, [2,0,33]=100, [2,1,33]=100
  // and excludes [3,*,*] because that is outside the edge.
  // Average=125, weighted will be much closer to 200 because this
  // is a very rare value.
  const double value = lowres[1*(1*51)+0*(51)+16]; // [1,0,16]
  if (verbose())
    std::cout << "  " << value << "\n";
  switch (decimation) {
  case OpenZGY::DecimationType::LowPass:
    TEST_EQUAL_FLOAT(value, 146.73, 0.02);
    break;
  case OpenZGY::DecimationType::WeightedAverage:
    // WeightedAverage as lod 1 isn't really supported.
    // With the new plan A the result will not be deterministic,
    // because bricks are processed in parallel in an arbitrary
    // order, and "histogram so far" will vary. The same problem
    // occurs even for lod 2 if we implement incremental compute.
    // But at least that will only depend on the order of writes
    // from the application. So, slightly more deterministic.
    TEST_EQUAL_FLOAT(value, 199.5, 1.20);
    break;
  case OpenZGY::DecimationType::Average:
    TEST_EQUAL_FLOAT(value, 125.0, 0.02);
    break;
  default:
    TEST_CHECK(false && "missing test case");
    break;
  }

  reader->close();
}

static void
test_lod_lowpass()
{
  test_lod(OpenZGY::DecimationType::LowPass);
}

static void
test_lod_weighted()
{
  test_lod(OpenZGY::DecimationType::WeightedAverage);
}

static void
test_lod_average()
{
  test_lod(OpenZGY::DecimationType::Average);
}

static void
test_copy()
{
  std::string fname = get_testdata("Fancy-int8.zgy");
  LocalFileAutoDelete lad("testcopy.zgy");
  std::string oname = lad.name();

  // read
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> size = reader->size();
  std::unique_ptr<float[]>buf(new float[size[0] * size[1] * size[2]]);
  reader->read(orig, size, buf.get(), 0);

  // write
  ZgyWriterArgs args = ZgyWriterArgs()
    .metafrom(reader)
    //.datatype(SampleDataType::float32)
    .filename(oname);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  writer->write(origin, size, buf.get());
  writer->finalize(std::vector<OpenZGY::DecimationType>{
        OpenZGY::DecimationType::LowPass,
        OpenZGY::DecimationType::WeightedAverage
      }, nullptr);
  reader->close();
  writer->close();

  // check
  // TODO-High: FAILS: I can't find any LOD algorithm that makes the copy
  // match the input. The input was probably written from Python with
  // default finalize(), see "fancy-1.zgy" in test/black.py. That code can
  // also produce "fancy-5.zgy" which is the same file but written from
  // the old accessor. So; do the algorithms differ?
  // TODO-Test: Using a float file would have been better for comparing
  // lod generation between OpenZGY/Python and OpenZGY/C++.

  Test_Utils::compare_files(fname, oname, 1.0e-5, 1.0e+30);
}

static void
test_enums()
{
  std::stringstream ss;
  ss << SampleDataType::int8 << ", "
     << UnitDimension::length << ", "
     << DecimationType::Median;
  TEST_CHECK(ss.str() == "SampleDataType::int8, UnitDimension::length, DecimationType::Median");
}

namespace {
  /**
   * Compressor function that always returns "cannot compress", so the
   * file remains uncompressed. The instance keeps track of how many
   * times the low level code attempted a compression. So the unit
   * tests can check that the code did in fact try.
   *
   * The instance itself is noncopyable; otherwise it would be
   * pointless to keep state. This means that although the signature
   * is correct it will need to be wrapped in a lambda that can be
   * copied while still pointing to the same DummyCompressor instance.
   * get() will return a suitable lambda.
   */
  class DummyCompressor
  {
  public:
    int called;
    DummyCompressor(const DummyCompressor&) = delete;
    const DummyCompressor& operator=(const DummyCompressor&) = delete;
    DummyCompressor() : called(0)
    {
    }
    IZgyWriter::rawdata_t operator()(const IZgyWriter::rawdata_t&, const IZgyWriter::size3i_t&) {
      ++called;
      return IZgyWriter::rawdata_t{nullptr, 0};
    }
    IZgyWriter::compressor_t get()
    {
      return [this](const IZgyWriter::rawdata_t& data, const IZgyWriter::size3i_t& shape) {return (*this)(data, shape);};
    }
  };
} // namespace

static void
test_dummy_compress()
{
  LocalFileAutoDelete lad("dummycompress.zgy");
  DummyCompressor compressor, lodcompressor;
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename(lad.name())
    .size(33, 28, 92)
    .datatype(SampleDataType::float32)
    .compressor(compressor.get())
    .lodcompressor(lodcompressor.get());
  std::shared_ptr<IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  std::vector<float> data = Test_Utils::random_vector(33*28*92LL);
  const IZgyWriter::size3i_t origin{0,0,0};
  const IZgyWriter::size3i_t size = writer->size();
  writer->write(origin, size, data.data());
  TEST_CHECK(compressor.called == 2);
  writer->finalize();
  TEST_CHECK(lodcompressor.called == 1);
  writer->close();
}

struct onevalue_t
{
  double range_lo, range_hi;
  double stats_lo, stats_hi;
  double histo_lo, histo_hi;
  std::int64_t stats_count;
  std::int64_t histo_count;
  std::vector<std::int64_t> bins;

  onevalue_t(const std::array<float,2>& range,
             const SampleStatistics& stats,
             const SampleHistogram& histo)
    : range_lo(range[0]), range_hi(range[1])
    , stats_lo(stats.min), stats_hi(stats.max)
    , histo_lo(histo.minvalue), histo_hi(histo.maxvalue)
    , stats_count(stats.cnt)
    , histo_count(histo.samplecount)
    , bins(histo.bins)
  {
  }
};

static onevalue_t
test_histo_onevalue(SampleDataType dtype, float value, bool fill, const std::array<float,2>& datarange)
{
  if (verbose())
    std::cout << "\nTest dtype " << (int)dtype
              << " value " << value
              << (fill ? " only" : " and unwritten bricks")
              << "\n";
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("testhisto.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .filename(lad.name())
    .size(64,64,3*64)
    .datatype(dtype) // if float, datarange will be ignored.
    .datarange(datarange[0], datarange[1]);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  if (std::isfinite(value)) {
    std::vector<float> buf(64*64*64, value);
    writer->write(size3i_t{0,0,0}, size3i_t{64,64,64}, buf.data());
  }
  if (fill && std::isfinite(value)) {
    std::vector<float> buf(64*64*128, value);
    writer->write(size3i_t{0,0,64}, size3i_t{64,64,128}, buf.data());
  }
  writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate}, nullptr, OpenZGY::FinalizeAction::BuildFull, /*force=*/true);
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(lad.name());

  if (verbose()) {
    std::cout << "Data range " << reader->datarange()[0]
              << " " << reader->datarange()[1] << "\n";
    std::cout << "Statistics " << reader->statistics() << "\n";
    std::cout << "Histogram  " << reader->histogram() << "\n";
    std::cout << std::flush;
  }
  return onevalue_t(reader->datarange(),
                    reader->statistics(),
                    reader->histogram());
}

static onevalue_t
test_histo_onevalue(SampleDataType dtype, float value, bool fill)
{
  float center = std::isfinite(value) ? value : -0.25f;
  return test_histo_onevalue(dtype, value, fill,
                             std::array<float,2>{center-1, center+1});
}

#if 0
/*
 * SINGLEPASS: Most of these subtle issues work differently with dynamic
 * hisogram. This means that the test is not useful for testing regressions,
 * and since there are other unit tests for the new dynamic histograms
 * I'll just disable this test. If you really need it:...
 *  - Histogram ranges can be wider
 *  - Statistics and/or histogram might include the value "zero".
 *  - Disabling trimRange() might have make a difference.
 */
static void
test_histo_cornercase_float()
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const std::int64_t BRICK = 64*64*64;

  // Float: datarange with zero size is valid on input,
  // in fact the data range isn't specified by the user.
  // Reading back data gives the statistical range
  // which for float may include defaultvalue.
  // The histogram will use the fuzzy algorithm.

  // The numbers in brackets correspond to the ones in
  // GenLodImpl::suggestHistogramRange().

  onevalue_t r = test_histo_onevalue(SampleDataType::float32, nan, false);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==    0 && r.stats_hi ==    0);
  TEST_CHECK(r.histo_lo == -128 && r.histo_hi == +127);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[128] == r.histo_count);

  // [4] one all zero brick, two never written.
  // Expected result same as for nothing written.
  r = test_histo_onevalue(SampleDataType::float32, 0, false);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==    0 && r.stats_hi ==    0);
  TEST_CHECK(r.histo_lo == -128 && r.histo_hi == +127);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[128] == r.histo_count);

  // [4] three all zero bricks.
  // Expected result same as for nothing written.
  r = test_histo_onevalue(SampleDataType::float32, 0, true);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==    0 && r.stats_hi ==    0);
  TEST_CHECK(r.histo_lo == -128 && r.histo_hi == +127);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[128] == r.histo_count);

  // [6] single negative value, plus two never written bricks.
  // The statistics and histogram include the never-written
  // samples as if they were zero.
  // Note: I won't be testing the "some never written" scenario
  // for every remaining case; it is hopefully enough to
  // confirm once that never-written is treated as written-zero.
  r = test_histo_onevalue(SampleDataType::float32, -42, false);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==  -42 && r.stats_hi ==    0);
  TEST_EQUAL_FLOAT(r.histo_lo, -42, 0.001);
  TEST_EQUAL_FLOAT(r.histo_hi, 0, 0.001);
  TEST_EQUAL(r.stats_count, 3*BRICK);
  // Histogram gets fully utilized, so all our samples end up in first/last.
  TEST_EQUAL(r.bins[0], BRICK);
  TEST_EQUAL(r.bins[255], 2*BRICK);

  // [6] single negative value in all three bricks.
  // The value range and the statistics should have the true
  // range i.e. low==high and the histogram range should be wider.
  r = test_histo_onevalue(SampleDataType::float32, -42, true);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==  -42 && r.stats_hi ==  -42);
  TEST_EQUAL_FLOAT(r.histo_lo, -42, 0.001);
  TEST_EQUAL_FLOAT(r.histo_hi, 0, 0.001);
  TEST_EQUAL(r.stats_count, 3*BRICK);
  // Histogram gets fully utilized, so all our samples end up in first/last.
  TEST_EQUAL(r.bins[0], 3*BRICK);
  TEST_EQUAL(r.bins[255], 0);

  // [6] single positive value in all three bricks.
  // Result similar to the above but the ranges are swapped.
  r = test_histo_onevalue(SampleDataType::float32, +42, true);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==   42 && r.stats_hi ==   42);
  TEST_EQUAL_FLOAT(r.histo_lo,  0, 0.001);
  TEST_EQUAL_FLOAT(r.histo_hi, 42, 0.001);
  TEST_EQUAL(r.stats_count, 3*BRICK);
  // Histogram gets fully utilized, so all our samples end up in first/last.
  TEST_EQUAL(r.bins[0], 0);
  TEST_EQUAL(r.bins[255], 3*BRICK);
}
#endif

static void
test_histo_cornercase_int()
{
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const std::int64_t BRICK = 64*64*64;

  // Integral data.
  // Histogram range should always match the user provided range,
  // which for never-written is -1.25 to +0.75 and for the
  // remaining cases value +/- 1. This means that value won't be
  // exactly representable as an integer (it maps to -0.5) and
  // this will be noticeable in the statistics. The 0.5 factor
  // may also lead to numerical instability. The samples end up
  // either in bin 127 or bin 128.
  // Also, range might be wider then statistics (unlike the float
  // case) if not all possible storage values have been used.
  onevalue_t r = test_histo_onevalue(SampleDataType::int8, nan, false);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 0) < 0.25);
  TEST_CHECK(std::abs(r.stats_lo - 0) > 0.001); // 0.0 not representable.
  TEST_CHECK(r.histo_lo == -1.25 && r.histo_hi == 0.75); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  // I don't really care where the "0" samples end up. It won't be the center.
  TEST_CHECK(r.bins[127] + r.bins[128] == 0);

  r = test_histo_onevalue(SampleDataType::int8, 0, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 0) < 0.25);
  TEST_CHECK(std::abs(r.stats_lo - 0) > 0.001); // 0.0 not representable.
  TEST_CHECK(r.histo_lo == 0-1 && r.histo_hi == 0+1); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  r = test_histo_onevalue(SampleDataType::int8, -42, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  // SINGLEPASS: Values that used to exist might still show in min/max range.
  // Especially is trimRange() is disabled. Remove the check completely,
  // so it doesn't break again if/when trimRange() gets put back.
  //TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo + 42) < 0.25);
  TEST_CHECK(std::abs(r.stats_lo + 42) > 0.001); // 42.0 not representable.
  TEST_CHECK(r.histo_lo == -42-1 && r.histo_hi == -42+1);// user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  r = test_histo_onevalue(SampleDataType::int8, +42, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  // SINGLEPASS: (see above) TEST_CHECK(r.stats_lo == r.stats_hi);
  //TEST_CHECK(std::abs(r.stats_lo - 42) < 0.25);
  //TEST_CHECK(std::abs(r.stats_lo - 42) > 0.001); // 42.0 not representable.
  TEST_CHECK(r.histo_lo == 42-1 && r.histo_hi == 42+1); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  // 16 bit not much different from 8 bit, but the statistics will be
  // closer to the supplied value because the quantization error is smaller.
  r = test_histo_onevalue(SampleDataType::int16, nan, false);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 0) < 0.25/256);
  TEST_CHECK(std::abs(r.stats_lo - 0) > 0.001/256); // 0.0 not representable.
  TEST_CHECK(r.histo_lo == -1.25 && r.histo_hi == 0.75); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  // I don't really care where the "0" samples end up. It won't be the center.
  TEST_CHECK(r.bins[127] + r.bins[128] == 0);

  r = test_histo_onevalue(SampleDataType::int16, 0, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 0) < 0.25/256);
  TEST_CHECK(std::abs(r.stats_lo - 0) > 0.001/256); // 0.0 not representable.
  TEST_CHECK(r.histo_lo == 0-1 && r.histo_hi == 0+1); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  r = test_histo_onevalue(SampleDataType::int16, -42, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  // SINGLEPASS: (see above) TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo + 42) < 0.25/256);
  TEST_CHECK(std::abs(r.stats_lo + 42) > 0.001/256); // 42.0 not representable.
  TEST_CHECK(r.histo_lo == -42-1 && r.histo_hi == -42+1);// user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  r = test_histo_onevalue(SampleDataType::int16, +42, true);

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  // SINGLEPASS: (see above) TEST_CHECK(r.stats_lo == r.stats_hi);
  //TEST_CHECK(std::abs(r.stats_lo - 42) < 0.25/256);
  //TEST_CHECK(std::abs(r.stats_lo - 42) > 0.001/256); // 42.0 not representable.
  TEST_CHECK(r.histo_lo == 42-1 && r.histo_hi == 42+1); // user choice exactly.
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[127] + r.bins[128] == 3*BRICK);

  // Behavior when all explicitly written values get clipped.
  // Expect both the histogram and the statistics to only reflect
  // the clipped value (-5) as if that value and not -42 had been
  // written.
  r = test_histo_onevalue(SampleDataType::int8, -42, true,
                          std::array<float,2>{-5, +760});

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  // SINGLEPASS: (see above) TEST_CHECK(r.stats_lo == -5 && r.stats_hi == -5);
  TEST_CHECK(r.histo_lo == -5 && r.histo_hi == +760);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[0] == 3*BRICK);

  // As above, all explicitly written values get clipped but now
  // there are a few unwritten bricks. Expect both the histogram
  // and the statistics to only reflect the clipped value (-5) as
  // if that value and not -42 had been written.
  // Defaultvalue is +1 because the range does not give a zero
  // centric histogram. The statistics should also reflect that.
  // I.e. expect +1 to be part of the range.
  r = test_histo_onevalue(SampleDataType::int8, -42, false,
                          std::array<float,2>{-5, +760});

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
    // Data dependent
  TEST_CHECK(r.stats_lo == -5 && r.stats_hi == +1);
  TEST_CHECK(r.histo_lo == -5 && r.histo_hi == +760);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[0] == BRICK);
  TEST_CHECK(r.bins[2] == 2*BRICK);

  // Similar to the above but no values written at all.
  // Defaultvalue is still 1 due to missing zero-centric propery
  // so this is what should be reflected in the statistics.
  r = test_histo_onevalue(SampleDataType::int8, nan, false,
                          std::array<float,2>{-5, +760});

  // Invariants for the integer case
  TEST_CHECK(r.range_lo <= r.stats_lo && r.range_hi >= r.stats_hi);
  TEST_CHECK(r.histo_lo == r.range_lo && r.histo_hi == r.range_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  // Data dependent
  TEST_CHECK(r.stats_lo == +1 && r.stats_hi == +1);
  TEST_CHECK(r.histo_lo == -5 && r.histo_hi == +760);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[2] == 3*BRICK);
}

/**
 * filestats() is used to extract human-readable information about a file.
 * Also numerical values to be used for display purposes only.
 */
static void
test_filestats()
{
  std::string fname = get_testdata("Fancy-int8.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  std::shared_ptr<const FileStatistics> stat = reader->filestats();

  std::stringstream ss;
  stat->dump(ss);
  std::string msg = ss.str();
  if (verbose()) {
    std::cout << msg;
    std::cout << "Id:    " << reader->verid() << std::endl;
  }

  TEST_CHECK(stat->fileVersion() == 3);
  TEST_CHECK(stat->fileSize() == 1310720);
  TEST_CHECK(msg.find("Size:  1310720") != std::string::npos);
  TEST_CHECK(msg.find("Brick: 2 missing") != std::string::npos);
  TEST_CHECK(stat->headerSize() == 2569);
  TEST_CHECK(stat->alphaNormalCount() == 0);
  TEST_CHECK(stat->alphaNormalSizePerEntry() == 64*64);
  TEST_CHECK(stat->alphaCompressedCount() == 0);
  TEST_CHECK(stat->alphaCcompressedSize() == 0);
  TEST_CHECK(stat->alphaMissingCount() == 4); // 2 lod0, 1 lod1, 1 lod2
  TEST_CHECK(stat->alphaConstantCount() == 0);
  TEST_CHECK(stat->brickNormalCount() == 4); // 2 lod0, 1 lod1, 1 lod2
  TEST_CHECK(stat->brickNormalSizePerEntry() == 64*64*64*sizeof(std::uint8_t));
  TEST_CHECK(stat->brickCompressedCount() == 0);
  TEST_CHECK(stat->brickCompressedSize() == 0);
  TEST_CHECK(stat->brickMissingCount() == 2); // 2 in lod0.
  TEST_CHECK(stat->brickConstantCount() == 3); // 1 in lod1 (from the empties)
  TEST_CHECK(stat->usedSize() == stat->headerSize() + stat->brickNormalCount() * stat->brickNormalSizePerEntry());
  TEST_CHECK(stat->usedIfUncompressed() == stat->usedSize());
  TEST_CHECK(stat->compressionFactor() == 1);
  TEST_CHECK(stat->isCompressed() == false);

  // Check some other meta information while we are here.
  TEST_CHECK(reader->verid().size() == 36);
  TEST_CHECK(reader->verid() != "00000000-0000-0000-0000-000000000000");

  TEST_CHECK(reader->datarange()[0] == reader->raw_datarange()[0]);
  TEST_CHECK(reader->datarange()[1] == reader->raw_datarange()[1]);

  reader->close();
}

/**
 * Test coordiname transforms.
 * Can be done both on files open for read and open for write.
 * Both cases should be tested because in the write case all
 * the calls will go via ZgySafeWriter.
 */
static void
do_test_transform(const IZgyTools& reader)
{
  typedef OpenZGY::IZgyReader::corners_t corners_t;

  const double size[2]{112,64};
  const double annotbeg[2]{1234, 5678};
  const double annotinc[2]{5, 2};
  const double annotend[2]{
    annotbeg[0] + (size[0]-1) * annotinc[0],
    annotbeg[1] + (size[1]-1) * annotinc[1]};

  const corners_t expect_index{{
    {        0,         0},
    {size[0]-1,         0},
    {        0, size[1]-1},
    {size[0]-1, size[1]-1}}};
  const corners_t expect_annot{{
    {annotbeg[0], annotbeg[1]},
    {annotend[0], annotbeg[1]},
    {annotbeg[0], annotend[1]},
    {annotend[0], annotend[1]}}};
  const corners_t expect_world{{
    {1000, 1000},
    {3775, 1000},
    {1000, 2890},
    {3775, 2890}}};

  const corners_t actual_index = reader.indexcorners();
  const corners_t actual_annot = reader.annotcorners();
  const corners_t actual_world = reader.corners();

  for (int ii=0; ii<4; ++ii) {
    for (int jj=0; jj<2; ++jj) {
      TEST_CHECK(expect_index[ii][jj] == actual_index[ii][jj]);
      TEST_CHECK(expect_annot[ii][jj] == actual_annot[ii][jj]);
      TEST_CHECK(expect_world[ii][jj] == actual_world[ii][jj]);
    }
  }

  typedef std::array<double,2> point_t;
  // Annot increments (2,5), world increments (25,30), azimuth 0.
  TEST_CHECK((reader.indexToAnnot(point_t{{1,2}}) == point_t{1239,5682}));
  TEST_CHECK((reader.annotToIndex(point_t{{1239,5682}}) == point_t{1,2}));

  TEST_CHECK((reader.indexToWorld(point_t{{1,2}}) == point_t{1025,1060}));
  TEST_CHECK((reader.worldToIndex(point_t{{1025,1060}}) == point_t{1,2}));

  TEST_CHECK((reader.annotToWorld(point_t{{1239,5682}}) == point_t{1025,1060}));
  TEST_CHECK((reader.worldToAnnot(point_t{{1025,1060}}) == point_t{1239,5682}));
}

static void
test_transform_r()
{
  std::string fname = get_testdata("Fancy-int8.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  do_test_transform(*reader);
  reader->close();
}

static void
test_transform_w()
{
  std::string fname = get_testdata("Fancy-int8.zgy");
  LocalFileAutoDelete lad("testfile.zgy");
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(fname);
  ZgyWriterArgs args = ZgyWriterArgs()
    .metafrom(reader)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = OpenZGY::IZgyWriter::open(args);
  do_test_transform(*writer);
  writer->finalize();
  writer->close();
  reader->close();
}

static void
test_all_exceptions()
{
  using namespace OpenZGY::Errors;

  must_throw("Exception test", [&](){
    throw ZgyFormatError("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyCorruptedFile("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyUserError("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyInternalError("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyEndOfFile("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgySegmentIsClosed("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyAborted("Exception test");
  });
  must_throw("Exception test", [&](){
    throw ZgyMissingFeature("Exception test");
  });
  must_throw("Bogus file name:", [&](){
    throw ZgyIoError("Bogus file name", 2);
  });
  must_throw("Exception test", [&](){
    throw ZgyNotReadOnlyError("Exception test");
  });
}

static void
test_ambig1()
{
  LocalFileAutoDelete lad("ambig1.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(128,129,4)
    .datatype(SampleDataType::int16)
    .datarange(-32768, +32767)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  const std::int16_t data1[4]{0,0,0,0};
  writer->write(IZgyWriter::size3i_t{0,0,0}, IZgyWriter::size3i_t{1,1,4}, &data1[0]);
  writer->finalize();
  writer->close();
}

/**
 * In this test the ambiguity occurs because 4*short = 8 bytes are written,
 * which is the size of a double. This makes it look like a scalar.
 *
 * This test will create a huge file with almost all constant values.
 * If it starts running very slowly, this is a regression.
 * See also test_genlod2(), which is more focused on that problem.
 */
static void
test_ambig2()
{
  LocalFileAutoDelete lad("ambig2.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(131313,13,131)
    .datatype(SampleDataType::int16)
    .datarange(-32768, +32767)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  const std::int16_t fortytwo{42};
  writer->writeconst(IZgyWriter::size3i_t{0,0,0}, IZgyWriter::size3i_t{1,1,4}, &fortytwo);
  writer->finalize();
  writer->close();
}

/**
 * In this test the ambiguity occurs when a low resolution brick gets
 * so small that size in the inline direction becomes 1 and the data
 * then looks like it is 2d.
 */
static void
test_ambig3()
{
  LocalFileAutoDelete lad("ambig3.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(2,20,256)
    .datatype(SampleDataType::int16)
    .datarange(-32768, +32767)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  const std::int16_t fortytwo{42};
  writer->writeconst(IZgyWriter::size3i_t{0,0,0}, IZgyWriter::size3i_t{1,1,4}, &fortytwo);
  writer->finalize();
  writer->close();
}

static bool
do_2d(const IZgyWriter::size3i_t& bs, const IZgyWriter::size3i_t& size, bool compress = false)
{
  typedef IZgyWriter::size3i_t size3i_t;
  bool ok = true;
  LocalFileAutoDelete lad("2d.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .filename(lad.name());
  if (compress)
    args.zfp_compressor(65);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  const float fortytwo{42};
  const float onehundred{100};
  writer->writeconst(size3i_t{0,0,0}, size, &fortytwo);
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &onehundred);
  writer->finalize();
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  ok = TEST_CHECK(reader->nlods() > 1) && ok;

  float data[4]{0};
  size3i_t readsize{1,1,1};
  if      (size[2]>=4) readsize[2]=4;
  else if (size[1]>=4) readsize[1]=4;
  else if (size[0]>=4) readsize[0]=4;
  reader->read(size3i_t{0,0,0}, readsize, data, 0);
  ok = TEST_EQUAL_FLOAT(data[0], 100.0f, 1.0f) && ok;
  if (readsize[0] * readsize[1] * readsize[2] >= 4) {
    ok = TEST_EQUAL_FLOAT(data[1],  42.0f, 1.0f) && ok;
    ok = TEST_EQUAL_FLOAT(data[2],  42.0f, 1.0f) && ok;
    ok = TEST_EQUAL_FLOAT(data[3],  42.0f, 1.0f) && ok;
  }
  return ok;
}

static void
test_2d()
{
  typedef IZgyWriter::size3i_t size3i_t;
  const bool allow_2d = InternalZGY::Environment::getNumericEnv("OPENZGY_ALLOW_2D", 0) > 0;

  if (verbose()) {
    if (allow_2d)
      std::cout << "\ntest_2d is running all tests.\n" << std::flush;
    else
      std::cout << "\ntest_2d skipping some tests.\n" << std::flush;
  }

  // The lowpass decimation algorithm might have problems with very short
  // traces, even if they aren't size 1.
  TEST_CHECK(do_2d(size3i_t{64, 64, 64}, size3i_t{256, 20, 2}, false));
  TEST_CHECK(do_2d(size3i_t{64, 64, 64}, size3i_t{256, 20, 2}, true));

  // One dimemsion is 1 only in size. The corresponding bricksize is 4,
  // to avoid potential corner cases especially in compressed files,
  TEST_CHECK(do_2d(size3i_t{ 4, 64, 64}, size3i_t{   1, 64+5, 64+7}, false));
  TEST_CHECK(do_2d(size3i_t{64,  4, 64}, size3i_t{64+3,    1, 64+7}, false));
  TEST_CHECK(do_2d(size3i_t{64, 64,  4}, size3i_t{64+3, 64+5,    1}, false));

  // Enable compression for the above
  TEST_CHECK(do_2d(size3i_t{ 4, 64, 64}, size3i_t{   1, 64+5, 64+7}, true));
  TEST_CHECK(do_2d(size3i_t{64,  4, 64}, size3i_t{64+3,    1, 64+7}, true));
  TEST_CHECK(do_2d(size3i_t{64, 64,  4}, size3i_t{64+3, 64+5,    1}, true));

  if (allow_2d) {
    // One dimemsion (the same one) is 1 for both size and bricksize.
    // This would be the typical case for storing 2d data.
    TEST_CHECK(do_2d(size3i_t{ 1, 64, 64}, size3i_t{   1, 64+5, 64+7}));
    TEST_CHECK(do_2d(size3i_t{64,  1, 64}, size3i_t{64+3,    1, 64+7}));
    TEST_CHECK(do_2d(size3i_t{64, 64,  1}, size3i_t{64+3, 64+5,    1}));

    // Does compression actually work in this case?
    TEST_CHECK(do_2d(size3i_t{ 1, 64, 64}, size3i_t{   1, 64+5, 64+7}, true));
    TEST_CHECK(do_2d(size3i_t{64,  1, 64}, size3i_t{64+3,    1, 64+7}, true));
    TEST_CHECK(do_2d(size3i_t{64, 64,  1}, size3i_t{64+3, 64+5,    1}, true));

    // One dimemsion is 1 only in bricksize. This might be used not for
    // 2d data but for creating a slice-optimized cube. Not really useful
    // because tweaking the order of bricks written to the file will likely
    // give better results.
    TEST_CHECK(do_2d(size3i_t{ 1, 64, 64}, size3i_t{  11, 64+5, 64+7}));
    TEST_CHECK(do_2d(size3i_t{64,  1, 64}, size3i_t{64+3,   11, 64+7}));
    TEST_CHECK(do_2d(size3i_t{64, 64,  1}, size3i_t{64+3, 64+5,   11}));

    // One dimemsion (not the same one) is 1 for both size and bricksize.
    // Completely silly but useful as a monkey test.
    TEST_CHECK(do_2d(size3i_t{64,  1, 64}, size3i_t{   1, 64+5, 64+7}));
    TEST_CHECK(do_2d(size3i_t{64, 64,  1}, size3i_t{64+3,    1, 64+7}));
    TEST_CHECK(do_2d(size3i_t{ 1, 64, 64}, size3i_t{64+3, 64+5,    1}));

    // One-dimensional data, anybody?
    TEST_CHECK(do_2d(size3i_t{64,   1,  1}, size3i_t{64+3,    1,    1}));
    TEST_CHECK(do_2d(size3i_t{ 1,  64,  1}, size3i_t{   1, 64+5,    1}));
    TEST_CHECK(do_2d(size3i_t{ 1,   1, 64}, size3i_t{   1,    1, 64+7}));

    // Zero dimensions.
    TEST_CHECK(do_2d(size3i_t{ 1,   1,  1}, size3i_t{   1,    1,    1}));
  }
}

static void
test_decimate_edge()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("decimate_edge.zgy");
  const size3i_t size{19, 39, 61};
  const size3i_t   bs{16, 32, 64};
  ZgyWriterArgsV2 args = ZgyWriterArgsV2()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(SampleDataType::int8)
    .datarange(-128, +127)
    .decimation(std::vector<OpenZGY::DecimationType>{
        OpenZGY::DecimationType::Average})
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  std::vector<std::int8_t> data(size[0]*size[1]*size[2], 10);
  for (int ii=0; ii<size[0]; ii += 2)
    for (int jj=0; jj<size[1]; jj += 2)
      for (int kk=0; kk<size[2]; kk += 2)
        data[ii*size[1]*size[2] + jj*size[2] + kk] = 90;
  writer->write(size3i_t{0,0,0}, size, data.data());
  writer->finalize();
  writer->close();

  const size3i_t lod1size{(size[0]+1)/2, (size[1]+1)/2, (size[2]+1)/2};
  std::vector<std::int8_t> check(lod1size[0]*lod1size[1]*lod1size[2], -1);
  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  reader->read(size3i_t{0,0,0}, lod1size, check.data(), 1);
  const auto offset = [](const size3i_t& size, std::int64_t ii, std::int64_t jj, std::int64_t kk) {
                        return ii*size[1]*size[2]+jj*size[2]+kk;
                   };
  // Each lod1 sample was computed as the simple average of one "90" sample
  // and 7, 3, 1, or 0 "10" samples depending on how close we are to the edge.
  // Expected average values: 20, 30, 50, 90.
  const size3i_t last{lod1size[0]-1, lod1size[1]-1, lod1size[2]-1};
  TEST_EQUAL((int)check[offset(lod1size,       1,       2,       3)], 20);
  TEST_EQUAL((int)check[offset(lod1size, last[0],       2,       3)], 30);
  TEST_EQUAL((int)check[offset(lod1size,       1, last[1],       3)], 30);
  TEST_EQUAL((int)check[offset(lod1size,       1,       2, last[2])], 30);
  TEST_EQUAL((int)check[offset(lod1size,       1, last[1], last[2])], 50);
  TEST_EQUAL((int)check[offset(lod1size, last[0],       2, last[2])], 50);
  TEST_EQUAL((int)check[offset(lod1size, last[0], last[1],       3)], 50);
  TEST_EQUAL((int)check[offset(lod1size, last[0], last[1], last[2])], 90);
}

/**
 * Check that you can read from a file while it is still open for writing.
 */
static void
do_test_readwrite(const std::string& filename, const IOContext *context = nullptr)
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;

  std::shared_ptr<OpenZGY::IZgyWriter> writer =
    OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(context)
     .filename(filename)
     .size(192,192,192)
     .datatype(SampleDataType::int16)
     .datarange(-32768,+32767));
  const float value1{1000}, value2{2000}, value3{3000};
  writer->writeconst(size3i_t{0,0,0}, size3i_t{64,64,64}, &value1);
  writer->writeconst(size3i_t{64,64,64}, size3i_t{64,64,64}, &value2);
  writer->writeconst(size3i_t{96,96,96}, size3i_t{64,64,64}, &value3);

  // Read from the still open writer and check the results.
  std::vector<float> check(192*192*192, -999);
  writer->read(size3i_t{0,0,0}, size3i_t{192,192,192}, check.data());
  int numzero{0}, num1000{0}, num2000{0}, num3000{0};
  for (float value : check) {
    switch ((int)(value+0.001)) {
    case 0:    ++numzero; break;
    case 1000: ++num1000; break;
    case 2000: ++num2000; break;
    case 3000: ++num3000; break;
    }
  }
  TEST_EQUAL(num1000, 64*64*64);
  TEST_EQUAL(num2000, 64*64*64 - 32*32*32);
  TEST_EQUAL(num3000, 64*64*64);
  TEST_EQUAL(numzero + num1000 + num2000 + num3000, (int)check.size());
  writer->close_incomplete();
}

/**
 * Check that you can read from a local file while it is still open for writing.
 */
static void
test_readwrite_local()
{
  Test_Utils::LocalFileAutoDelete lad("readwrite_local.zgy");
  do_test_readwrite(lad.name(), nullptr);
}

#ifdef HAVE_SD
/**
 * Check that you can read from a cloud file while it is still open for writing.
 * This includes satisfying reads from the cache where appliccable.
 */
static void
test_readwrite_cloud()
{
  auto context = Test_Utils::default_sd_context();
  Test_Utils::CloudFileAutoDelete cad("readwrite_cloud.zgy", context);
  do_test_readwrite(cad.name(), context);
}

static void
read_hammer(std::shared_ptr<IZgyReader> reader, int repeats)
{
  const std::array<std::int64_t,3> start{0, 0, 0};
  const std::array<std::int64_t,3> size{64, 64, 64};
  std::unique_ptr<float[]> data(new float[size[0]*size[1]*size[2]]);
  for (int ii = 0; ii < repeats; ++ii)
    reader->read(start, size, data.get(), 0);
}

/**
 * Hammer the same brick in the same file again and again and ...
 * To use this for a real performance test you may need to increase
 * the number of threads and repeats.
 */
static void
test_hammer()
{
  using InternalZGY::Environment;
  const int repeats = Environment::getNumericEnv("OPENZGY_HAMMER_REPEATS", 8);
  const int threads = Environment::getNumericEnv("OPENZGY_HAMMER_THREADS", 8);
  const std::string name= Environment::getStringEnv("OPENZGY_HAMMER_NAME",
                                                    cloud_synt2_name().c_str());
  std::shared_ptr<IZgyReader> reader =
    IZgyReader::open(name, Test_Utils::default_sd_context());
  std::vector<std::thread> workers;
  for (int ii = 0; ii < threads; ++ii)
    workers.push_back(std::thread(read_hammer, reader, repeats));
  for (int ii = 0; ii < threads; ++ii)
    workers[ii].join();
}

/**
 * Most errors thrown by SDAPI now get caught and re-thrown as OpenZGY
 * exceptions. In some cases OpenZGY is able to provide a better error
 * message than SDAPI.
 */
static void
test_sderrors()
{
  const std::string filename = cloud_synt2_name();

  SeismicStoreIOContext context(*Test_Utils::default_sd_context());

  // SDAPI should not retry a missing token, but currently there
  // is a bug that does just that. Explicitly turn off retries.
  context.retryCount(0);

  // Neither token not token callback were provided.
  context.sdtoken("", "");
  must_throw("Missing access token or callback in iocontext", [&](){
    OpenZGY::IZgyReader::open(filename, &context);
  });

  // Token callback provided but it always returns empty.
  std::string token;
  std::function<std::string()> functor = [&token]() {return token;};
  context.sdtokencb(functor, "");
  must_throw("Seismic Store", [&](){
    OpenZGY::IZgyReader::open(filename, &context);
  });

}

#endif

namespace {
  template<typename T>
  static std::string
  formatMe(std::int64_t pos, std::int64_t ii, std::int64_t jj, std::int64_t kk, const T* value)
  {
    std::stringstream ss;
    ss << "Pos " << pos << " brick ("
       << ii << "," << jj << "," << kk
       << ") value ("
       << (double)value[0] << "," << (double)value[1] << "," << (double)value[2]
       << ")";
    return ss.str();
  }
}

/**
 * Make sure the contents of the padding area is not visible to the
 * application, and that the behavior is predictable and the same
 * regardless of whether the one brick at a time shortcut is active.
 * The contents on disk is an implementation detail, currently it will
 * have replication up to the next multiple of 4 and defaultvalue
 * after that.
 *
 * See ZgyInternalBulk::expeditedRead() and deliverOneBrick().
 * Currently the feature is NOT enabled, and the current behavior
 * (replication up to a multiple of 4, then defaultvalue) is
 * explicitly tested for.
 *
 * Test data size 2*3*7 bricks, (1*64)+43, (2*64)+34, (6*64)+21
 * In the last column, Brick 2 is missing and brick 3 is a constvalue
 * Use an int16 file with real 0 mapping to (default-)storage 10000
 */
static void
test_edgebricks()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("edgebricks.zgy");
  const size3i_t   bs{64, 64, 64};
  const size3i_t size{(1*64)+43, (2*64)+34, (6*64)+21};
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(SampleDataType::int16)
    .datarange(-32768-10000, +32767-10000)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  if (!TEST_CHECK(bool(writer)))
    return;
  std::vector<std::int16_t> data(bs[0]*bs[1]*bs[2]);
  for (std::size_t pos=0, end=data.size(); pos < end; ++pos)
    data[pos] = ((pos * 947) % 4057) + 42;
  writer->write(size3i_t{1*64, 2*64, 0*64}, bs, data.data());
  writer->writeconst(size3i_t{1*64, 2*64, 1*64}, bs, data.data()); // const
  //writer->write(size3i_t{1*64, 2*64, 2*64}, bs, data.data()); // skip
  writer->write(size3i_t{1*64, 2*64, 3*64}, bs, data.data());
  writer->write(size3i_t{1*64, 2*64, 4*64}, bs, data.data());
  writer->write(size3i_t{1*64, 2*64, 5*64}, bs, data.data());
  writer->write(size3i_t{1*64, 2*64, 6*64}, bs, data.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{
         OpenZGY::DecimationType::Average
      }, nullptr);
  writer->close();
  writer.reset();

  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  if (!TEST_CHECK(bool(reader)))
    return;
  if (verbose())
    reader->filestats()->dump(std::cout, "");

  // Read the last brick column.
  std::vector<std::int16_t> check(64*64*7*64, 888);
  std::function<std::size_t(int, int, int)> offset;
  offset = [](int ii, int jj, int kk){return (ii*64 + jj) * 7*64 + kk;};
  reader->read(size3i_t{1*64, 2*64, 0}, size3i_t{64, 64, 7*64},check.data());
  TEST_EQUAL(check[offset(0,0,1)], ((1 * 947) % 4057) + 42); // inside
  TEST_EQUAL(check[offset(0,0,67)], 42); // inside in const
  TEST_EQUAL(check[offset(0,0,131)], 10000); // inside in miss

  // The normal, general API does *not* hide the real contents of the brick.
  // The chech should return 1000 for all three samples, instead it will
  // return 2754 (replicating last sample inside) for the first three.
  TEST_EQUAL(check[offset(0,0,6*64+21)], 2754); // outside in Z
  TEST_EQUAL(check[offset(0,0,6*64+22)], 2754); // outside in Z
  TEST_EQUAL(check[offset(0,0,6*64+23)], 2754); // outside in Z
  TEST_EQUAL(check[offset(0,0,6*64+24)], 10000); // well outside in Z

  TEST_EQUAL(check[offset(0,34,0)], 4062); // outside in J
  TEST_EQUAL(check[offset(0,35,0)], 4062); // outside in J
  TEST_EQUAL(check[offset(0,36,0)], 10000); // well outside in J
  TEST_EQUAL(check[offset(0,36,0)], 10000); // well outside in J

  TEST_EQUAL(check[offset(43,0,0)], 1454); // outside in I
  TEST_EQUAL(check[offset(44,0,0)], 10000); // well outside in I
  TEST_EQUAL(check[offset(45,0,0)], 10000); // well outside in I
  TEST_EQUAL(check[offset(46,0,0)], 10000); // well outside in I

  TEST_EQUAL(check[offset(43,34,1*64)], 42); // outside but in constvalue
  TEST_EQUAL(check[offset(43,34,2*64)], 10000); // outside but in missing

  // Read one brick at a time: constvalue, missing, normal.
  std::vector<std::int16_t> check1(64*64*64);
  std::vector<std::int16_t> check2(64*64*64);
  std::vector<std::int16_t> check3(64*64*64);
  reader->read(size3i_t{1*64, 2*64, 1*64}, bs, check1.data());
  reader->read(size3i_t{1*64, 2*64, 2*64}, bs, check2.data());
  reader->read(size3i_t{1*64, 2*64, 6*64}, bs, check3.data());
  offset = [](int ii, int jj, int kk){return (ii*64 + jj) * 64 + kk;};

  TEST_EQUAL(check3[offset(0,0,0)], 42); // inside
  TEST_EQUAL(check3[offset(0,0,1)], ((1 * 947) % 4057) + 42); // inside

  TEST_EQUAL(check3[offset(0,0,21)], 2754); // outside in Z
  TEST_EQUAL(check3[offset(0,0,22)], 2754); // outside in Z
  TEST_EQUAL(check3[offset(0,0,23)], 2754); // outside in Z
  TEST_EQUAL(check3[offset(0,0,24)], 10000); // outside in Z

  TEST_EQUAL(check3[offset(0,34,0)], 4062); // outside in J
  TEST_EQUAL(check3[offset(0,35,0)], 4062); // outside in J
  TEST_EQUAL(check3[offset(0,36,0)], 10000); // outside in J
  TEST_EQUAL(check3[offset(0,37,0)], 10000); // outside in J

  TEST_EQUAL(check3[offset(43,0,0)], 1454); // outside in I
  TEST_EQUAL(check3[offset(44,0,0)], 10000); // outside in I
  TEST_EQUAL(check3[offset(45,0,0)], 10000); // outside in I
  TEST_EQUAL(check3[offset(46,0,0)], 10000); // outside in I

  for (const auto it : check1) // padding in const brick is the constant.
    if (it != 42)
      if (!TEST_EQUAL(it, 42))
        break;
  for (const auto it : check2) // empty brick.
    if (it != 10000)
      if (!TEST_EQUAL(it, 10000))
        break;
}

/**
 * Test that reading a file BAT i.e. one brick at a time works.
 * These now trigger some performance tweaks in the accessor.
 *
 * The test file has 4*5*6 = 120 bricks and is written in optimal order.
 * most bricks store (i/64,j/64,k/64) in first three bytes to make it
 * simple to check that the contents are correct. Every 7th brick will
 * be missing. Every 13th brick will be filled with a constant value
 * -(i+j+k)/64
 */
template<typename T, int bricksize, bool compression>
static void
do_testbat(const std::string& filename)
{
  typedef IZgyWriter::size3i_t size3i_t;
  const size3i_t   bs{bricksize, bricksize, bricksize};
  const size3i_t bsm1{bricksize-1, bricksize, bricksize};
  const size3i_t size{4*bs[0]-10, 5*bs[1]-20, 6*bs[2]-23};
  const SampleDataType dt = (sizeof(T) == 4 ? SampleDataType::float32 :
                             sizeof(T) == 2 ? SampleDataType::int16 :
                             SampleDataType::int8);
#ifdef HAVE_SD
  SeismicStoreIOContext sd_context(*Test_Utils::default_sd_context());
  sd_context.segsize(3); // Try to trigger read crossing seg boundary.
  sd_context.segsplit(7);
  const IOContext * const context = &sd_context;
#else
  const IOContext * const context = nullptr;
#endif
  ZgyWriterArgs args = ZgyWriterArgs()
    .iocontext(context)
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(dt)
    .datarange(-1, +1)
    .filename(filename);
  if (compression)
    args.zfp_compressor(99);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  if (!TEST_CHECK(bool(writer)))
    return;
  std::vector<T> data(bs[0]*bs[1]*bs[2], 42);
  {
    std::int64_t mypos{ 0 };
    for (std::int64_t ii = 0; ii < size[0]; ii += bs[0]) {
      for (std::int64_t jj = 0; jj < size[1]; jj += bs[1]) {
        for (std::int64_t kk = 0; kk < size[2]; kk += bs[2]) {
          if ((mypos % 7) == 0) {
            // missing brick
          }
          else if ((mypos % 13) == 0) {
            // Constant-value
            data[0] = static_cast<T>(-(ii / bs[0]) - (jj / bs[1]) - (kk / bs[2]));
            writer->writeconst(size3i_t{ ii,jj,kk }, bs, data.data());
          }
          else {
            data[0] = static_cast<T>(ii / bs[0]);
            data[1] = static_cast<T>(jj / bs[1]);
            data[2] = static_cast<T>(kk / bs[2]);
            writer->write(size3i_t{ ii,jj,kk }, bs, data.data());
          }
          ++mypos;
        }
      }
    }
  }

  writer->finalize(std::vector<OpenZGY::DecimationType>{
         OpenZGY::DecimationType::Average
      }, nullptr);
  writer->close();
  writer.reset();

  // Test plan:
  // Read one brick at a time using 10 threads. Can trigger both cases.
  // Not tested: Second open is for read/write.
  // Compression should trigger malloc hack, not the shortcut.
  // Compressed reads might also cross segment boundary.

  {
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      IZgyReader::open(filename, context);
    if (!TEST_CHECK(bool(reader)))
      return;
    if (verbose())
      reader->filestats()->dump(std::cout, "");
    InternalZGY::MTGuard guard;
#pragma omp parallel num_threads(10)
    {
      std::vector<T> check1(bs[0]*bs[1]*bs[2], 0);
      std::vector<T> check2(bs[0]*bs[1]*bs[2], 0);
#pragma omp for schedule(dynamic,1)
      for (int pos = 0; pos < 4*5*6; ++pos) {
        int tmppos = pos;
        const int kk = tmppos % 6; tmppos /= 6;
        const int jj = tmppos % 5; tmppos /= 5;
        const int ii = tmppos % 4;
        std::fill(check1.begin(), check1.end(), static_cast<T>(-88));
        std::fill(check2.begin(), check2.end(), static_cast<T>(-66));
        guard.run([&](){
          // Will trigger both tweaks, "brick shortcut" has precedence.
          reader->read(size3i_t{ii*bs[0],jj*bs[1],kk*bs[2]}, bs, check1.data());
          // Will trigger only the malloc tweak.
          reader->read(size3i_t{ii*bs[0],jj*bs[1],kk*bs[2]}, bsm1, check2.data());
          std::array<float,3> expect;
          if ((pos % 7) == 0) {
            expect = std::array<float,3>{0,0,0};
          }
          else if ((pos % 13) == 0) {
            T value = static_cast<T>(-(ii+jj+kk));
            expect = std::array<float,3>{(float)value, (float)value, (float)value};
          }
          else {
            expect = std::array<float,3>{(float)ii, (float)jj, (float)kk};
          }
          if (check1[0]!=expect[0]||check1[1]!=expect[1] ||check1[2]!=expect[2]) {
            std::string expect_str = formatMe(pos, ii, jj, kk, expect.data());
            std::string actual_str = formatMe(pos, ii, jj, kk, check1.data());
            TEST_EQUAL(actual_str, expect_str);
            throw std::runtime_error("Mismatch in check 1");
          }
          if (check2[0]!=expect[0]||check2[1]!=expect[1] ||check2[2]!=expect[2]) {
            std::string expect_str = formatMe(pos, ii, jj, kk, expect.data());
            std::string actual_str = formatMe(pos, ii, jj, kk, check2.data());
            TEST_EQUAL(actual_str, expect_str);
            throw std::runtime_error("Mismatch in check 2");
          }
          // Too much hassle to check of the buffer if it is an edge brick.
          // I want to check the others in case the read went across
          // a segment boundary. That triggers soecuak case handling.
          if (ii != 3 && jj != 4 && kk != 5 && (pos%7) != 0 && (pos%13) != 0) {
            for (auto it = check1.begin() + 3; it != check1.end(); ++it)
              if (*it != 42)
                if (!TEST_EQUAL((double)*it, 42))
                  throw std::runtime_error("Mismatch in check 1 last part");
            for (auto it = check2.begin() + 3; it != check2.end() - 1*bsm1[1]*bsm1[2]; ++it)
              if (*it != 42)
                if (!TEST_EQUAL((double)*it, 42))
                  throw std::runtime_error("Mismatch in check 2 last part");
          }
        });
      } // loop
    } // parallel region
    try {
      guard.finished();
    }
    catch (const std::exception& ex) {
      TEST_EQUAL(std::string(ex.what()), std::string("success"));
    }
  }
}

static void
test_bat_local_1()
{
  LocalFileAutoDelete lad("bat1.zgy");
  do_testbat<std::int8_t, 32, false>(lad.name());
}

static void
test_bat_local_2()
{
  LocalFileAutoDelete lad("bat2.zgy");
  do_testbat<std::int16_t, 32, false>(lad.name());
}

static void
test_bat_local_4()
{
  LocalFileAutoDelete lad("bat4.zgy");
  do_testbat<float, 64, false>(lad.name());
}

static void
test_bat_local_zfp()
{
  LocalFileAutoDelete lad("bat.zgy");
  do_testbat<float, 64, true>(lad.name());
}

/**
 * Tests for whether dead traces are included in histogram and statistics.
 *
 * Current behavior, but not critical.
 * - Never-written bricks show up in the histogram and the statistics.
 * - The behavior did not change after the big rewrite.
 * - The behavior does not match the buggy and inconsistent ZGY handling
 * - Dead traces inside a not-dead brick do show up.
 *
 * Required behavior
 * - Padding samples must not be included.
 * - Never-written and const bricks agree on the behavior.
 * - Valuetype and fixed vs. dynamic histogram agree on the behavior.
 * - Clearing the entire survey starts all over, including stats.min/max
 * - Incremental and full rebuild agree on the behavior. (*)
 * - Genlod vs. keep track while writing agree on the behavior. (*)
 *
 * (*) I am considering to remove statistics and histogram generation from
 * genlod and always keep track of changes while writing. That would make
 * the last two checks N/A. It would also make float histograms uglier when
 * doing a full build of compressed file histogram. I.e. as ugly as the
 * build-as-you-go histograms. So, not a big deal.
 */
static void
test_stats_empty()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("stats_empty.zgy");
  const size3i_t   bs{64, 64, 64};
  const size3i_t size{(1*64)+43, (2*64)+34, (6*64)+21};
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(SampleDataType::int8)
    .datarange(-128, +127)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  if (!TEST_CHECK(bool(writer)))
    return;
  writer->finalize();
  writer->close();
  writer.reset();

  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  if (!TEST_CHECK(bool(reader)))
    return;
  SampleStatistics stats = reader->statistics();
  SampleHistogram  histo = reader->histogram();
  TEST_EQUAL(stats.min, 0.0);
  TEST_EQUAL(stats.max, 0.0);
  TEST_EQUAL(stats.cnt, size[0]*size[1]*size[2]);
  TEST_EQUAL_FLOAT(histo.minvalue, -128, 1e-4);
  TEST_EQUAL_FLOAT(histo.maxvalue, +127, 1e-4);
  TEST_EQUAL(histo.samplecount, stats.cnt);
  if (!TEST_EQUAL(histo.bins.size(), 256))
    return;
  TEST_EQUAL(histo.bins[128+0], size[0]*size[1]*size[2]);
}

static void
do_stats_common(const std::string& filename, int testnumber)
{
  typedef IZgyWriter::size3i_t size3i_t;
  const size3i_t   bs{64, 64, 64};
  const size3i_t size{(1*64)+43, (2*64)+34, (6*64)+21};
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(SampleDataType::int8)
    .datarange(-128, +127)
    .filename(filename);
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  if (!TEST_CHECK(bool(writer)))
    return;

  std::vector<std::int8_t> data(bs[0]*bs[1]*bs[2], 0);
  data[0] = 42;
  // Edge brick, non-constant, 43*64*64 samples.
  writer->write(size3i_t{1*64, 0*64, 0*64}, bs, data.data());
  // Full brick explicitly constant = 42, 64*64*64 samples.
  writer->writeconst(size3i_t{0*64, 0*64, 5*64}, bs, data.data());
  // Edge brick explicitly constant = 15,  64*64*21 samples.
  data[0] = 15;
  writer->writeconst(size3i_t{0*64, 0*64, 6*64}, bs, data.data());

  switch (testnumber) {
  case 0:
    break;
  case 1:
    {
      // Overwrite the entire survey with all-constant 1.
      const std::int8_t one{1};
      writer->writeconst(size3i_t{0,0,0}, size, &one);
    }
    break;
  case 2:
    {
      // Splitting the overwrite won't completely reset the file.
      const std::int8_t one{1};
      writer->writeconst(size3i_t{0,0,0}, size3i_t{size[0], size[1], 64}, &one);
      writer->writeconst(size3i_t{0,0,64}, size3i_t{size[0], size[1], size[2]-64}, &one);
    }
    break;
  default:
    TEST_CHECK(false && "unrecognized test number");
  }
  writer->finalize();
  writer->close();
  writer.reset();
}

/**
 * Verify the fix for a particularly obscure bug. When a brick was
 * written as all-constant 1, it would be read back as all-constant 0.
 * Unless that brick had already been allocated on the file.
 * Other constant values worked fine.
 */
static void
test_stats_one()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("stats_one.zgy");
  do_stats_common(lad.name(), 1);

  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  if (!TEST_CHECK(bool(reader)))
    return;
  const size3i_t   bs = reader->bricksize();
  const size3i_t size = reader->size();
  const SampleStatistics stats = reader->statistics();
  const SampleHistogram  histo = reader->histogram();

  for (std::int64_t ii = 0; ii < size[0]; ii += bs[0]) {
    for (std::int64_t jj = 0; jj < size[1]; jj += bs[1]) {
      for (std::int64_t kk = 0; kk < size[2]; kk += bs[2]) {
        std::int8_t value{99};
        reader->read(size3i_t{ii,jj,kk}, size3i_t{1,1,1}, &value, 0);
        TEST_EQUAL((float)value, 1.0f);
      }
    }
  }

  std::pair<bool,double> cvalue =
    reader->readconst(size3i_t{0,0,0}, bs, 0, false);
  if (TEST_CHECK(cvalue.first))
    TEST_EQUAL_FLOAT(cvalue.second, 1.0, 1e-4);

  // This brick will have the correct (1) value, but the brick
  // will not have been deallocated. So the is_const hint will be false.
  //std::pair<bool,double> cvalue2 =
  //  reader->readconst(size3i_t{1*64, 0*64, 0*64}, bs, 0, false);
  // UPDATE: In the new code being developed, there is code to
  // leak any already allocated bricks in this case.
  // TODO-WIP-BrickedAPI: Decide what the code ought to do
  // Decision in ZgyInternalBulk::writeAlignedBrickList().
  //FAILS! TEST_CHECK(!cvalue2.first);

  TEST_EQUAL_FLOAT(stats.min, 1, 1e-6);
  TEST_EQUAL_FLOAT(stats.max, 1, 1e-6);
  TEST_EQUAL(stats.cnt, size[0]*size[1]*size[2]);

  TEST_EQUAL_FLOAT(histo.minvalue, -128, 1e-4);
  TEST_EQUAL_FLOAT(histo.maxvalue, +127, 1e-4);
  TEST_EQUAL(histo.samplecount, stats.cnt);
  if (!TEST_EQUAL(histo.bins.size(), 256))
    return;
  TEST_EQUAL(histo.bins[128+0], 0);
  TEST_EQUAL(histo.bins[128+1], size[0]*size[1]*size[2]);
}

/**
 * If a file has been written to, and then is cleared by a writeconst()
 * of the entire survey, then statistics and histogram must reflect
 * the current contents.
 *
 * The statistics may or may not include the defaultvalue (0).
 * This is a known issue, amd is not considered important.
 */
static void
test_stats_many()
{
  typedef IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("writeconst.zgy");
  do_stats_common(lad.name(), 0);

  std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(lad.name());
  if (!TEST_CHECK(bool(reader)))
    return;
  //const size3i_t   bs = reader->bricksize();
  const size3i_t size = reader->size();
  const SampleStatistics stats = reader->statistics();
  const SampleHistogram  histo = reader->histogram();

  if (verbose()) {
    std::stringstream ss;
    ss << "size=" << size[0]*size[1]*size[2];
    for (std::size_t ii=0; ii<histo.bins.size(); ++ii)
      if (histo.bins[ii] != 0)
        ss << " h[" << ii - 128 << "]=" << histo.bins[ii];
    ss << " stats min " << stats.min << " max " << stats.max << "\n";
    std::cout << "\n" << ss.str() << std::endl;
  }

  TEST_EQUAL_FLOAT(stats.min, 0, 1e-6);
  TEST_EQUAL_FLOAT(stats.max, 42, 1e-6);
  TEST_EQUAL(stats.cnt, size[0]*size[1]*size[2]);

  TEST_EQUAL_FLOAT(histo.minvalue, -128, 1e-4);
  TEST_EQUAL_FLOAT(histo.maxvalue, +127, 1e-4);
  TEST_EQUAL(histo.samplecount, stats.cnt);
  if (!TEST_EQUAL(histo.bins.size(), 256))
    return;
  TEST_EQUAL(histo.bins[128+42], 64*64*64 + 1);
  TEST_EQUAL(histo.bins[128+15], 64*64*21);
  TEST_EQUAL(histo.bins[128+0], histo.samplecount - histo.bins[128+42] - histo.bins[128+15]);
}

#ifdef HAVE_SD

static void
test_bat_sd_1()
{
  CloudFileAutoDelete cad("bat1.zgy", Test_Utils::default_sd_context());
  do_testbat<std::int8_t, 32, false>(cad.name());
}

static void
test_bat_sd_2()
{
  CloudFileAutoDelete cad("bat2.zgy", Test_Utils::default_sd_context());
  do_testbat<std::int16_t, 32, false>(cad.name());
}

static void
test_bat_sd_4()
{
  CloudFileAutoDelete cad("bat4.zgy", Test_Utils::default_sd_context());
  do_testbat<float, 64, false>(cad.name());
}

static void
test_bat_sd_zfp()
{
  CloudFileAutoDelete cad("bat4.zgy", Test_Utils::default_sd_context());
  do_testbat<float, 64, true>(cad.name());
}

static void
showex(const std::exception& /*ex*/)
{
  //if (verbose())
  //  std::cout << "Caught " << typeid(ex).name() << ": " << ex.what() << std::endl;
}

static bool
ok_alturl(const std::string& filename)
{
  std::shared_ptr<IZgyUtils> utils =
    IZgyUtils::utils("sd://", Test_Utils::default_sd_context());
  try {
    std::string alturl = utils->alturl(filename);
    if (alturl.size() < 3*filename.size()) {
      if (verbose())
        std::cout << "No alturl for \"" << filename << "\". Too short.\n";
      return false;
    }
  }
  catch (const OpenZGY::Errors::ZgyNotReadOnlyError& ex)
  {
    showex(ex);
    if (verbose())
      std::cout << "No alturl for \"" << filename << "\". " << ex.what()<< "\n";
    return false;
  }
  if (verbose())
    std::cout << "Good alturl for \"" << filename << "\".\n";
  return true;
}

static bool
ok_writable(const std::string& filename)
{
  try {
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      IZgyWriter::reopen(ZgyWriterArgs()
                         .iocontext(Test_Utils::default_sd_context())
                         .filename(filename));
    if (!TEST_CHECK(bool(writer)))
      return false;
    if (verbose())
      std::cout << "Writable \"" << filename << "\".\n";
    return true;
  }
  catch (const std::exception& ex)
  {
    showex(ex);
    if (verbose())
      std::cout << "Read-only \"" << filename << "\". " << ex.what()<< "\n";
    return false;
  }
}

/**
 * Check handling of the read-only state, given the 3 boolean settings
 * in the IOContext.
 *
 * The following is NOT tested:
 *
 *  - Confirm that opening a file does not set a lock if it is
 *    read-only and vice versa. Problematic to test here because we
 *    have no direct SDAPI access. Also, that test rightly belongs in
 *    SDAPI, not here.
 *
 *  - Investigate whether toggling the read-only flag takes effect
 *    immediately or (more likely) requires a close and re-open. When
 *    a file is open for write and set to read-only mode, can the
 *    writer still output to it? Test is problematic for the same
 *    reasons as above.
 */
static void
test_roflag()
{
  CloudFileAutoDelete cad("roflag.zgy", Test_Utils::default_sd_context());
  SeismicStoreIOContext ctx(*Test_Utils::default_sd_context());

  // Create a read-only file. AltUrl should work, update should not.
  ctx.setRoAfterWrite(true).forceRoBeforeRead(false).forceRwBeforeWrite(false);
  {
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      IZgyWriter::open(ZgyWriterArgs()
                       .iocontext(&ctx)
                       .size(128, 42, 555)
                       .filename(cad.name()));
    if (!TEST_CHECK(bool(writer)))
      return;
    writer->close();
    writer.reset();
    TEST_CHECK(ok_alturl(cad.name()));
    TEST_CHECK(!ok_writable(cad.name()));
  }

  // Create a read/write file. Uses the same file name but this will be
  // a delete and re-create. AltUrl should not work, update should.
  // Assumption: No need to use the forceRwBeforeWrite hack, i.e.
  // SDAPI allows TRUNCATE open even when the file already exists
  // and is read-only. The test also verifies that assumption.
  ctx.setRoAfterWrite(false).forceRoBeforeRead(false).forceRwBeforeWrite(false);
  {
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      IZgyWriter::open(ZgyWriterArgs()
                       .iocontext(&ctx)
                       .size(128, 42, 555)
                       .filename(cad.name()));
    if (!TEST_CHECK(bool(writer)))
      return;
    writer->close();
    writer.reset();
    TEST_CHECK(!ok_alturl(cad.name()));
    TEST_CHECK(ok_writable(cad.name()));
  }

  // The file is now writable. See if we can automatically make it readable
  // when needed because we want to open it without any locking.
  ctx.setRoAfterWrite(false).forceRoBeforeRead(true).forceRwBeforeWrite(false);
  {
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      IZgyReader::open(cad.name(), &ctx);
    if (!TEST_CHECK(bool(reader)))
      return;
    reader->close();
    reader.reset();
    TEST_CHECK(ok_alturl(cad.name()));
    TEST_CHECK(!ok_writable(cad.name()));
  }

  // The file is now readable. See if we can automatically make it writable
  // when the applicaton wants to update it. Leave it writable when done.
  ctx.setRoAfterWrite(false).forceRoBeforeRead(false).forceRwBeforeWrite(true);
  {
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      IZgyWriter::reopen(ZgyWriterArgs()
                         .iocontext(&ctx)
                         .filename(cad.name()));
    if (!TEST_CHECK(bool(writer)))
      return;
    writer->close();
    writer.reset();
    TEST_CHECK(!ok_alturl(cad.name()));
    TEST_CHECK(ok_writable(cad.name()));
  }

  // The file is now writable. See if we can automatically make it readable
  // when needed because we want to get an AltUrl.
  ctx.setRoAfterWrite(false).forceRoBeforeRead(true).forceRwBeforeWrite(false);
  {
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      IZgyReader::open(cad.name(), &ctx);
    if (!TEST_CHECK(bool(reader)))
      return;
    reader->close();
    reader.reset();
    TEST_CHECK(ok_alturl(cad.name()));
    TEST_CHECK(!ok_writable(cad.name()));
  }
}

/**
 * Follow up on api.nolock, see if SDAPI sets locks even though we
 * sabotage as much as possible of that feature.
 */
static void
test_roflag_nolock()
{
  CloudFileAutoDelete cad("roflag_nolock.zgy", Test_Utils::default_sd_context());
  SeismicStoreIOContext ctx(*Test_Utils::default_sd_context());
  ctx.setRoAfterWrite(true).forceRoBeforeRead(true).forceRwBeforeWrite(true);
  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  std::shared_ptr<OpenZGY::IZgyReader> reader;

  // Create a new file and mark it readonly.
  writer = IZgyWriter::open(ZgyWriterArgs()
                            .filename(cad.name())
                            .iocontext(&ctx)
                            .size(128, 42, 555));
  if (!TEST_CHECK(bool(writer)))
    return;
  writer->close();
  writer.reset();

  // Open the newly created file. Had it not been readonly already,
  // it would have been closed, marked readonly, and opened again.
  reader = IZgyReader::open(cad.name(), &ctx);
  if (!TEST_CHECK(bool(reader)))
    return;

  // Try to trick SDAPI to allow opening the file for write, which is
  // supposed to require an exclusive lock. But the read lock was not
  // set since the file was supposedly readonly. Which also means that
  // we should have been able to open the readonly file, mark it as
  // writable, and then open it for update.
  writer = IZgyWriter::reopen(ZgyWriterArgs()
                              .filename(cad.name())
                              .iocontext(&ctx)
                              );
  if (!TEST_CHECK(bool(writer)))
    return;

  // What happens if the file gets marked readonly, perhaps by
  // reader->close, while still being written to? or by the writer
  // itself, right after opening? Possibly academic since neither of
  // the three "force" flags will allow this. Had it been possible, my
  // guess is that the already acquired locks would remain.

  writer->close();
  reader->close();
}

/**
 * Follow up on api.nolock, see if SDAPI will immediately drop all
 * acquired read locks when the readonly flag is set. I guess it is more likely that readonly (i.e. don't set any read locks) only takes effect for new files.
 */
static void
test_roflag_unlock()
{
  CloudFileAutoDelete cad("roflag_unlock.zgy", Test_Utils::default_sd_context());
  SeismicStoreIOContext ctx(*Test_Utils::default_sd_context());
  ctx.setRoAfterWrite(false).forceRoBeforeRead(false).forceRwBeforeWrite(false);

  // Create a new file and leave it writable.
  auto writer = IZgyWriter::open(ZgyWriterArgs()
                                 .filename(cad.name())
                                 .iocontext(&ctx)
                                 .size(128, 42, 555));
  if (!TEST_CHECK(bool(writer)))
    return;
  writer->close();
  writer.reset();

  // Open the newly created file for read. The file itself is
  // writable, so this will set a read a.k.a shared lock.
  auto reader1 = IZgyReader::open(cad.name(), &ctx);
  if (!TEST_CHECK(bool(reader1)))
    return;

  // Mark the file readonly, then open it in a separate handle.
  ctx.setRoAfterWrite(true).forceRoBeforeRead(true).forceRwBeforeWrite(true);
  auto reader2 = IZgyReader::open(cad.name(), &ctx);
  if (!TEST_CHECK(bool(reader2)))
    return;

  // Allow time to run "sdutil stat" from a terminal window.
  //std::cerr << "\nRUN SDUTIL NOW\n" << std::flush;
  //std::this_thread::sleep_for(std::chrono::seconds(42));

  // Did both read locks get dropped? It should then be possible to
  // open for write after setting the file back to read/write.
  must_throw("is locked for read", [&](){
      writer = IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(cad.name())
                                  .iocontext(&ctx)
                                  );});
  if (writer)
    writer->close();
  reader2->close();
  reader1->close();
}

/**
 * On the build server the normal case for access tokens is that the
 * client credentials grant "$CLIENT_ID" and "$CLIENT_SECRET" is used
 * to obtain a token outside of the actual tests. See private/grabtoken.sh.
 * Now OpenZGY can handle client credentials grant directly, *if* SDAPI was
 * built with the SLB-internal extensions enabled. Test this feature.
 */
static void
test_client_cred()
{
  using InternalZGY::Environment;
  std::string name = cloud_synt2_name();
  std::string client_id = Environment::getStringEnv("CLIENT_ID");
  std::string client_secret = Environment::getStringEnv("CLIENT_SECRET");
  std::string audience = Environment::getStringEnv("SD_TARGET_AUDIENCE");
  std::string token_service = Environment::getStringEnv("SDAPI_SAUTH_TOKEN_SERVICE_URL");
  if (client_id.empty() || client_secret.empty()) {
    // This is normal when running outside the build servers.
    if (verbose())
      std::cerr << "Skipping test. No client credentials.\n";
  }
  else {
    // The $SD_TARGET_AUDIENCE is picked up inside SDAPI.
    // Unlike $CLIENT_ID and $CLIENT_SECRET it cannot be
    // specified explicitly. If unset the open will fail
    // with an obscure error message, so test it here.
    TEST_CHECK(!audience.empty());
    TEST_CHECK(!token_service.empty());
    SeismicStoreIOContext ctxt = SeismicStoreIOContext()
      .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
      .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
      .sdtoken("SlbAuthServiceAccount:" + client_id + ":" + client_secret);
    std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(name, &ctxt);
    if (!TEST_CHECK(bool(reader)))
      return;
    reader->close();
    if (verbose())
      std::cerr << "Success in using client credentials grant.\n";
  }
}

#if 0 // Only for ad-hoc debugging

namespace {
SeismicStoreIOContext getContext()
{
  using InternalZGY::Environment;
  return SeismicStoreIOContext()
    .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
    .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
    .sdtoken(Environment::getStringEnv("OPENZGY_TOKEN") != "" ?
             Environment::getStringEnv("OPENZGY_TOKEN") :
             "FILE:carbon.slbapp.com", "");
}
}

/**
 * Not usable as an automated test.
 * Partly because it takes two hours to run, waiting for the access
 * token to expire, and partly because somebody needs to check by hand
 * that garbage data was not output to the console.
 */
static void test_bug_671969()
{
  auto ctxt = getContext();
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open("sd://opendes/slb-testdata/Synt2.zgy", &ctxt);
  const std::array<std::int64_t,3> orig{0,0,0};
  const std::array<std::int64_t,3> more{64,192,64};
  const std::array<std::int64_t,3> size{64,64,64};
  std::unique_ptr<float[]>buf(new float[64*64*64]);
  try {
    std::cerr << "1st read" << std::endl;
    reader->read(orig, size, buf.get(), 0);
    for (int ii=10; ii>0; --ii) {
      std::cerr << ii << "... " << std::flush;
      std::this_thread::sleep_for (std::chrono::minutes(12));
    }
    std::cerr << "0!" << std::endl;
    std::cerr << "2nd read" << std::endl;
    reader->read(more, size, buf.get(), 0);
    std::cerr << "close" << std::endl;
    reader->close();
    std::cerr << "done" << std::endl;
  }
  catch(const std::exception& ex) {
    std::cerr << "Exception: " << ex.what() << std::endl;
  }
}

#endif // Ad-hoc debugging

/**
 * Not usable as an automated test.
 */
static void test_list_segments()
{
  const std::string filename = InternalZGY::Environment::getStringEnv
    ("FILENAME", "sd://opendes/slb-testdata/Synt2.zgy");

  std::cout << "\n\nUsing iterator over  \"" << filename << "\"\n";
  auto list = Test_Utils::get_segments(filename);
  std::sort(list.begin(), list.end());
  for (const auto& it : list)
    if ((it.second % (1024*1024)) != 0)
      std::cout << std::setw(16) << it.first << "  " << it.second << "\n";
    else
      std::cout << std::setw(16) << it.first << "  " << it.second/(1024*1024) << " MB\n";

  std::cout << "\nAssuming linear layout \"" << filename << "\"\n";
  std::vector<std::int64_t> sizes = Test_Utils::get_segsizes(filename);
  for (const auto& it : sizes)
    std::cout << it << "\n";
}

/**
 * No-op when used as an automated test.
 * Set $INAME and $ONAME to make an actual copy.
 * When $ONAME is a folder, it must already exist and be empty.
 * Yes, all this is a major hack. "sdutil cp" ought to have
 * supported upload and download with folders.
 */
static void test_cp()
{
  using InternalZGY::Environment;
  const std::string iname = Environment::getStringEnv("INAME");
  const std::string oname = Environment::getStringEnv("ONAME");
  const bool sd_src = iname.substr(0, 5) == "sd://";
  const bool sd_dst = oname.substr(0, 5) == "sd://";
  if (sd_src && oname.empty()) {
    if (verbose())
      std::cout << "\nList \"" << iname << "\"\n";
    for (const auto& it : Test_Utils::get_segments(iname)) {
      if (verbose())
        std::cout << "Item " << std::setw(10) << it.second
                  << " \"" << it.first << "\"" << std::endl;
    }
    //const std::uint64_t nblocks = src.getBlockNum();

  }
  else if (iname.empty() || oname.empty()) {
    if (verbose())
      std::cout << "\nSet INAME and optionally ONAME to run api.cp\n";
  }
  else if (sd_src && !sd_dst)
    Test_Utils::copy_sd_to_folder(iname, oname, verbose());
  else if (!sd_src && sd_dst)
    Test_Utils::copy_folder_to_sd(iname, oname, verbose());
  else if (sd_src && sd_dst)
    throw std::runtime_error("sd-to-sd copies not implemented yet");
  else if (!sd_src && !sd_dst)
    throw std::runtime_error("local copies not implemented yet");
}

#endif

namespace {
  static void
  printBuilder(const char *msg, const SampleHistogram& h, const SampleStatistics& s, bool details)
  {
    const double width = h.bins.size() <= 1 ? 0 : (h.maxvalue - h.minvalue) / (h.bins.size() - 1);
    const double begin = h.minvalue;
    printf("%sHisto range %+g..%+g in %d bins width %g, stats count %d inf N/A range %+g..%+g\n",
           msg, h.minvalue, h.maxvalue, (int)h.bins.size(), width,
           (int)s.cnt, (float)s.min, (float)s.max);
    if (details) {
      for (std::size_t ii=0; ii<h.bins.size(); ++ii) {
        const double center = begin+ii*width;
        const double lo = begin+(ii-0.5)*width;
        const double hi = begin+(ii+0.5)*width;
        // Print first, last, value==0, and all non-empty bins.
        if (h.bins[ii] != 0 || lo*hi <= 0 || ii == 0 || ii == h.bins.size()-1)
          printf("  bin[%d]: %d value %+g (%+g..%+g)\n", (int)ii, (int)h.bins[ii], center, lo, hi);
      }
    }
  }
}

/**
 * Test the public wrapper for the histogram Builder.
 *
 * The focus is on just the API. White-box testing of the functionality
 * is done in test_histobuilder.cpp. See also dynamicExample() and
 * StaticExample() in test_histobuilder.cpp.
 *
 * Collect sample statistics and create a 4096-bin intermediate histogram.
 * When done collecting, convert the histogram to 256-bin zero-centric.
 * The intermediate histogram can also be accessed. But that will have many
 * empty bins and won't be zero-centric.
 *
 * When the builder is not fixed-width, the intermediate size a.k.a. bin count
 * must be at least twice that of the desired result. Choosing a too large
 * size might have a performance impact. Especially on the cost of operator+=.
 *
 * For multi-threaded computation, use one builder for each thread and combine
 * the results at the end by using operator+=.
 */
static void
test_build_hist_dynamic()
{
  const std::vector<float> data1 {-99, 56, 23};
  const std::vector<float> data2 {42};
  const std::vector<float> data3 {0};

  // Add samples from two very short traces.
  SampleHistogramBuilder builder(4096);
  builder.add(data1.data(), data1.data() + data1.size());
  builder.add(data2.data(), data2.data() + data2.size());

  // Add 1001 samples with value zero.
  SampleHistogramBuilder tmp(4096);
  tmp.add(data3.data(), data3.data() + data3.size());
  tmp *= 1001;
  builder += tmp;

  // And subtract two of them.
  SampleHistogramBuilder tmp2(4096);
  tmp2.add(data3.data(), data3.data() + data3.size());
  tmp2 *= 2;
  builder -= tmp2;

  const SampleStatistics& stats = builder.getstats();
  const SampleHistogram&  histo = builder.finalize(256);

  if (verbose())
    printBuilder("\n", histo, stats, true);

  TEST_EQUAL(stats.cnt, 1003);
  TEST_EQUAL_FLOAT(stats.min, -99, 1e-7);
  TEST_EQUAL_FLOAT(stats.max, +56, 1e-7);
  if (!TEST_EQUAL(histo.bins.size(), 256)) return;
  TEST_EQUAL(histo.samplecount, 1003);
  TEST_EQUAL_FLOAT(histo.minvalue, -100.625, 0.005);
  TEST_EQUAL_FLOAT(histo.maxvalue,  +58.750, 0.005);
  TEST_EQUAL(histo.bins[0], 0);     // -100.938..-100.313
  TEST_EQUAL(histo.bins[3], 1);     // -99.0625..-98.4375
  TEST_EQUAL(histo.bins[161], 999); // -0.31250..+0.31250
  TEST_EQUAL(histo.bins[198], 1);   // +22.8125..+23.4375
  TEST_EQUAL(histo.bins[228], 1);   // +41.5625..+42.1875
  TEST_EQUAL(histo.bins[251], 1);   // +55.9375..+56.5625
  TEST_EQUAL(histo.bins[255], 0);   // +58.4375..+59.0625
}

/**
 * Test the public wrapper for the histogram Builder.
 *
 * The focus is on just the API. White-box testing of the functionality
 * is done in test_histobuilder.cpp. See also dynamicExample() and
 * StaticExample() in test_histobuilder.cpp.
 *
 * Generating a fixed-width histogram, limits given by caller.
 * The limits chosen here will not create a zero-centric histogram.
 * Making a better choice is left as an exercise for the student.
 */
static void
test_build_hist_static()
{
  const std::vector<float> data1 {-99, 56, 23};
  const std::vector<float> data2 {42};
  const std::vector<float> data3 {0};

  // Add samples from two very short traces.
  SampleHistogramBuilder builder(256, -100, +100);
  builder.add(data1.data(), data1.data() + data1.size());
  builder.add(data2.data(), data2.data() + data2.size());

  // Add 1001 samples with value zero.
  SampleHistogramBuilder tmp(256, -100, +100);
  tmp.add(data3.data(), data3.data() + data3.size());
  tmp *= 1001;
  builder += tmp;

  // And subtract two of them.
  SampleHistogramBuilder tmp2(256, -100, +100);
  tmp2.add(data3.data(), data3.data() + data3.size());
  tmp2 *= 2;
  builder -= tmp2;

  const SampleStatistics& stats = builder.getstats();
  const SampleHistogram&  histo = builder.finalize(256);

  if (verbose())
    printBuilder("\n", histo, stats, true);

  TEST_EQUAL(stats.cnt, 1003);
  TEST_EQUAL_FLOAT(stats.min, -99, 1e-7);
  TEST_EQUAL_FLOAT(stats.max, +56, 1e-7);
  if (!TEST_EQUAL(histo.bins.size(), 256)) return;
  TEST_EQUAL(histo.samplecount, 1003);
  TEST_EQUAL_FLOAT(histo.minvalue, -100, 1e-7);
  TEST_EQUAL_FLOAT(histo.maxvalue, +100, 1e-7);
  TEST_EQUAL(histo.bins[0], 0);     // -100.392..-99.6078
  TEST_EQUAL(histo.bins[1], 1);     // -99.6078..-98.8235
  TEST_EQUAL(histo.bins[127], 999); // -0.78431..+0.00000
  TEST_EQUAL(histo.bins[128], 0);   // +0.00000..+0.78431
  TEST_EQUAL(histo.bins[157], 1);   // +22.7451..+23.5294
  TEST_EQUAL(histo.bins[181], 1);   // +41.5686..+42.3529
  TEST_EQUAL(histo.bins[199], 1);   // +55.6863..+56.4706
  TEST_EQUAL(histo.bins[255], 0);   // +99.6078..+100.392
}


static void
test_build_hist_with_vt()
{
  // A histogram isn't tied to the storage valuetype,
  // and it doesn't matter which type is used to fill it.
  const std::vector<float>        data1 {99, 56, 23};
  const std::vector<std::int16_t> data2 {101, 109};
  const std::vector<std::int8_t>  data3 {5, 55};
  const std::vector<double>       alldata {99, 56, 23, 101, 109, 5, 55};

  SampleHistogramBuilder builder(256, 0, 255);
  builder.add(data1.data(), data1.data() + data1.size());
  builder.add(data2.data(), data2.data() + data2.size());
  builder.add(data3.data(), data3.data() + data3.size());

  const SampleStatistics& stats = builder.getstats();
  const SampleHistogram&  histo = builder.gethisto();

  TEST_EQUAL(stats.cnt, 7);
  TEST_EQUAL_FLOAT(stats.min, 5, 1e-7);
  TEST_EQUAL_FLOAT(stats.max, 109, 1e-7);
  TEST_EQUAL_FLOAT(stats.sum, 448, 1e-7);
  if (!TEST_EQUAL(histo.bins.size(), 256)) return;
  TEST_EQUAL(histo.samplecount, 7);
  TEST_EQUAL_FLOAT(histo.minvalue, 0, 1e-7);
  TEST_EQUAL_FLOAT(histo.maxvalue, 255, 1e-7);
  for (auto sample : alldata) {
    TEST_EQUAL(histo.bins[static_cast<std::int64_t>(sample)], 1);
  }

  builder.scale(0, 255, 0, 510);
  const SampleStatistics& stats2 = builder.getstats();
  const SampleHistogram&  histo2 = builder.gethisto();

  TEST_EQUAL(stats2.cnt, 7);
  TEST_EQUAL_FLOAT(stats2.min, 2*5, 1e-7);
  TEST_EQUAL_FLOAT(stats2.max, 2*109, 1e-7);
  TEST_EQUAL_FLOAT(stats2.sum, 2*448, 1e-7);
  if (!TEST_EQUAL(histo2.bins.size(), 256)) return;
  TEST_EQUAL(histo2.samplecount, 7);
  TEST_EQUAL_FLOAT(histo2.minvalue, 2*0, 1e-7);
  TEST_EQUAL_FLOAT(histo2.maxvalue, 2*255, 1e-7);
}

/**
 * Return the number of lods available when reading.
 * Also do some consistency checks.
 */
static int
apparent_nlods(const std::string& filename)
{
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(filename);
  const std::int32_t nlods = reader->nlods();
  const bool singlebrick = (reader->size()[0] <= reader->bricksize()[0] &&
                            reader->size()[1] <= reader->bricksize()[1] &&
                            reader->size()[2] <= reader->bricksize()[2]);
  // Consistency check: missing lod means reading lod>0 is invalid
  // and (if there are multiple bricks) the old library cannot read this.
  const std::int64_t version = reader->filestats()->fileVersion();
  bool can_read_lods = false;
  try {
    float data{-999.25};
    reader->read(std::array<std::int64_t,3>{0,0,0},
                 std::array<std::int64_t,3>{1,1,1},
                 &data, /*lod=*/1);
    can_read_lods = true;
  }
  catch (const OpenZGY::Errors::ZgyUserError&) {
    can_read_lods = false;
  }
  reader->close();
  bool ok{true};
  if (singlebrick) {
    ok = TEST_CHECK(!can_read_lods && "No lods in single-brick cube") && ok;
    ok = TEST_CHECK(version == 3 && "No lods, so old reader works") && ok;
  }
  else if (nlods == 1)
  {
    ok = TEST_CHECK(!can_read_lods && "Cannot read lod ?") && ok;
    ok = TEST_CHECK(version == 4 && "Old reader requires LOD to exist") && ok;
  }
  else {
    ok = TEST_CHECK(can_read_lods && "Why can't the code read any lods?") && ok;
    ok = TEST_CHECK(version == 3 && "Old reader can read normal files") && ok;
  }
  return ok ? nlods : -1;
}

/**
 * Test the check for incomplete lowres.
 *
 * The focus is only on whether fullres is available or not when
 * reading a file. Test added because the rules are about to change.
 * Also, available lowres usually implies that the old reader can open
 * this file, and vice versa. Compression would also exclude the old
 * accessor, but that is not tested here. There is also an open
 * question whether it is valid for the application to read lowres
 * bricks at their own risk even when nlods() claim there aren't any.
 *
 * See also the reopen.seversion test. There is considerble overlap.
 *
 * The old accessor reqires all LOD bricks to be correct.
 *
 *  - single brick empty or all-const with no lod (zero data bricks).  2 tests.
 *  - single brick allocated with no lod (one data brick).             1 test.
 *  - empty or all-const all lods (any survey size, zero real bricks). 2 tests.
 *  - [logically all-const but some inflated bricks in any lod.]       0 tests.
 *  - [logically all-const in lod0 but stale lowres exists.]           0 tests.
 *  - file with non-const data, no lowres at all (no lowres bricks).   1 test.
 *  - file with non-const data, partial lowres but not the top brick.  1 test.
 *  - file with non-const data, fully generated incr or at end.        2 tests.
 */
static void
test_nlods()
{
  LocalFileAutoDelete lad("nlods.zgy");
  const std::array<std::int64_t,3> zero{0,0,0};
  const std::array<std::int64_t,3> onesample{1,1,1};
  const float fortytwo{42};
  const ZgyWriterArgsV2 small = ZgyWriterArgsV2()
    .filename(lad.name())
    .decimation(std::vector<DecimationType>{DecimationType::Average})
    .size(12, 13, 14); // 1 lod, 1 brick.
  const ZgyWriterArgsV2 large = ZgyWriterArgsV2()
    .filename(lad.name())
    .decimation(std::vector<DecimationType>{DecimationType::Average})
    .size(100, 200, 300); // 4 lods, 2*4*5 = 40 lod0 bricks, 6 lod1, 2 lod2, 1.

  {
    // Single brick empty with no lod (zero real data bricks).
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(small));
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 1);
  }

  {
    // Single brick all-const with no lod (zero real data bricks).
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(small));
    writer->writeconst(zero, writer->size(), &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 1);
  }

  {
    // Single brick allocated with no lod (one real data brick).
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(small));
    writer->write(zero, onesample, &fortytwo);
    TEST_EQUAL(apparent_nlods(lad.name()), 1);
  }

  {
    // Empty all lods (any survey size, zero allocated data bricks).
    // Even though lods were not generated, they will appear to have been.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large).lodmode(LodMode::Never));
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 4);
  }

  {
    // All-const all lods (any survey size, zero allocated data bricks).
    // Even though lods were not generated, they will appear to have been.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large).lodmode(LodMode::Never));
    writer->writeconst(zero, writer->size(), &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 4);
  }

  {
    // File with non-const data, no lowres at all (no lowres bricks allocated).
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large).lodmode(LodMode::Never));
    writer->write(zero, onesample, &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 1);
  }

  {
    // File with non-const data, partial lowres but not the top brick.
    // Note that the ability to do this might be removed from the API.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large)
                                .lodIncrForce(1)
                                .lodIncrLevel(1)
                                .lodLastForce(2)
                                .lodLastLevel(1));
    writer->write(zero, onesample, &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 1);
  }

  {
    // File with non-const data, fully generated lowres.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large).lodmode(LodMode::Early));
    writer->write(zero, onesample, &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 4);
  }

  {
    // File with non-const data, fully generated lowres (at end).
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgsV2(large).lodmode(LodMode::Rebuild));
    writer->write(zero, onesample, &fortytwo);
    writer->close();
    TEST_EQUAL(apparent_nlods(lad.name()), 4);
  }
}

namespace TestIncrLod {
#if 0
}
#endif
static const std::vector<std::int64_t> contents
  {    3,   5,   7,  11,  13,  17,  19,  23,  29,
      31,  37,  41,  43,  47,  53,  59,  61,  67,
      71,  73,  79,  83,  89,  97, 101, 103, 107,
     109, 113, 127, 131, 137, 139, 149, 151, 157,
     163, 167, 173, 179, 181, 191, 193, 197, 199,
  };

// Describes one brick-column in the test.
struct PosAndValue
{
  std::array<std::int64_t,3> start;
  std::array<std::int64_t,3> count;
  std::int16_t sample;
  PosAndValue(int i0, int j0, int nk, std::int64_t val)
    : start(std::array<std::int64_t,3>{i0*64, j0*64,     0})
    , count(std::array<std::int64_t,3>{ 1*64,  1*64, nk*64})
    , sample((std::int16_t)val)
  {
  }
};

static std::vector<PosAndValue>
getWork(int ni, int nj, int nk)
{
  std::vector<PosAndValue> result;
  int c{0};
  const int wrap{(int)contents.size()};
  for (int ii=0; ii<ni; ++ii)
    for (int jj=0; jj<nj; ++jj, ++c)
      result.push_back(PosAndValue(ii, jj, nk, 100*contents[c % wrap]));
  return result;
}

/*
 * The "kludge" argument can be patched to "true" to adjust expected
 * results assuming that stat/hist was not collected from the first
 * write in api.incrlod3. Also patch newTrackedBricksTryEnable to
 * pretend that the file being re-opened did not contain stats/histo.
 * Then run api.incrlod3. Could make a similar kludge for api.incrlod2.
 * Left as an exercise for the reader.
 */
static bool
check_incrlod(const std::string& filename, bool kludge = false)
{
  const std::vector<PosAndValue> work = getWork(9, 5, 3);
  std::shared_ptr<IZgyReader> reader = IZgyReader::open(filename);
  bool ok{true};

  // The unit test is about low resolution bricks, but it shouldn't
  // hurt to check statistics and histogram as well.
  const SampleStatistics stat = reader->statistics();
  const SampleHistogram  hist = reader->histogram();

  ok = TEST_EQUAL(stat.cnt, (9*5*(kludge?2:3)) * (64*64*64)) && ok;
  ok = TEST_EQUAL(stat.max, 100*contents.back()) && ok;
  ok = TEST_EQUAL(hist.samplecount, stat.cnt) && ok;
  ok = TEST_EQUAL(hist.minvalue, -32768) && ok;
  ok = TEST_EQUAL(hist.maxvalue, +32767) && ok;

  // Check the histogram itself. Assume no two sample values share a bin.
  auto getcount = [&hist](int value) -> std::int64_t
    {
      if (hist.maxvalue <= hist.minvalue)
        return 0;
      const double B = (hist.bins.size() - 1)/(hist.maxvalue - hist.minvalue);
      const double A = -hist.minvalue * B;
      const int slot = (int)std::floor(A + B*value + 0.5);
      if (slot < 0 || slot >= (int)hist.bins.size())
        return 0;
      //std::cout << "value " << value
      //          << " slot " << slot
      //          << " count " << hist.bins[slot] << std::endl;
      return hist.bins[slot];
    };
  for (const PosAndValue& it : work) {
    if (it.sample <= 500) // two sample values share the same bin
      ok = TEST_EQUAL(getcount(it.sample), (kludge?4:6)*64*64*64) && ok;
    else
      ok = TEST_EQUAL(getcount(it.sample), (kludge?2:3)*64*64*64) && ok;
  }

  ok = TEST_EQUAL(reader->nlods(), 5) && ok;
  ok = TEST_EQUAL(work.size(), 45) && ok;
  for (std::int32_t lod = 0; lod < reader->nlods(); ++lod) {
    std::int64_t factor = std::int64_t(1) << lod;
    for (const PosAndValue& it : work) {
      std::array<std::int64_t,3> pos
        {it.start[0] / factor,
         it.start[1] / factor,
         it.start[2] / factor};
      std::array<std::int64_t,3> siz{1, 1, 1};
      std::int16_t sample{0};
      reader->read(pos, siz, &sample, lod);
      ok = TEST_EQUAL(sample, it.sample) && ok;
    }
  }
  return ok;
}

static std::shared_ptr<IZgyWriter>
open_incrlod(const std::string& filename)
{
  return IZgyWriter::open
    (ZgyWriterArgsV2()
     .filename(filename)
     .size(9*64, 5*64, 3*64)
     .datatype(SampleDataType::int16)
     .datarange(-32768, +32767)
     .decimation(std::vector<DecimationType>{DecimationType::Average})
     );
}

static std::shared_ptr<IZgyWriter>
reopen_incrlod(const std::string& filename)
{
  return IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .filename(filename)
     .decimation(std::vector<DecimationType>{DecimationType::Average})
     );
}
} // namespace TestIncrLod

/**
 * Test that low resolution bricks get correctly written even if the
 * file was closed and re-opened during the writing part and
 * regardless of whether read/modify/wride was needed due to
 * misaligned writes. As a secondary test, also check that the
 * statistics and histogram are correct.
 *
 * api.incrlod1 writes in a single session with no r/m/w.
 * api.incrlod2 writes in two sessions with no r/m/w.
 * api.incrlod3 writes in two sessions (for upper/lower) and with r/m/w.
 *
 * TODO-WIP-BrickedAPI: Possible additional checks: Measure the I/O
 * during finalize and (more difficult) during incremental write.
 * Case 1 should write each brick just once. Case 2 more, case 3 most.
 */
static void
test_incrlod1()
{
  using namespace TestIncrLod;
  LocalFileAutoDelete lad("incrlod1.zgy");
  std::shared_ptr<IZgyWriter> writer = open_incrlod(lad.name());
  const std::vector<PosAndValue> work = getWork(9, 5, 3);
  for (const PosAndValue& it : work)
    writer->writeconst(it.start, it.count, &it.sample);
  writer->finalize(); // TODO minor, count the bricks handled here.
  writer->close();
  TEST_CHECK(check_incrlod(lad.name()));
}

static void
test_incrlod2()
{
  using namespace TestIncrLod;
  LocalFileAutoDelete lad("incrlod2.zgy");
  std::shared_ptr<IZgyWriter> writer = open_incrlod(lad.name());
  const std::vector<PosAndValue> work = getWork(9, 5, 3);
  std::vector<PosAndValue> work1(work.begin(), work.begin() + 32);
  std::vector<PosAndValue> work2(work.begin() + 32, work.end());
  for (const PosAndValue& it : work1)
    writer->writeconst(it.start, it.count, &it.sample);
  writer->finalize();
  writer->close();
  writer = reopen_incrlod(lad.name());
  for (const PosAndValue& it : work2)
    writer->writeconst(it.start, it.count, &it.sample);
  writer->finalize();
  writer->close();
  TEST_CHECK(check_incrlod(lad.name()));
}

static void
test_incrlod3()
{
  using namespace TestIncrLod;
  LocalFileAutoDelete lad("incrlod3.zgy");
  std::shared_ptr<IZgyWriter> writer = open_incrlod(lad.name());
  const std::vector<PosAndValue> work = getWork(9, 5, 3);
  for (PosAndValue it : work) {
    // "it" is deliberately made a mutable, temporary copy.
    it.count[2] = std::min((std::int64_t)70, it.count[2]);
    writer->writeconst(it.start, it.count, &it.sample);
  }
  writer->finalize();
  writer->close();
  writer = reopen_incrlod(lad.name());
  for (PosAndValue it : work) {
    it.start[2] = 42;
    it.count[2] -= 42;
    writer->writeconst(it.start, it.count, &it.sample);
  }
  writer->finalize();
  writer->close();
  TEST_CHECK(check_incrlod(lad.name()));
}

/**
 * A single-brick file opened for update can use any LodMode,
 * regardless of whether it contains data or not. Check with/without
 * data and check modes Never and Early.
 */
static void
test_lodmode1()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode1.zgy");

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(22, 33, 44)
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  w->close();

  // First Never, Second Never.
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  w->close();

  // All previous Never, now Early.
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Early)
     .filename(lad.name()));
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &fortytwo);
  w->close();

  // Previous Early, now Never. Still ok because no actual
  // low resolution data was generated.
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  w->close();
}

/**
 * A file opened for update that still has no data can be reopened in
 * any mode. Regardless of whether the previous mode was Never, Early,
 * or Rebuild. "No data" can be either never written or reset to
 * constant.
 */
static void
test_lodmode2()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode2.zgy");

  const std::vector<LodMode> modes{LodMode::Never, LodMode::Rebuild};

  for (const LodMode mode : modes) {
    std::shared_ptr<IZgyWriter> w = IZgyWriter::open
      (ZgyWriterArgsV2()
       .size(64, 64, 512)
       .lodmode(mode)
       .filename(lad.name()));
    w->close();

    w = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .lodmode(LodMode::Never)
       .filename(lad.name()));
    w->close();

    w = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .lodmode(LodMode::Rebuild)
       .filename(lad.name()));
    const float sixteen{16};
    w->writeconst(index3_t{0,0,0}, w->size(), &sixteen);
    w->close();

    w = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .lodmode(LodMode::Never)
       .filename(lad.name()));
    w->close();

    w = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .lodmode(LodMode::Rebuild)
       .filename(lad.name()));
    w->close();
  }
}

/**
 * A file opened for update that has data and lowres cannot use
 * LodMode::Never, but it can use Early. Or Rebuild, but that is
 * unnecessary to test.
 */
static void
test_lodmode3()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode3.zgy");

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(64, 64, 512)
     .lodmode(LodMode::Early)
     .filename(lad.name()));
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &fortytwo);
  w->close();

  must_throw("LodMode::Never not allowed on a file that", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .filename(lad.name()));
      w->close();
    });

  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Early)
     .filename(lad.name()));
  w->close();
}

/**
 * A file with real data written with LodMode::Never can only be
 * reopened with Never (keep the file with no lowres) or Rebuild.
 */
static void
test_lodmode4()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode4.zgy");

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(64, 64, 512)
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &fortytwo);
  w->close();

  const std::vector<LodMode> badmodes{LodMode::Early, LodMode::Early1, LodMode::Default};
  for (const LodMode mode : badmodes) {
    must_throw("LodMode::Rebuild or Never needed when the file has no lowres.", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(mode)
         .filename(lad.name()));
      w->close();
    });
  }

  const std::vector<LodMode> okmodes{LodMode::Never, LodMode::Rebuild};
  for (const LodMode mode : okmodes) {
    w = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .lodmode(mode)
       .filename(lad.name()));
    w->close();
  }
}

/**
 * A file with real data and real lowres, but with a decimation
 * algorithm that ended up with only zeros in the higer lods, might
 * behave oddly. Try to run some specific tests.
 */
static void
test_lodmode5()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode5.zgy");
  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(64, 64, 512)
     .lodmode(LodMode::Rebuild)
     .decimation(std::vector<DecimationType>{DecimationType::Decimate})
     .filename(lad.name()));
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,2}, index3_t{1,1,1}, &fortytwo);
  // The trace in level 0: {0, 0, 42, 0, ...}
  // Decimated in level 1, every other sample: {0, 42, 0, ...}
  // Decimated in lever 2, {0, 0, ...}
  w->close();

  // Should not be allowed to reopen in LodMode::Never, but it may look
  // like we have no lowres.
  must_throw("LodMode::Never not allowed on a file that already has lowres.", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .filename(lad.name()));
      w->close();
    });

  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Early1)
     .filename(lad.name()));
  w->close();
}

/**
 * A file that used to contain real data, with or without lowres, that
 * was subsequently cleared, might behave oddly because the checks
 * don't account for inflated bricks i.e. data written to disk that
 * ought to have been stored as a const brick.
 *
 * Currently, clearing the survey in a single call will delberately
 * leak any existing bricks and result in a truly cleared survey.
 *
 * Instead the test will clear the survey in two parts; this leaves
 * inflated bricks and also the statistics will be too wide.
 */
static void
test_lodmode6()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode6.zgy");

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(64, 64, 512)
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  const float sixteen{16};
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &fortytwo);
  //w->writeconst(index3_t{0,0,0}, w->size(), &sixteen);
  w->writeconst(index3_t{0,0,0}, index3_t{64,64,64}, &sixteen);
  w->writeconst(index3_t{0,0,64}, index3_t{64,64,512-64}, &sixteen);
  w->close();

  // The file now has inflated lod0 bricks, no lowres bricks. This
  // will look like real data was written in LodMode::Never, which is
  // sort of was. Looking at the statistics would show that the file
  // was actually cleared. Early will not be allowed, while Never
  // will be.

  std::shared_ptr<IZgyReader> r = IZgyReader::open(lad.name());
  TEST_EQUAL(r->statistics().min, 0); // 16 if clearing in one step.
  TEST_EQUAL(r->statistics().max, 42); // 16 if clearing in one step.
  r->close();

  must_throw("LodMode::Rebuild or Never needed when the file has no lowres.", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Early)
         .filename(lad.name()));
      w->close();
    });

  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Never)
     .filename(lad.name()));
  w->close();

  w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(64, 64, 512)
     .lodmode(LodMode::Rebuild)
     .filename(lad.name()));
  w->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &fortytwo);
  w->close();
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Rebuild)
     .filename(lad.name()));
  //w->writeconst(index3_t{0,0,0}, w->size(), &sixteen);
  w->writeconst(index3_t{0,0,0}, index3_t{64,64,64}, &sixteen);
  w->writeconst(index3_t{0,0,64}, index3_t{64,64,512-64}, &sixteen);
  w->close();

  // The file now has inflated both lod0 bricks and lowres bricks.
  // This will look like real data was written including lowres data.
  // sort of is. Looking at the statistics would show that the file
  // was actually cleared. LodMode::Never would have worked, but gets
  // refused. While LodMode::Early will be accepted. So, the opposite
  // of the previous case.

  r = IZgyReader::open(lad.name());
  TEST_EQUAL(r->statistics().min, 0); // 16 if clearing in one step.
  TEST_EQUAL(r->statistics().max, 42); // 16 if clearing in one step.
  r->close();

  must_throw("already has lowres", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .filename(lad.name()));
      w->close();
    });

  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Early)
     .filename(lad.name()));
  w->close();
}

/**
 * A file that is finalized and contains compressed data should not be
 * opened for update at all. LodMode::Never is disallowed because it
 * would leave stale lowres behind, and the others are disallowed
 * because the lowres cannot be updated. It is not important exactly
 * where this is caught and what the error message is, but it should
 * happen before any writes to the file. So the file doesn't get
 * corrupted. Caught during lowres write, e.g. exception "*in update
 * mode*" is definitely to late. For cloud I/O it might get even
 * trickier because the last segment is reopened very early.
 */
static void
test_lodmode7()
{
  typedef std::array<std::int64_t,3> index3_t;
  Test_Utils::LocalFileAutoDelete lad("lodmode7.zgy");

  if (verbose())
    std::cout << std::endl;

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(192, 640, 100)
     .zfp_compressor(70)
     .zfp_lodcompressor(30)
     .lodmode(LodMode::Early)
     .filename(lad.name()));

  const float sixteen{16};
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, w->size(), &sixteen); // 60 bricks
  w->writeconst(index3_t{0,0,0}, index3_t{60,100,70}, &fortytwo); // 4 bricks
  const std::int64_t size4c = w->filestats()->brickCompressedSize();
  w->finalize();
  const std::int64_t size8c = w->filestats()->brickCompressedSize();
  w->close();

  // Numbers from a test run, seems legit. ~21 KB/brick for lod0,
  // ~6 KB/brick for the higher compressed lowres bricks.
  TEST_EQUAL(size4c, 84560);
  TEST_EQUAL(size8c, 25264+84560);

  if (verbose()) {
    std::cout << "Size of 4 lod0 " << size4c
              << ", 4 lowres " << size8c - size4c
              << std::endl;
  }

  // lod1 size in bricks: 2 x 5 x 1 = 10
  // lod2 size in bricks: 1 x 3 x 1 = 3
  // lod3 = 2, lod4 = 1, total 16 + the 60 fullres = 76
  // Expect just 1 dirty lod per level, 4 total dirty lods and 12 clean.

  std::shared_ptr<IZgyReader> r = IZgyReader::open(lad.name());
  std::shared_ptr<const FileStatistics> stat = r->filestats();
  r->close();

  TEST_EQUAL(stat->fileVersion(),               4);
  TEST_EQUAL(stat->fileSize(),            1158400); // Little more than 1 brick.
  TEST_EQUAL(stat->dataStart(),       1*1024*1024); // Headers always padded
  TEST_EQUAL(stat->brickNormalCount(),          0);
  TEST_EQUAL(stat->brickCompressedCount(),      8);
  TEST_EQUAL(stat->brickCompressedSize(),  109824); // ~13 KB / brick
  TEST_EQUAL(stat->brickMissingCount(),         0);
  TEST_EQUAL(stat->brickConstantCount(),       68);
  TEST_CHECK(stat->isCompressed());

  // May get either of these errors; which one doesn't matter:
  // A finalized file cannot have compressed data appended.
  // LodMode::Never not allowed on a file that already has lowres.
  must_throw("finalized file cannot have compressed data", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .zfp_compressor(70)
         .filename(lad.name()));
    });

  // May get either of these errors; which one doesn't matter:
  // A finalized compressed file cannot be opened for update."
  // LodMode::Never not allowed on a file that already has lowres.
  must_throw("finalized compressed file cannot be", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .filename(lad.name()));
    });

    w = nullptr;
    must_throw("finalized compressed file cannot be", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Early1)
         .filename(lad.name()));
    });

  // Should not have managed to open the file above. But if we did,
  // (remove the last check in ZgyInternalMeta::initFromReopen)
  // then make sure the actual write is vetoed because we are
  // trying to overwrite an existing compressed brick.
  // A similar test is done in "api.lodmode9" and that one doesn't
  // require patching the code.
  if (w) {
    const float twenty{20};
    w->writeconst(index3_t{40,0,0}, index3_t{1,1,1}, &twenty);
    // Should not have managed to write, above. If it worked, then:
    // Currently opened with no compression, so the rewritten first
    // brick and the rewritten 4 lowres bricks will all be full size.
    // The 5 compressed bricks will get leaked. Which is NOT ALLOWED.
    // New file size: pad to one brick (1MB), add 5 bricks, total 7.
    w->finalize();
    w->close();
  }

  // The next few checks are only meaningful if the two checks that
  // ought to have thrown, didn't.
  if (w) {
    r = IZgyReader::open(lad.name());
    stat = r->filestats();
    r->close();

    TEST_EQUAL(stat->fileVersion(),               4);
    TEST_EQUAL(stat->fileSize(),        7*1024*1024);
    TEST_EQUAL(stat->dataStart(),       1*1024*1024 + 23056); // Leaked some...
    TEST_EQUAL(stat->brickNormalCount(),          5);
    TEST_EQUAL(stat->brickCompressedCount(),      3);
    // Issue with the heuristic for brick sizes. The padding added to
    // make the first uncompressed brick become aligned is assumed to
    // belong to the last uncompressed one. filestats() isn't expected
    // to be 100% correct, so this is not a bug.
    TEST_EQUAL(stat->brickCompressedSize(), 1025520);
    TEST_EQUAL(stat->brickMissingCount(),         0);
    TEST_EQUAL(stat->brickConstantCount(),       68);
    TEST_CHECK(stat->isCompressed());
  }
}

/**
 * An uncompressed file that is finalized should not be reopened
 * with the lodcompressor turned on. Because the lowres will need
 * to be rebuilt, and it will not be allowed to overwrite any
 * regular brick with a compressed one. Technically, leaving the
 * lowres uncompressed and the fullres compressed is valid.
 * But silly.
 *
 * TODO-WIP-BrickedAPI: Low: Ideally an exception should be thrown as
 * soon as possible, before any write is done. This is particularly
 * difficult for cloud I/O. And on the flip side, it *is* allowed to
 * reopen if the only purpose is to update certain metadata. But why
 * bother turning on compression in that case?
 */
static void
test_lodmode8()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("lodmode8.zgy");

  if (verbose())
    std::cout << std::endl;

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(192, 640, 100)
     .lodmode(LodMode::Early)
     .filename(lad.name()));

  const float sixteen{16};
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, w->size(), &sixteen); // 60 bricks
  w->writeconst(index3_t{0,0,0}, index3_t{60,100,70}, &fortytwo); // 4 bricks
  w->finalize();
  w->close();
  w = nullptr;

  // Arguably this throws too early, because if the sole purpose
  // was to update metadata then LodMode::Never is the intuitive
  // choice. Maybe I need a LodMode::IWillNotWriteAnything?
  must_throw("finalized file cannot have compressed data appended", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Never)
         .zfp_compressor(70)
         .filename(lad.name()));
    });

  // After enabling the last check in ZgyInternalMeta::initFromReopen,
  // it is not possible to update meta in a compressed finalized file.
  // This is a (missing) feature, not a bug.
  must_throw("finalized file cannot have compressed data appended", [&]()
    {
      w = IZgyWriter::reopen
        (ZgyWriterArgsV2()
         .lodmode(LodMode::Early1)
         .zfp_compressor(70)
         .filename(lad.name()));
    });

  // And this might be too late?
  if (w) {
    must_throw("in update mode", [&]()
      {
        const float twenty{20};
        w->writeconst(index3_t{40,0,0}, index3_t{1,1,1}, &twenty);
      });
    w->finalize();
    w->close();
  }
}

/**
 * A file that is NOT finalized and contains compressed data can be
 * updated, should not allow updating any of the compressed bricks
 * regardless of whether the new data is compressed or not. Because
 * the old bulk data would be leaked.
 *
 * Current BUG:
 */
static void
test_lodmode9()
{
  typedef std::array<std::int64_t,3> index3_t;
  Test_Utils::LocalFileAutoDelete lad("lodmode9.zgy");

  if (verbose())
    std::cout << std::endl;

  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .size(192, 640, 100)
     .zfp_compressor(70)
     .zfp_lodcompressor(30)
     .lodmode(LodMode::Never)
     .filename(lad.name()));

  const float sixteen{16};
  const float fortytwo{42};
  w->writeconst(index3_t{0,0,0}, w->size(), &sixteen); // 60 bricks
  w->writeconst(index3_t{0,0,0}, index3_t{60,100,70}, &fortytwo); // 4 bricks
  w->finalize();
  w->close();

  // Reopen with compression.
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Never)
     .zfp_compressor(50)
     .filename(lad.name()));

  must_throw("compressed brick with compressed data in update mode constant", [&]()
    {
      const float twenty{20};
      w->writeconst(index3_t{40,0,0}, index3_t{1,1,1}, &twenty);
    });
  w->finalize();
  w->close();

  // Reopen without compression.
  w = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .lodmode(LodMode::Never)
     .filename(lad.name()));

  // FAILS! Ought to have thrown. See bug description in the
  // OpenForUpdate overload of the ZgyWriter constructor.
  // If this does not throw, space on the file will leak.
  //must_throw("something else", [&]()
  //{
      const float seventeen{20};
      w->writeconst(index3_t{40,0,0}, index3_t{1,1,1}, &seventeen);
  //});
  w->finalize();
  w->close();
}

/**
 * Test for Issue 2044116.
 *
 * Allowing a cube to be written without low resolution data is a new
 * feature in OpenZGY Some fuzzy logic is needed to decide whether low
 * resolution data exists or not. Issue 2044116 found that the fuzzy
 * logic wasn't quite correct.
 *
 * Run a test writing and reading back the single, provided trace
 * using the provided lodmode. Check the results.
 *
 * Variations to run:
 *
 * (1) All constant cube, lowres built:
 *     Can read lowres AND can re-finalize.
 * (2) All constant cube, lowres not built:
 *     can read lowres AND can finalize.
 * (3) Survey cleared, data written, no lowres generated:
 *     Cannot read lowres, can finalize.
 * (4) Real data written and finalized but top level ended up all const:
 *     Can read lowres, cannot finalize.
 * (5) Real data written and finalized but all lods became all const:
 *     Can read lowres, AND can finalize.
 * (6) Real data written and not finalized. Nothing special:
 *     Cannot read lowres, can finalize.
 * (7) Real data written and finalized. Nothing special:
 *     Can read lowres, cannot finalize.
 *
 * Note that these tests partly overlap api.lodmode{1..9} above.
 */
static bool
do_nlods_cornercase(
     const std::string& filename,
     const LodMode lodmode,
     const bool clear_first,
     const std::vector<std::int8_t>& data,
     const std::string& expect_lowres_const,
     const bool expect_need_lowres)
{
  typedef std::array<std::int64_t,3> index3_t;
  std::shared_ptr<IZgyWriter> w = IZgyWriter::open
    (ZgyWriterArgsV2()
     .filename(filename)
     .lodmode(lodmode)
     .bricksize(16, 16, 32)
     .size(16, 16, 7*32)
     .datatype(SampleDataType::int8)
     .datarange(-128, +127)
     .decimation(std::vector<DecimationType>{DecimationType::Decimate}));
  if (clear_first) {
    std::int8_t dflt{42};
    w->writeconst(index3_t{0,0,0}, w->size(), &dflt);
  }
  if (data.size() != 0)
    w->write(index3_t{0,0,0}, index3_t{1, 1, (std::int64_t)data.size()}, data.data());
  w->finalize();
  w->close();
  w.reset();

  // Read back fullres and each lowres levels. Make a list, as a string,
  // where a number indicates that this lod level in the survey is all
  // constant, with samples set to "number". Levels that are not all
  // constant are represented by an asterisk.
  std::shared_ptr<IZgyReader> r = IZgyReader::open(filename, nullptr);
  std::stringstream actual_lowres_const;
  actual_lowres_const << "[";
  for (std::int32_t lod = 0; lod < (std::int32_t)r->nlods(); ++lod) {
    const index3_t size_at_lod
      {r->brickcount()[lod][0] * r->bricksize()[0],
       r->brickcount()[lod][1] * r->bricksize()[1],
       r->brickcount()[lod][2] * r->bricksize()[2]};
    std::pair<bool,double> isconst =
      r->readconst(index3_t{0,0,0}, size_at_lod, lod, false);
    actual_lowres_const << (lod==0 ? "" : ", ")
                        << (isconst.first ?
                            std::to_string((int)isconst.second) :
                            std::string("*"));
  }
  actual_lowres_const << "]";
  r->close();
  r.reset();

  bool actual_need_lowres = false;
  try {
    std::shared_ptr<IZgyWriter> u = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .filename(filename)
       .lodmode(LodMode::Never));
    u->close();
    u.reset();
  }
  catch (const OpenZGY::Errors::ZgyUserError& ex) {
    const std::string message{ex.what() ? ex.what() : "(null)"};
    if (message.find("LodMode::Never not allowed") == std::string::npos)
      throw;
    actual_need_lowres = true;
  }

  // This checks both the reported number of lods that can be read
  // (just count the number of elements) and the constant-value,
  // if any, for each level. The latter check is mainly to ensure
  // that I am testing what I think I am.
  const bool ok1 = TEST_EQUAL(actual_lowres_const.str(), expect_lowres_const);

  // It is illegal to change a file from having lowres to not having it.
  // Normally this means that if expect_lowres_const has more than one
  // element, trying to open with lowres disabled should throw. But there
  // are special cases (tiny survey, all-const survey) where reading
  // lowres is allowed and clearing the lowres is also valid.
  const bool ok2 = TEST_EQUAL(actual_need_lowres, expect_need_lowres);

  if (ok1 && ok2 && verbose())
    std::cout << "\nConstant flag by lod level: " << actual_lowres_const.str();

  return ok1 && ok2;
}

/**
 * Just clear the survey, don't write any real data.
 * Request lowres generation; that should not matter either way.
 *
 * (1) All constant cube, lowres built:
 *     Can read lowres AND can re-finalize.
 */
static void
test_nlods_cornercase_1()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_1.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Early, /*clear_first=*/true,
     std::vector<std::int8_t>{},
     std::string{"[42, 42, 42, 42]"},
     /*expect_need_lowres=*/false);
  TEST_CHECK(ok);
}

/**
 * Just clear the survey, don't write any real data.
 * Don't request lowres generation; that should not matter either way.
 *
 * (2) All constant cube, lowres not built:
 *     can read lowres AND can finalize.
 */
static void
test_nlods_cornercase_2()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_2.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Never, /*clear_first=*/true,
     std::vector<std::int8_t>{},
     std::string{"[42, 42, 42, 42]"},
     /*expect_need_lowres=*/false);
  TEST_CHECK(ok);
}

/**
 * Try to trick the code to believe there is valid lowres because the
 * top level is Constant, not Missing. Since there is no valid lowres,
 * opening with LodMode::Never should be allowed. Remember that the
 * opposite is not true.
 *
 * (3) Survey cleared, data written, no lowres generated:
 *     Cannot read lowres, can finalize.
 *
 * Logging at level 1 should say "lodN switched to Missing",
 * meaning that we hit the special case in _close_internal()
 * where the "phantom" lowres was explicitly deleted.
 */
static void
test_nlods_cornercase_3()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_3.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Never, /*clear_first=*/true,
     std::vector<std::int8_t>{0,0,1,0},
     std::string{"[*]"},
     /*expect_need_lowres=*/false);
  TEST_CHECK(ok);
}

/**
 * Test specifically for Issue 2044116.
 *
 * This survey is so sparse that decimation at higher levels end up
 * as constant zero. On read, the application still needs to be told
 * that lowres is available.
 *
 * (4) Real data written and finalized but top level ended up all const:
 *     Can read lowres, cannot finalize.
 */
static void
test_nlods_cornercase_4()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_4.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Early, /*clear_first=*/false,
     std::vector<std::int8_t>{0,0,1,0},
     std::string{"[*, *, 0, 0]"},
     /*expect_need_lowres=*/true);
  TEST_CHECK(ok);
}

/**
 * Check that the fix for Issue 2044116 didn't break something else.
 *
 * This survey is so sparse that *all* decimation ends up constant zero.
 * On read, the application still needs to be told that lowres is available.
 * But since there is no decimation yet, the application is still allowed
 * to open with lowres turned off.
 *
 * (5) Real data written and finalized but all lods became all const:
 *     Can read lowres, AND can finalize.
 *     Or not... It doesn't really matter if the app in this case
 *     isn't allowed to use Never because the file was in fact finalized
 *     by the application. It is an implementation detail that we mightis sort of finalized.
 *
 * Logging at level 1 says "lodN switched to Missing",
 * meaning that we hit the special case in _close_internal()
 * where a "phantom" lowres was explicitly deleted. The reson
 * for this is subtle. The common test routine tries to open
 * the file for update, with no lod generation, just to see if
 * this was allowed. This will actually convert the file from
 * one where lowres can be accessed into one that it cannot.
 * Nobody should complain because the app did in fact ask for
 * LodMode::Never. But obvious, it is not.
 */
static void
test_nlods_cornercase_5()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_5.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Early, /*clear_first=*/false,
     std::vector<std::int8_t>{0,1,0,0},
     std::string{"[*, 0, 0, 0]"},
     /*expect_need_lowres=*/true);
  TEST_CHECK(ok);
}

/**
 * Test this common case. The top brick will be missing because neither
 * survey clear nor lowres generation was done.
 *
 * (6) Real data written and not finalized. Nothing special:
 *     Cannot read lowres, can finalize.
 */
static void
test_nlods_cornercase_6()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_6.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Never, /*clear_first=*/false,
     std::vector<std::int8_t>{1,2,3,4,5,6,7,8,9},
     std::string{"[*]"},
     /*expect_need_lowres=*/false);
  TEST_CHECK(ok);
}

/**
 * Test this common case. The top brick will be missing because neither
 * survey clear nor lowres generation was done.
 *
 * (7) Real data written and finalized. Nothing special:
 *     Can read lowres, cannot finalize.
 */
static void
test_nlods_cornercase_7()
{
  Test_Utils::LocalFileAutoDelete lad("nlods_cornercase_7.zgy");
  const bool ok = do_nlods_cornercase
    (lad.name(), LodMode::Early, /*clear_first=*/false,
     std::vector<std::int8_t>{1,2,3,4,5,6,7,8,9},
     std::string{"[*, *, *, *]"},
     /*expect_need_lowres=*/true);
  TEST_CHECK(ok);
}

/**
 * OpenZGY normally manages guids itself, but the application can
 * explicitly ask to copy guids (actually just "verid") from an open
 * file. The test is on-prem only, with no bulk data needed.
 *
 * OpenZGY will not change guids (or most other metadata) after
 * a file has been open. So, it won't try things like not making a
 * new verid if nothing actually changed. Test this as well.
 */
static void
test_setguids()
{
  LocalFileAutoDelete lad1("setguids_1.zgy");
  LocalFileAutoDelete lad2("setguids_2.zgy");
  LocalFileAutoDelete lad3("setguids_3.zgy");

  // Create a file that we will be copying the guid from.
  // First, verify the existing rules for guids.
  // They change on any open for create or update, and don't
  // change on open for read.
  std::shared_ptr<IZgyReader> reader;
  {
    auto w1 = IZgyWriter::open(ZgyWriterArgsV2()
                               .filename(lad1.name())
                               .size(7, 9, 113));
    std::string verid = w1->verid();
    w1->close();
    // Safe to read the verid from the writer; it would already be correct.
    reader = IZgyReader::open(lad1.name());
    TEST_EQUAL(reader->verid(), verid);
    reader->close();
    reader.reset();
    w1 = IZgyWriter::reopen(ZgyWriterArgsV2().filename(lad1.name()));
    TEST_CHECK(w1->verid() != verid);
    verid = w1->verid();
    w1->close();
    // Another check that the verid is decided immediately on opening.
    reader = IZgyReader::open(lad1.name());
    TEST_EQUAL(reader->verid(), verid);
  }

  // Create file #2, copying guids from file #1. Verify before & after close.
  {
    auto w2 = IZgyWriter::open(ZgyWriterArgsV2()
                               .filename(lad2.name())
                               .size(7, 9, 113)
                               .guidsfrom(reader));
    TEST_EQUAL(w2->verid(), reader->verid());
    w2->close();
    auto r2 = IZgyReader::open(lad2.name());
    TEST_EQUAL(r2->verid(), reader->verid());
    r2->close();
  }

  // Test updating the guids later.
  // Create and close file #3 with no special guid handling. Then,
  // Open file #3 for update & copy file #1's guids. Verify before
  // and after file #3 close.
  {
    auto w3 = IZgyWriter::open(ZgyWriterArgsV2()
                               .filename(lad3.name())
                               .size(7, 9, 113));
    w3->close();
    w3 = IZgyWriter::reopen(ZgyWriterArgsV2()
                            .filename(lad3.name())
                            .guidsfrom(reader));
    TEST_EQUAL(w3->verid(), reader->verid());
    w3->close();
    auto r3 = IZgyReader::open(lad3.name());
    TEST_EQUAL(r3->verid(), reader->verid());
    r3->close();
  }
  reader->close();
}

static std::uint64_t
filetime(const std::string& filename)
{
#ifdef _WIN32
  struct _stat64 st;
  if (!TEST_CHECK(::_stat64(filename.c_str(), &st) >= 0))
    return 0;
  return static_cast<std::uint64_t>(st.st_mtime);
#else
  struct stat st;
  if (!TEST_CHECK(::stat(filename.c_str(), &st) >= 0))
    return 0;
  return static_cast<std::uint64_t>(st.st_mtime);
#endif
}


/**
 * If a file is opened for write but never written to and never
 * explicitly finalized, will it be touched at all? This can be
 * important if an application opens two write handles on the same
 * file. Even if just one of them is used for writing, the other might
 * decide it wants to flush the unchanged metadata and lookup tables
 * anyway.
 *
 * Answer: Yes. A new verid is set immediately on reopen, and that
 * is visible to the application. So it is unsafe to then decide
 * we didn't want the verid to be updated anyway.
 *
 * For sdms files, there is an additional reason. The last segment
 * will be copied to memory and deleted in the file. Ideally that
 * ought to be deferred util the first write. It isn't.
 */
static void
test_notouch()
{
  LocalFileAutoDelete lad("notouch.zgy");
  auto w = IZgyWriter::open(ZgyWriterArgsV2()
                            .filename(lad.name())
                            .size(7, 9, 113));
  std::string oldid = w->verid();
  w->close();
  std::uint64_t oldtime = filetime(lad.name());
  std::this_thread::sleep_for(std::chrono::seconds(2));
  w = IZgyWriter::reopen(ZgyWriterArgsV2().filename(lad.name()));
  std::string newid = w->verid();
  w->close();
  std::uint64_t newtime = filetime(lad.name());
  //TEST_EQUAL(newtime, oldtime);
  //TEST_EQUAL(newid, oldid);
  // It does actually change...
  TEST_CHECK(newtime != oldtime);
  TEST_CHECK(newid != oldid);
}

static void
do_patched_file(
     const std::string& expected_error,
     const std::function<void(char *data)>& patch)
{
  LocalFileAutoDelete lad("maliciously_patced.zgy");
  std::string dstname = lad.name();
  std::string srcname = get_testdata("Fancy-int8.zgy");
  const std::int64_t bufsize(5*256*1024); // 4 bricks + 1 for headers
  std::unique_ptr<char[]> data(new char[bufsize]);
  errno = 0;
  std::ifstream src(srcname, std::ios_base::in|std::ios_base::binary);
  const int src_errno = errno;
  if (!TEST_CHECK(src.good())) {
    std::cout << "ERROR: copy_file cannot read \"" << srcname << "\""
              << " errno " << src_errno << std::endl;
    return;
  }
  errno = 0;
  std::ofstream dst(dstname, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
  const int dst_errno = errno;
  if (!TEST_CHECK(dst.good())) {
    std::cout << "ERROR: copy_file cannot write \"" << dstname << "\""
              << " errno " << dst_errno << std::endl;
    return;
  }
  src.read(data.get(), bufsize);
  patch(data.get());
  dst.write(data.get(), src.gcount());
  TEST_CHECK(dst.good());
  TEST_CHECK(src.good());
  dst.close();
  src.close();

  try {
    auto reader = OpenZGY::IZgyReader::open(dstname);
    TEST_EQUAL(std::string("<no error>"), expected_error);
  }
  catch (const OpenZGY::Errors::ZgyError& ex) {
    if (std::string(ex.what()).find(expected_error) == std::string::npos)
      TEST_EQUAL(std::string(ex.what()), expected_error);
  }
}

static void
test_malicious_neg_bs()
{
  // bricksize[0] @ byte 9..12 set msb to a negative number.
  do_patched_file("bricksize must be >= 4", [](char* data)
                  {
                    data[12] = -100;
                  });
}

static void
test_malicious_zer_bs()
{
  // bricksize[0] @ byte 9..12 set to zero. Rest left alone.
  do_patched_file("bricksize must be >= 4", [](char* data)
                  {
                    data[9] = 0;
                    data[10] = 0;
                    data[11] = 0;
                    data[12] = 0;
                  });
}

static void
test_malicious_zero_bs()
{
  // bricksize[*] set to zero.
  do_patched_file("bricksize must be >= 4", [](char* data)
                  {
                    data[9] = 0;
                    data[10] = 0;
                    data[11] = 0;
                    data[12] = 0;
                    data[13] = 0;
                    data[14] = 0;
                    data[15] = 0;
                    data[16] = 0;
                    data[17] = 0;
                    data[18] = 0;
                    data[19] = 0;
                    data[20] = 0;
                  });
}

static void
test_malicious_huge_bs()
{
  // bricksize[] too big for the product to fit in an int32
  do_patched_file("Total brick size in bytes must be <= 0x7fffffff", [](char* data)
                  {
                    data[9] = 0;
                    data[10] = 0;
                    data[11] = 16;
                    data[12] = 0;
                    data[13] = 0;
                    data[14] = 0;
                    data[15] = 8;
                    data[16] = 0;
                    data[17] = 0;
                    data[18] = 0;
                    data[19] = 32;
                    data[20] = 0;
                  });
}

static void
test_malicious_neg_size()
{
  // size[0] @ byte 103..106 set msb to a negative number.
  do_patched_file("survey size[0] is -1677721488 but must be between 1 and 0x7fffffff", [](char* data)
                  {
                    data[106] = -100;
                  });
}

static void
test_malicious_zer_size()
{
  // bricksize[0] @ byte 9..12 set to zero. Rest left alone.
  do_patched_file("survey size[0] is 0 but must be between 1 and 0x7fffffff", [](char* data)
                  {
                    std::fill(&data[103], &data[103+4], (char)0);
                  });
}

static void
test_malicious_zero_size()
{
  // bricksize[*] set to zero.
  do_patched_file("survey size[0] is 0 but must be between 1 and 0x7fffffff", [](char* data)
                  {
                    std::fill(&data[103], &data[103+12], (char)0);
                  });
}

static void
test_malicious_huge_size()
{
  // size[] too big for the product to fit in an int64
  do_patched_file("size in bytes must be <= 0x7fffffffffffffff", [](char* data)
                  {
                    data[103] = 1;
                    data[104] = 0;
                    data[105] = 0;
                    data[106] = 64;
                    data[107] = 2;
                    data[108] = 0;
                    data[109] = 0;
                    data[110] = 64;
                    data[111] = 42;
                    data[112] = 0;
                    data[113] = 0;
                    data[114] = 64;
                  });
}

static void
test_malicious_oom_size()
{
  // size technically valid, but might give an oom?
  // If there has been more explicit tests for size of lut
  // vs. size of file then this would be more dependable.
  do_patched_file("is past EOF", [](char* data)
                  {
                    data[103] = 0;
                    data[104] = 0;
                    data[105] = 16;
                    data[106] = 0;
                    data[107] = 0;
                    data[108] = 0;
                    data[109] = 16;
                    data[110] = 0;
                    data[111] = 0;
                    data[112] = 0;
                    data[113] = 32;
                    data[114] = 0;
                  });
}

/**
 * An entry in the brick lookup table is past eof.
 */
static void
test_malicious_bad_brick()
{
  // Note: unit test .verbose meta.open_fancy os useful to extract
  // which offset to use.
  //
  // An entry in the brick lookup table is past eof.
  // In this case, patch the single brick as the highest lod.
  //
  // The library cannot know whether the problem is a truncated file
  // or a corrupt pointer. It will report the former. Because that is
  // more likely.
  do_patched_file("file is truncated", [](char* data)
                  {
                    // Set file pointer to zero first
                    std::fill(&data[0x09c1], &data[0x09c1+8], (char)0);
                    data[0x09c1+6] = 42;
                  });
}

/**
 * As test_malicious_bad_brick, but for the alpha lookup.
 * Do we even care? Is the alpha table still read in?
 */
static void
test_malicious_bad_alpha()
{
  do_patched_file("file is truncated", [](char* data)
                  {
                    // Set file pointer to zero first
                    std::fill(&data[0x09a1], &data[0x09a1+8], (char)0);
                    data[0x09a1+6] = (char)42;
                  });
}

/**
 * Two brick- or alpha pointers refer to the same physical location.
 * Or not identical but still overlapping.
 */
static void
test_malicious_overlap()
{
  // Tip: run "test_all .verbose meta.open_fancy" to get offsets.

  do_patched_file("Bricks overlap", [](char* data)
                  {
                    char *bricklut = data + 0x09c1;
                    std::copy(bricklut, bricklut+8, bricklut+16);
                  });

  // The code won't check overlaps in the alpha table or between tables.
  // Alpha isn't actually used. So, we don't care.
  do_patched_file("<no error>", [](char* data)
                  {
                    char *alphalut = data + 0x09a1;
                    std::copy(alphalut, alphalut+8, alphalut+16);
                  });

  do_patched_file("<no error>", [](char* data)
                  {
                    char *alphalut = data + 0x09a1;
                    char *bricklut = data + 0x09a1;
                    std::copy(bricklut, bricklut+8, alphalut);
                  });
}

/**
 * String buffer problems
 */
static void
test_malicious_slbufsize()
{
  do_patched_file("is past EOF", [](char* data)
                  {
                    // offset of info header plus offset inside info header.
                    char *slbufsize_ptr = data + 0x156; //0x9 + 0x14d;
                    slbufsize_ptr[2] = 42;
                  });

  // slbufsize is normally treated as unsigned, so it might be caught
  // by the "is past EOF" check instead of making the string table overlap
  // something else. And since it is an uint32 it shouldn't overflow
  // the offset calculation either.
  do_patched_file("is past EOF", [](char* data)
                  {
                    std::uint8_t *slbufsize_ptr = (std::uint8_t*)data + 0x156;
                    // slbufsize = -16
                    slbufsize_ptr[0] = 0xf0;
                    slbufsize_ptr[1] = 0xff;
                    slbufsize_ptr[2] = 0xff;
                    slbufsize_ptr[3] = 0xff;
                  });
}

#if 0 // Dec 2024, why was this never enabled?
/**
 * Strings not null terminated.
 */
static void
test_malicious_asciiz()
{
  do_patched_file("?", [](char* data)
                  {
                    // offset of info header plus offset inside info header.
                    char *slbufsize_ptr = data + 0x156; //0x9 + 0x14d;
                    slbufsize_ptr[2] = 42;
                  });
}
#endif

class Register
{
public:
  Register()
  {
    register_test("api.dump",                test_dump);
    register_test("api.zgywriterargs",       test_ZgyWriterArgs);
    // Locking is disabled, set to complain-only mode.
    //register_test("api.locks",               test_locks);
    register_test("api.readmeta_r",          test_readmeta_r);
    register_test("api.readmeta_w",          test_readmeta_w);
    register_test("api.readmeta_v1_r",       test_readmeta_v1_r);
    // Test file is not checked in yet.
    //register_test("api.readcmeta",           test_readcmeta);
    register_test("api.readconst",           test_readconst);
    register_test("api.readbulk",            test_readbulk);
    register_test("api.readwrite_b",         test_readwrite_int8);
    register_test("api.readwrite_s",         test_readwrite_int16);
    register_test("api.readsubtiles",        test_readsubtiles);
    register_test("api.readbadvt",           test_readbadvt);
    register_test("api.readbadpos",          test_readbadpos);
    register_test("api.readnotopen",         test_readnotopen);
    register_test("api.writenotopen_b",      test_writenotopen<std::int8_t>);
    register_test("api.writenotopen_s",      test_writenotopen<std::int16_t>);
    register_test("api.writenotopen_f",      test_writenotopen<float>);
    register_test("api.createargs",          test_createargs);
    register_test("api.ioerror",             test_ioerror);
    register_test("api.finalize_1",          test_finalize<1>);
    register_test("api.finalize_2",          test_finalize<2>);
    register_test("api.finalize_3",          test_finalize<3>);
    register_test("api.finalize_4",          test_finalize<4>);
    register_test("api.finalize_5",          test_finalize<5>);
    register_test("api.finalize_6",          test_finalize<6>);
    register_test("api.finalize_7",          test_finalize<7>);
    register_test("api.genlod",              test_genlod);
    register_test("api.genlod2",             test_genlod2);
    register_test("api.write",               test_write);
    register_test("api.compress_noop",       test_compress_noop);
    register_test("api.compress_zfp",        test_compress_zfp);
    register_test("api.compress_off",        test_compress_off);
    register_test("api.inflated_constant",   test_inflated_constant);
    register_test("api.stats_after_clip",    test_stats_after_clip);
#ifdef HAVE_SD
    register_sd_test("api.write_cloud",         test_write_cloud);
    register_sd_test("api.write_cloud_mt",      test_write_cloud_mt);
    register_sd_test("api.alturl",              test_alturl);
    register_sd_test("api.idtoken",             test_idtoken);
    register_sd_test("api.basicinfo",           test_basicinfo);
    register_sd_test("api.sharecred",           test_sharecred);
#endif
    register_test("api.noinfo",              test_noinfo);
    register_test("api.historange",          test_historange);
    register_test("api.lod_lowpass",         test_lod_lowpass);
    register_test("api.lod_weighted",        test_lod_weighted);
    register_test("api.lod_average",         test_lod_average);
    register_test("api.copy",                test_copy);
    register_test("api.enums",               test_enums);
    register_test("api.dummy_compress",      test_dummy_compress);
    // SINGLEPASS: No longer useful as a regression test.
    //register_test("api.histo_cornercase_f",  test_histo_cornercase_float);
    register_test("api.histo_cornercase_i",  test_histo_cornercase_int);
    register_test("api.filestats",           test_filestats);
    register_test("api.transform_r",         test_transform_r);
    register_test("api.transform_w",         test_transform_w);
    register_test("api.all_exceptions",      test_all_exceptions);
    register_test("api.ambig1",              test_ambig1);
    register_test("api.ambig2",              test_ambig2);
    register_test("api.ambig3",              test_ambig3);
    register_test("api.2d",                  test_2d);
    register_test("api.decimate_edge",       test_decimate_edge);
#ifdef HAVE_SD
    register_sd_test("api.tokencb2",            test_tokencb2);
#endif
    register_test("api.readwrite",           test_readwrite_local);
#ifdef HAVE_SD
    register_sd_test("api.readwrite_cloud",     test_readwrite_cloud);
    register_sd_test("api.hammer",              test_hammer);
    register_sd_test("api.sderrors",            test_sderrors);
#endif
    register_test("api.edgebricks",          test_edgebricks);
    register_test("api.bat_local_1",         test_bat_local_1);
    register_test("api.bat_local_2",         test_bat_local_2);
    register_test("api.bat_local_4",         test_bat_local_4);
    register_test("api.bat_local_zfp",       test_bat_local_zfp);
    register_test("api.stats_empty",         test_stats_empty);
    register_test("api.stats_one",           test_stats_one);
    register_test("api.stats_many",          test_stats_many);
#ifdef HAVE_SD
    register_sd_test("api.bat_sd_1",            test_bat_sd_1);
    register_sd_test("api.bat_sd_2",            test_bat_sd_2);
    register_sd_test("api.bat_sd_4",            test_bat_sd_4);
    register_sd_test("api.bat_sd_zfp",          test_bat_sd_zfp);
    register_sd_test("api.roflag",              test_roflag);
    register_sd_test("api.roflag_nolock",       test_roflag_nolock);
    register_sd_test("api.roflag_unlock",       test_roflag_unlock);
    register_sd_test("api.client_cred",         test_client_cred);
#if 0 // Not usable as an automated test
    register_sd_test("bug.671969",              test_bug_671969);
#endif
    register_sd_test("api.list_segments",       test_list_segments);
    register_sd_test("api.cp",                  test_cp);
#endif
    register_test("api.build_hist_dynamic",  test_build_hist_dynamic);
    register_test("api.build_hist_static",   test_build_hist_static);
    register_test("api.build_hist_with_vt",  test_build_hist_with_vt);
    register_test("api.nlods",               test_nlods);
    register_test("api.incrlod1",            test_incrlod1);
    register_test("api.incrlod2",            test_incrlod2);
    register_test("api.incrlod3",            test_incrlod3);
    register_test("api.lodmode1",            test_lodmode1);
    register_test("api.lodmode2",            test_lodmode2);
    register_test("api.lodmode3",            test_lodmode3);
    register_test("api.lodmode4",            test_lodmode4);
    register_test("api.lodmode5",            test_lodmode5);
    register_test("api.lodmode6",            test_lodmode6);
    register_test("api.lodmode7",            test_lodmode7);
    register_test("api.lodmode8",            test_lodmode8);
    register_test("api.lodmode9",            test_lodmode9);
    register_test("api.nlods_cornercase_1",  test_nlods_cornercase_1);
    register_test("api.nlods_cornercase_2",  test_nlods_cornercase_2);
    register_test("api.nlods_cornercase_3",  test_nlods_cornercase_3);
    register_test("api.nlods_cornercase_4",  test_nlods_cornercase_4);
    register_test("api.nlods_cornercase_5",  test_nlods_cornercase_5);
    register_test("api.nlods_cornercase_6",  test_nlods_cornercase_6);
    register_test("api.nlods_cornercase_7",  test_nlods_cornercase_7);
    register_test("api.setguids",            test_setguids);
    register_test("api.notouch",             test_notouch);
    register_test("api.malicious.neg_bs",    test_malicious_neg_bs);
    register_test("api.malicious.zer_bs",    test_malicious_zer_bs);
    register_test("api.malicious.zero_bs",   test_malicious_zero_bs);
    register_test("api.malicious.huge_bs",   test_malicious_huge_bs);
    register_test("api.malicious.neg_size",  test_malicious_neg_size);
    register_test("api.malicious.zer_size",  test_malicious_zer_size);
    register_test("api.malicious.zero_size", test_malicious_zero_size);
    register_test("api.malicious.huge_size", test_malicious_huge_size);
    register_test("api.malicious.oom_size",  test_malicious_oom_size);
    register_test("api.malicious.bad_brick", test_malicious_bad_brick);
    register_test("api.malicious.bad_alpha", test_malicious_bad_alpha);
    register_test("api.malicious.overlap",   test_malicious_overlap);
    register_test("api.malicious.slbufsize", test_malicious_slbufsize);
  }
} dummy;

} // namespace
