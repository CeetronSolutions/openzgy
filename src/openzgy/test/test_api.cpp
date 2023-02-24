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
do_test_readmeta(const OpenZGY::IZgyTools& r, bool complete, int version)
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
  TEST_EQUAL(r.nlods(), (complete ? 3 : 1));
  TEST_EQUAL_FLOAT(r.datarange()[0], -10038.5, 0.5);
  TEST_EQUAL_FLOAT(r.datarange()[1], +9761.62, 0.5);
  if (complete) {
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
  if (complete) {
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
  do_test_readmeta(*reader, true, 3);
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
  do_test_readmeta(*writer, false, 3);
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
  do_test_readmeta(*reader, true, 1);
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
  TEST_CHECK(buf[0] == 0);
  TEST_CHECK(buf[1] == 1);
  TEST_CHECK(buf[2] == 2);
  TEST_CHECK(buf[8*64+0] == 10);
  TEST_CHECK(buf[8*64+1] == 11);
  TEST_CHECK(buf[8*64+2] == 12);
  TEST_CHECK(buf[64*64+0] == 20);
  TEST_CHECK(buf[64*64+1] == 21);
  TEST_CHECK(buf[64*64+2] == 22);
  reader->close();

  reader = OpenZGY::IZgyReader::open(get_testdata("Empty-v1.zgy"));
  memset(buf.get(), 0, 64*64*64);
  reader->read(orig, size, buf.get(), 0);
  TEST_CHECK(buf[0] == 0);
  TEST_CHECK(buf[1] == 1);
  TEST_CHECK(buf[2] == 2);
  TEST_CHECK(buf[8*64+0] == 10);
  TEST_CHECK(buf[8*64+1] == 11);
  TEST_CHECK(buf[8*64+2] == 12);
  TEST_CHECK(buf[64*64+0] == 20);
  TEST_CHECK(buf[64*64+1] == 21);
  TEST_CHECK(buf[64*64+2] == 22);
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
      TEST_CHECK(reader->nlods() == 2);
      std::pair<bool,double> constval =
        reader->readconst(size3i_t{0,0,0}, size3i_t{64,64,64}, 1, true);
      TEST_CHECK(constval.first);
      TEST_CHECK(constval.second == 42);
      TEST_CHECK(stat.cnt == 128*128*128LL);
      TEST_CHECK(stat.sum == 128*128*128*42.0);
      TEST_CHECK(stat.ssq == 128*128*128*42.0*42.0);
      TEST_CHECK(stat.min == 42.0);
      TEST_CHECK(stat.max == 42.0);
      TEST_CHECK(hist.samplecount == 128*128*128LL);
      //FAILS! TEST_CHECK_(hist.minvalue == 42, "minvalue expect 42 got %lg", hist.minvalue);
      TEST_CHECK_(hist.maxvalue == 42, "maxvalue expect 42 got %lg", hist.maxvalue);
    }
    break;
  case 5:
  case 6:
  case 7: // Lowres data, histogram, statistics not available.
    {
      TEST_CHECK(reader->nlods() == 1);
      must_throw("outside the valid range", [&](){
        reader->readconst(size3i_t{0,0,0}, size3i_t{2,2,2}, 1, true);
      });
      TEST_CHECK(stat.cnt == 0);
      TEST_CHECK(stat.sum == 0);
      TEST_CHECK(stat.ssq == 0);
      TEST_CHECK(stat.min == std::numeric_limits<double>::infinity());
      TEST_CHECK(stat.max == -std::numeric_limits<double>::infinity());
      TEST_CHECK(hist.samplecount == 0);
      TEST_CHECK_(hist.minvalue == 42, "minvalue expect 42 got %lg", hist.minvalue);
      TEST_CHECK_(hist.maxvalue == 42, "maxvalue expect 42 got %lg", hist.maxvalue);
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
  int ii{0};
  for (auto& it : live)
    it = ++ii;
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
  std::string testfolder = InternalZGY::Environment::getStringEnv("OPENZGY_SDTESTDATA", "sd://sntc/testdata");
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
  TEST_CHECK(checkdata[0] == -1000);
  TEST_CHECK(checkdata[63] == 42);
  const SampleHistogram h = reader->histogram();
  TEST_CHECK(similar(h.minvalue, -32768, 1e-5));
  TEST_CHECK(similar(h.maxvalue, +32767, 1e-5));
  TEST_CHECK(h.samplecount == 33*28*92);
  TEST_CHECK(h.bins[124] == 2*3*4);
  TEST_CHECK(h.bins[128] == 33*28*92 - 2*3*4);
  const ZgyWriterArgs::corners_t corners = reader->corners();
  const double eps = 1.0e-10;
  TEST_CHECK(fabs(corners[0][0] -   5) <= eps);
  TEST_CHECK(fabs(corners[0][1] -   7) <= eps);
  TEST_CHECK(fabs(corners[1][0] -   5) <= eps);
  TEST_CHECK(fabs(corners[1][1] - 107) <= eps);
  TEST_CHECK(fabs(corners[2][0] - 205) <= eps);
  TEST_CHECK(fabs(corners[2][1] -   7) <= eps);
  TEST_CHECK(fabs(corners[3][0] - 205) <= eps);
  TEST_CHECK(fabs(corners[3][1] - 107) <= eps);
  TEST_CHECK(reader->size() == (OpenZGY::IZgyWriter::size3i_t{33, 28, 92}));
  TEST_CHECK(reader->datatype() == SampleDataType::int16);
  //TEST_CHECK(reader->hunitname() == "m");
  //TEST_CHECK(reader->zunitname() == "ms");
  TEST_CHECK(reader->annotstart()[0] == 1);
  TEST_CHECK(reader->annotstart()[1] == 500);
  TEST_CHECK(reader->zstart() == 100);
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
 *  \li Compressed files migh (but currently do not) skip the header padding.
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
    must_throw("compressed data is illegal", [&](){
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

  TEST_CHECK(stats->fileVersion() == 3); // v4 only set if actual compressed
  TEST_CHECK(stats->isCompressed() == false); // looked for compressed bricks.
  TEST_CHECK(stats->fileSize() == 4*64*64*64*4); // header plus 3 normal bricks
  TEST_CHECK(stats->brickNormalCount() == 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_CHECK(stats->brickCompressedCount() == 0);
  TEST_CHECK(stats->brickMissingCount() == 0);
  TEST_CHECK(stats->brickConstantCount() == 42);
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
  TEST_CHECK(stats->brickNormalCount() == 0);
  TEST_CHECK(stats->brickCompressedCount() == 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_CHECK(stats->brickMissingCount() == 0);
  TEST_CHECK(stats->brickConstantCount() == 42);
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
  TEST_CHECK(stats->fileSize() == 4*64*64*64*4); // header plus 3 normal bricks
  TEST_CHECK(stats->brickNormalCount() == 3); // 1 lod0, 1 lod1, 1 lod2
  TEST_CHECK(stats->brickCompressedCount() == 0);
  TEST_CHECK(stats->brickMissingCount() == 0);
  TEST_CHECK(stats->brickConstantCount() == 42);
  TEST_CHECK(stats->compressionFactor() > 0.99);
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
  TEST_CHECK(token == token_in_context);
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
  TEST_CHECK(similar(h.minvalue, -1000, 1e-5));
  TEST_CHECK(similar(h.maxvalue, +42, 1e-5));
  TEST_CHECK(h.samplecount == 33*28*92);
  reader->close();
}

static void
test_lod(OpenZGY::DecimationType decimation)
{
  LocalFileAutoDelete lad("testlods.zgy");
  ZgyWriterArgs args = ZgyWriterArgs()
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
  writer->finalize(std::vector<OpenZGY::DecimationType>{decimation}, nullptr);
  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(lad.name());
  std::unique_ptr<float[]> lowres(new float[2*1*51]);
  reader->read(origin, lod1size, lowres.get(), 1);

  // The only input samples that differ from 100 are [0,0,16] and [2,0,32]
  // corresponding to [0,0,8] and [1,0,16] in LOD 1. So the asserts
  // check those two and (for the first value) some samples above and below.

  if (verbose()) {
    std::cout << "LOD 1, algorithm " << int(decimation) << "\n";
    for (int ii=0; ii<16; ++ii)
      std::cout << "  " << lowres[ii] << "\n";
  }

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

  for (int ii=0; ii<16; ++ii)
    TEST_CHECK(similar(lowres[ii], expect[ii], 0.02));

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

  // [3] nothing written.
  // Note that the writer might need to pass force=true to finalize()
  // to get the histogram- and statistics information written out even
  // when no actual data has been written. I am unsure about how the
  // principle of least surprise applies here. As of Oct 2020 the force
  // is required. See the ZgyWriter constructor setting _dirty(false).

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
  TEST_CHECK(r.histo_lo ==  -42 && r.histo_hi ==    0);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[0] == BRICK);
  TEST_CHECK(r.bins[255] == 2*BRICK);

  // [6] single negative value in all three bricks.
  // The value range and the statistics should have the true
  // range i.e. low==high and the histogram range should be wider.
  r = test_histo_onevalue(SampleDataType::float32, -42, true);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==  -42 && r.stats_hi ==  -42);
  TEST_CHECK(r.histo_lo ==  -42 && r.histo_hi ==    0);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[0] == 3*BRICK);
  TEST_CHECK(r.bins[255] == 0);

  // [6] single positive value in all three bricks.
  // Result similar to the above but the ranges are swapped.
  r = test_histo_onevalue(SampleDataType::float32, +42, true);

  TEST_CHECK(r.range_lo == r.stats_lo && r.range_hi == r.stats_hi);
  TEST_CHECK(r.histo_count == r.stats_count);
  TEST_CHECK(r.stats_lo ==   42 && r.stats_hi ==   42);
  TEST_CHECK(r.histo_lo ==    0 && r.histo_hi ==   42);
  TEST_CHECK(r.stats_count == 3*BRICK);
  TEST_CHECK(r.bins[0] == 0);
  TEST_CHECK(r.bins[255] == 3*BRICK);
}

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
  TEST_CHECK(r.stats_lo == r.stats_hi);
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
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 42) < 0.25);
  TEST_CHECK(std::abs(r.stats_lo - 42) > 0.001); // 42.0 not representable.
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
  TEST_CHECK(r.stats_lo == r.stats_hi);
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
  TEST_CHECK(r.stats_lo == r.stats_hi);
  TEST_CHECK(std::abs(r.stats_lo - 42) < 0.25/256);
  TEST_CHECK(std::abs(r.stats_lo - 42) > 0.001/256); // 42.0 not representable.
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
  TEST_CHECK(r.stats_lo == -5 && r.stats_hi == -5);
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
  ZgyWriterArgs args = ZgyWriterArgs()
    .size(size[0], size[1], size[2])
    .bricksize(bs[0], bs[1], bs[2])
    .datatype(SampleDataType::int8)
    .datarange(-128, +127)
    .filename(lad.name());
  std::shared_ptr<OpenZGY::IZgyWriter> writer = IZgyWriter::open(args);
  std::vector<std::int8_t> data(size[0]*size[1]*size[2], 10);
  for (int ii=0; ii<size[0]; ii += 2)
    for (int jj=0; jj<size[1]; jj += 2)
      for (int kk=0; kk<size[2]; kk += 2)
        data[ii*size[1]*size[2] + jj*size[2] + kk] = 90;
  writer->write(size3i_t{0,0,0}, size, data.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{
         OpenZGY::DecimationType::Average
      }, nullptr);
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
  std::shared_ptr<float> data(new float[size[0]*size[1]*size[2]]);
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
 * See ZgyInternalBulk::expeditedRead() and _deliverOneBrick().
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
  std::int64_t pos{0};
  for (std::int64_t ii=0; ii<size[0]; ii += bs[0]) {
    for (std::int64_t jj=0; jj<size[1]; jj += bs[1]) {
      for (std::int64_t kk=0; kk<size[2]; kk += bs[2]) {
        if ((pos % 7) == 0) {
          // missing brick
        }
        else if ((pos % 13) == 0) {
          // Constant-value
          data[0] = static_cast<T>(-(ii/bs[0]) - (jj/bs[1]) - (kk/bs[2]));
          writer->writeconst(size3i_t{ii,jj,kk}, bs, data.data());
        }
        else {
          data[0] = static_cast<T>(ii/bs[0]);
          data[1] = static_cast<T>(jj/bs[1]);
          data[2] = static_cast<T>(kk/bs[2]);
          writer->write(size3i_t{ii,jj,kk}, bs, data.data());
        }
        ++pos;
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
showex(const std::exception& ex)
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
    SeismicStoreIOContext ctxt = SeismicStoreIOContext()
      .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
      .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
      .sdtoken("ServiceSauthV2:" + client_id + ":" + client_secret);
    try {
      std::shared_ptr<OpenZGY::IZgyReader> reader = IZgyReader::open(name, &ctxt);
      if (!TEST_CHECK(bool(reader)))
        return;
      reader->close();
      if (verbose())
        std::cerr << "Success in using client credentials grant.\n";
    }
    catch (const OpenZGY::Errors::ZgyUserError& ex) {
      if (0!=strcmp(ex.what(), "Credentials \"ServiceSauthV2\" not supported")) {
        throw;
      }
      else {
        // This is normal. SDAPI compiled from the public OSDU source.
        // TODO-Low: Might want to verify this from config.sh
        if (verbose())
          std::cerr << "Skipping test. Sauth extension not present.\n";
      }
    }
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
    OpenZGY::IZgyReader::open("sd://sntc/testdata/Synt2.zgy", &ctxt);
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
#ifdef HAVE_SD
    register_sd_test("api.write_cloud",         test_write_cloud);
    register_sd_test("api.write_cloud_mt",      test_write_cloud_mt);
    register_sd_test("api.alturl",              test_alturl);
    register_sd_test("api.idtoken",             test_idtoken);
    register_sd_test("api.sharecred",           test_sharecred);
#endif
    register_test("api.historange",          test_historange);
    register_test("api.lod_lowpass",         test_lod_lowpass);
    register_test("api.lod_weighted",        test_lod_weighted);
    register_test("api.lod_average",         test_lod_average);
    register_test("api.copy",                test_copy);
    register_test("api.enums",               test_enums);
    register_test("api.dummy_compress",      test_dummy_compress);
    register_test("api.histo_cornercase_f",  test_histo_cornercase_float);
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
#ifdef HAVE_SD
    register_sd_test("api.bat_sd_1",            test_bat_sd_1);
    register_sd_test("api.bat_sd_2",            test_bat_sd_2);
    register_sd_test("api.bat_sd_4",            test_bat_sd_4);
    register_sd_test("api.bat_sd_zfp",          test_bat_sd_zfp);
    register_sd_test("api.roflag",              test_roflag);
    register_sd_test("api.client_cred",         test_client_cred);
#if 0 // Not usable as an automated test
    register_sd_test("bug.671969",              test_bug_671969);
#endif
#endif
  }
} dummy;

} // namespace
