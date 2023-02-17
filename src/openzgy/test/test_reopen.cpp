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

/**
 * \file test_reopen.cpp
 *
 * High level tests as in test_api.cpp, using only the publoc API.
 * Tests relating to update of existing files have been moved here.
 */

#include "test_all.h"
#include "test_utils.h"
#include "test_sdutils.h"
#include "../api.h"
#include "../iocontext.h"
//#include "../exception.h"
#include "../impl/environment.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <memory>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>

using namespace OpenZGY;
//using namespace OpenZGY::Formatters;
using Test_Utils::LocalFileAutoDelete;
using Test_Utils::CloudFileAutoDelete;
using Test_Utils::must_throw;

namespace Test_API {
#if 0
}
#endif

namespace {
  // All methods in this file that might be testing SD
  // should explicitly specify that we don't want a
  // read-only file at the end.
#ifdef HAVE_SD
  const OpenZGY::SeismicStoreIOContext*
  default_sd_context_rw()
  {
    using InternalZGY::Environment;
    static OpenZGY::SeismicStoreIOContext instance =
      OpenZGY::SeismicStoreIOContext()
      .setRoAfterWrite(false)
      .forceRoBeforeRead(false)
      .forceRwBeforeWrite(false)
      .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
      .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
      .sdtoken(Environment::getStringEnv("OPENZGY_TOKEN"), "");
    return &instance;
  }

  const OpenZGY::IOContext*
  default_context_rw()
  {
    return default_sd_context_rw();
  }
#else
  const OpenZGY::IOContext*
  default_context_rw()
  {
    return nullptr;
  }
#endif
}







namespace {
  /**
   * Pass a SilentProgress instance to finalize() to get back the "total"
   * field that will be the same for all invocations. I happen to know
   * that this is the expected number of bricks read or written.
   */
  class SilentProgress
  {
  private:
  std::mutex _mutex;
    std::int64_t _done;
    std::int64_t _total;
  public:
    SilentProgress(const SilentProgress&) = delete;
    SilentProgress& operator=(const SilentProgress&) = delete;
    SilentProgress(): _done(0), _total(0) {}
    std::int64_t done() const { return _done; }
    std::int64_t total() const { return _total; }
    void reset() { _done = _total = 0; }
    bool operator()(std::int64_t done, std::int64_t total) {
      std::lock_guard<std::mutex> lk(_mutex);
      _done = std::max(_done, done);
      _total = total;
      return true;
    }
  };
} // anonymous namespace

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

static std::string
get_sdtestdata(const std::string& name)
{
  std::string testfolder = InternalZGY::Environment::getStringEnv("OPENZGY_SDTESTDATA", "sd://sntc/testdata");
  if (testfolder.back() != '/')
    testfolder += "/";
  return testfolder + name;
}

static void
do_write_twice(const std::string& filename, const IOContext* context = nullptr)
{
  ZgyWriterArgs firstargs = ZgyWriterArgs()
    .iocontext(context)
    .filename(filename)
    .size(33, 28, 92)
    .bricksize(64, 64, 64)
    .datatype(SampleDataType::int16)
    .datarange(-32768,+32767);

  ZgyWriterArgs secondargs = ZgyWriterArgs()
    .iocontext(context)
    .filename(filename)
    .ilstart(1).ilinc(2)
    .xlstart(500).xlinc(5)
    .zstart(100).zinc(4)
    //.hunit(UnitDimension::length, "m", 1)
    //.zunit(UnitDimension::time, "ms", 1000)
    .corners(ZgyWriterArgs::corners_t{5,7,5,107,205,7,205,107});

  ZgyWriterArgs thirdargs = ZgyWriterArgs()
    .iocontext(context)
    .filename(filename);

  // Create with basic information only.
  std::shared_ptr<OpenZGY::IZgyWriter> writer =
    OpenZGY::IZgyWriter::open(firstargs);
  writer->finalize(std::vector<OpenZGY::DecimationType>{}, nullptr);
  writer->close();
  writer.reset();

  // Try to re-open it. Should work since there are no data blocks.
  writer = IZgyWriter::reopen(secondargs);

  // Re-open again, re-specifying the basic info (no error unless it changes)
  // and not re-specifying the mutable part (should then be retained).
  // Nope, now it is illegal to set them, period.
  writer->close();
  writer.reset();
  writer = IZgyWriter::reopen(thirdargs);

  // Continue as in test_write(). Probably should have refactored to share code.
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
  std::unique_ptr<float[]> checkdata(new float[64*64*64]);
  const OpenZGY::IZgyWriter::size3i_t origin{0,0,0};
  const OpenZGY::IZgyWriter::size3i_t bsize{64,64,64};
  reader->read(origin, bsize, checkdata.get(), 0);
  TEST_CHECK(checkdata[0] == -1000);
  TEST_CHECK(checkdata[63] == 42);
  const SampleHistogram h = reader->histogram();
  TEST_EQUAL_FLOAT(h.minvalue, -32768, 0.1);
  TEST_EQUAL_FLOAT(h.maxvalue, +32767, 0.1);
  TEST_CHECK(h.samplecount == 33*28*92);
  TEST_CHECK(h.bins[124] == 2*3*4);
  TEST_CHECK(h.bins[128] == 33*28*92 - 2*3*4);
  const ZgyWriterArgs::corners_t corners = reader->corners();
  const double eps = 1.0e-5;
  TEST_EQUAL_FLOAT(corners[0][0],    5, eps);
  TEST_EQUAL_FLOAT(corners[0][1],    7, eps);
  TEST_EQUAL_FLOAT(corners[1][0],    5, eps);
  TEST_EQUAL_FLOAT(corners[1][1],  107, eps);
  TEST_EQUAL_FLOAT(corners[2][0],  205, eps);
  TEST_EQUAL_FLOAT(corners[2][1],    7, eps);
  TEST_EQUAL_FLOAT(corners[3][0],  205, eps);
  TEST_EQUAL_FLOAT(corners[3][1],  107, eps);
  TEST_CHECK(reader->size() == (OpenZGY::IZgyWriter::size3i_t{33, 28, 92}));
  TEST_CHECK(reader->datatype() == SampleDataType::int16);
  //TEST_CHECK(reader->hunitname() == "m");
  //TEST_CHECK(reader->zunitname() == "ms");
  TEST_CHECK(reader->annotstart()[0] == 1);
  TEST_CHECK(reader->annotstart()[1] == 500);
  TEST_CHECK(reader->zstart() == 100);
  reader->close();
}

/**
 * First test of the reopen() feature.
 * Written for a very limited reopen() that simply copied the file.
 * Should still work for the proper implementation, possibly with
 * a few tweaks.
 */
static void
test_reopen_plain()
{
  LocalFileAutoDelete lad("reopen_plain.zgy");
  do_write_twice(lad.name());
  do_check_written(lad.name());
}

#ifdef HAVE_SD

static void
test_reopen_plain_sd()
{
  Test_Utils::CloudFileAutoDelete cad("reopen_plain_sd.zgy", default_sd_context_rw());
  do_write_twice(cad.name(), default_sd_context_rw());
  do_check_written(cad.name(), default_sd_context_rw());
}

#endif // HAVE_SD

enum class TestTwiceFlags : int
{
 nothing         = 0,
 // Actions on first open:
 step1_vt_int8   = 1<<8,    // OFF/on. If on, vt range is ignored. No compress.
 step1_compress  = 1<<0,    // OFF/on.
 step1_write     = 1<<1,    // off/ON. Some +42 samples, a few +41.
 step1_finalize  = 1<<2,    // OFF/ON. Use "Decimate" algorithm.
 step1_keepopen  = 1<<11,   // OFF/on. Finalize but then continue writing.
 // After close and reopen:
 step2_nometa    = 1<<9,    // OFF/on. Second open don't change corners etc.
 step2_compress  = 1<<3,    // OFF/on.
 step2_write     = 1<<4,    // off/ON. Some -10 samples, a few -15.
 step2_rmw       = 1<<5,    // OFF/on. Only with both write, both uncomp.
 step2_replace   = 1<<6,    // OFF/on. All step1_write samples overwritten.
 step2_finalize  = 1<<7,    // off/ON. If on, use "Decimate" algorithm.
 step2_fin_incr  = 1<<10,   // OFF/on. BuildIncremental used.
};

inline TestTwiceFlags
operator&(TestTwiceFlags x, TestTwiceFlags y)
{
  return static_cast<TestTwiceFlags>
    (static_cast<int>(x) & static_cast<int>(y));
}

inline TestTwiceFlags
operator|(TestTwiceFlags x, TestTwiceFlags y)
{
  return static_cast<TestTwiceFlags>
    (static_cast<int>(x) | static_cast<int>(y));
}

namespace TestFormatters {
std::ostream& operator<<(std::ostream& os, TestTwiceFlags value)
{
  typedef TestTwiceFlags F;
  return os
    << (((value & F::step1_vt_int8) != F::nothing) ? "F::step1_vt_int8 | " : "")
    << (((value & F::step1_compress) != F::nothing) ? "F::step1_compress | " : "")
    << (((value & F::step1_write) != F::nothing) ? "F::step1_write | " : "")
    << (((value & F::step1_finalize) != F::nothing) ? "F::step1_finalize | " : "")
    << (((value & F::step1_keepopen) != F::nothing) ? "F::step1_keepopen | " : "")
    << (((value & F::step2_nometa) != F::nothing) ? "F::step2_nometa | " : "")
    << (((value & F::step2_compress) != F::nothing) ? "F::step2_compress | " : "")
    << (((value & F::step2_write) != F::nothing) ? "F::step2_write | " : "")
    << (((value & F::step2_rmw) != F::nothing) ? "F::step2_rmw | " : "")
    << (((value & F::step2_replace) != F::nothing) ? "F::step2_replace | " : "")
    << (((value & F::step2_finalize) != F::nothing) ? "F::step2_finalize | " : "")
    << (((value & F::step2_fin_incr) != F::nothing) ? "F::step2_fin_incr | " : "")
    << "F::nothing";
}
}

/**
 * Write a ZGY file in two parts. Create, close, reopen, close.
 * The TestTwiceFlags bitmask indictes what else to do.
 * There are currently 8 boolean parameters, so in theory 256
 * possible test cases. Only some of these make sense to test.
 * In the list below, upper case means the "preferred" setting
 * of this option.
 *
 * - Not allowed to update a version 1 on-prem file. (Need rawcopy utility)
 * - Not allowed to write a cloud file uploadad by sdutil.
 * - Overwrite might not shrink the value range.
 */
static bool
do_test_reopen(const std::string& filename, TestTwiceFlags flags)
{
  // Size (66,100,100) -> 660,000 samples
  // Bricksize (8,32,64) -> size in bricks 9,4,2 -> nlods 5
  // Valuerange (-100,+155) -> If int8: zero-centric. Ignored for float.
  // Writes will be a combination of:
  // - Full bricks.
  // - Bricks at the survey edge where all the usable samples
  //   and no padding samples are set
  // - Read/modify/write of a brick that has been fully written already.
  //   (Not allowed for compressed files).
  // First write touches all samples with xl 32..63 inclusive.
  // Second write touches all samples with xl 96..99 inclusive.
  // Second write r/m/w updates a small area from (0,32,0)
  // - This only makes sense if step1_write=ON, compression off.
  // Second write replace overwrites all samples from the first write.
  // - Extra credit, ensure this does not trigger a r/m/w.
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  using namespace TestFormatters;

  if (verbose())
    std::cout << "\n  *** test case *** " << flags << "\n";

  // ARRANGE

  const auto flagset = [flags](TestTwiceFlags f) {
                         return (flags & f) != TestTwiceFlags::nothing;
                       };

  constexpr double inf{std::numeric_limits<double>::infinity()};

  // If not finalized on the last close there will be no lowres, stats, histo.
  const int expect_nlods =
    (flagset(TestTwiceFlags::step2_finalize) ? 5 :
     1);

  const bool expect_finalized =
    (flagset(TestTwiceFlags::step2_finalize) ? true :
     false);

  const float expect_lod0first =
    (flagset(TestTwiceFlags::step2_replace) ? -10.0f :
     flagset(TestTwiceFlags::step2_rmw) ? -10.0f :
     flagset(TestTwiceFlags::step1_write) ?  42.0f :
     0);

  const float expect_lod0second =
    (flagset(TestTwiceFlags::step2_write) ? -10.0f :
     0);

  // first lod1 sample is normally identical to lod0 due to using simple decimation
  // unless of course there isn't any low resolution data available.
  const float expect_lod1first = expect_nlods==1 ? -1 : expect_lod0first;
  const float expect_lod1second = expect_nlods==1 ? -1 : expect_lod0second;

  const double expect_coord =
    (flagset(TestTwiceFlags::step2_nometa )|| flagset(TestTwiceFlags::step1_keepopen) ? 0 :
     5);

  const double expect_samplerate =
    (flagset(TestTwiceFlags::step2_nometa) || flagset(TestTwiceFlags::step1_keepopen) ? 6 :
     7);

  // step1_write writes mostly "42" but one in every 7*13 written samples,
  // rounded down, has a special value of "41" and not "42".
  // step2_write, step2_replace, and step2_rmw all write mostly "-10" samples
  // but as above, one in every 7*13 written samples is special. In this case
  // "-15" instead of "-10".
  //
  // Maybe I should have gone for a simpler test pattern?
  //
  // This gave me enough trouble to verify that I need examples:
  //
  // Only step2_rmw:
  //    3200 samples written to upper part of the first region,
  //      35 of those have the special value of -10.
  // Only step2_write:
  //    4000 written close to survey edge, full traces.
  //      43 of those are special.
  // Only step1_write:
  //   32000 written to second region
  //     351 of those are special.
  // step1_write plus write plus step2_rmw write to the same region
  // so the total number of nonzero samples should be the sane as above.
  //   32000 written
  //     351 special, of which
  //      35 should be -10 (same as the rmw-only case)
  //     316 (the rest of the specials) should be 41.
  // Adding step 2 write has no overlap, so just add the step 2 only numbers
  //   36000 (32000 + 4000) written
  //      78 (35 + 43) should be -10
  //     316 (316 + 0) should be 41

  const std::int64_t written_4142_r1 =
    (((flags & TestTwiceFlags::step1_write) == TestTwiceFlags::nothing) ? 0 :
     flagset(TestTwiceFlags::step2_replace) ? 0 :
     flagset(TestTwiceFlags::step2_rmw) ?
     10*32*100 - 1*32*100:
     10*32*100);

  const std::int64_t expect_41 = (written_4142_r1 ) / (7*13);
  const std::int64_t expect_42 = (written_4142_r1 - expect_41);

  const std::int64_t written_1015_r1 =
    (flagset(TestTwiceFlags::step2_write) ? 10*4*100 :
     0);

  const std::int64_t written_1015_r2 =
    (flagset(TestTwiceFlags::step2_replace) ? 10*32*100 :
     flagset(TestTwiceFlags::step2_rmw) ? 1*32*100 :
     0);

  const std::int64_t expect_15 = (written_1015_r1) / (7*13) + (written_1015_r2) / (7*13);
  const std::int64_t expect_10 = (written_1015_r1 + written_1015_r2 - expect_15);

  const std::int64_t expect_stat_cnt = expect_finalized ? 66*100*100 : 0;

  const double expect_stat_sum = !expect_finalized ? 0 :
    (expect_41 * 41.0) +
    (expect_42 * 42.0) +
    (expect_15 * -15.0) +
    (expect_10 * -10.0);

  const double expect_stat_ssq = !expect_finalized ? 0 :
    (expect_41 * 41 * 41.0) +
    (expect_42 * 42 * 42.0) +
    (expect_15 * -15 * -15.0) +
    (expect_10 * -10 * -10.0);

  const double expect_stat_min =  !expect_finalized ? inf :
    (expect_15 + expect_10 == 0 ? 0.0 : -15.0);

  // Note the quirk with incremental finalize. If there were ever
  // some 42's written that number remains in the statistical
  // range even if those samples were removed later.
  // Unless of course the caller asked for an incremental finalize
  // but didn't get it because step 1 did not finalize.
  const double expect_stat_max =  !expect_finalized ? -inf :
    (flagset(TestTwiceFlags::step1_write) &&
     flagset(TestTwiceFlags::step1_finalize) &&
     flagset(TestTwiceFlags::step2_replace) &&
     flagset(TestTwiceFlags::step2_finalize) &&
     flagset(TestTwiceFlags::step2_fin_incr)) ? 42 :
    (expect_41 + expect_42 == 0 ? 0 : 42);

  // Even if asking for incremental build, might not get it.
  // TODO-@@@ Verify by hand the expected block count.
  // Not useful right now because the algorithms are still
  // under development.
  SilentProgress p;
  const bool any_step2_write =
    ((flags & (TestTwiceFlags::step2_write |
               TestTwiceFlags::step2_rmw |
               TestTwiceFlags::step2_replace))
     != TestTwiceFlags::nothing);
  const std::int64_t expect_brickrw =
    !flagset(TestTwiceFlags::step2_finalize) ? 0 :
    !flagset(TestTwiceFlags::step1_finalize) ? 88 :
    !any_step2_write                         ? 0 :
    !flagset(TestTwiceFlags::step2_fin_incr) ? 88 :
    any_step2_write && !flagset(TestTwiceFlags::step1_write) ? 88 : // New: Empty file forces incr off.
    any_step2_write                          ? -21 : // means 1..21
    0;

  ZgyWriterArgs firstargs = ZgyWriterArgs()
    .iocontext(default_context_rw())
    .filename(filename)
    .size(66, 100, 100)
    .bricksize(8, 32, 64)
    .datatype((flags & TestTwiceFlags::step1_vt_int8) != TestTwiceFlags::nothing ?
              SampleDataType::int8 : SampleDataType::float32)
    .datarange(-100,+155)
    .ilstart(101).ilinc(1)
    .xlstart(500).xlinc(5)
    .zstart(100).zinc(6)
    .zunit(UnitDimension::time, "millisec", 0.001)
    .hunit(UnitDimension::length, "feet", 0.3048);

  if (flagset(TestTwiceFlags::step1_compress))
    firstargs.zfp_compressor(99);

  ZgyWriterArgs secondargs = ZgyWriterArgs()
    .iocontext(default_context_rw())
    .filename(filename);

  if ((flags & TestTwiceFlags::step2_nometa) == TestTwiceFlags::nothing)
    secondargs
      .corners(ZgyWriterArgs::corners_t{5,7,5,107,205,7,205,107})
      .zinc(7);

  if (flagset(TestTwiceFlags::step2_compress))
    secondargs.zfp_compressor(99);

  // ACTION
  // FIRST STEPS DONE ON A VIRGIN FILE
  std::shared_ptr<OpenZGY::IZgyWriter> writer =
    OpenZGY::IZgyWriter::open(firstargs);

  if (flagset(TestTwiceFlags::step1_write)) {
    std::vector<float> data(10*32*100, 42);
    for (std::size_t ii = 7*13-1; ii<data.size(); ii += 7*13)
      data[ii] = 41;
    writer->write(size3i_t{0,32,0}, size3i_t{10,32,100}, data.data());
  }

  if (flagset(TestTwiceFlags::step1_finalize)) {
    writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate}, nullptr);
  }

  // CLOSE AND RE-OPEN FILE IF REQUESTED
  if (!flagset(TestTwiceFlags::step1_keepopen)) {
    if (flagset(TestTwiceFlags::step1_finalize)) {
      writer->close();
    }
    else {
      writer->close_incomplete();
    }
    writer.reset();
    writer = IZgyWriter::reopen(secondargs);
  }

  // SECOND STEP IS TO UPDATE THE FILE THAT WAS JUST WRITTEN.

  // Write a region that does not overlap with the first one.
  if (flagset(TestTwiceFlags::step2_write)) {
    std::vector<float> data(10*4*100, -10);
    for (std::size_t ii = 7*13-1; ii<data.size(); ii += 7*13)
      data[ii] = -15;
    writer->write(size3i_t{0,96,0}, size3i_t{10,4,100}, data.data());
  }

  // Overwrite part of the area step 1, meaning there will now be
  // less "41" and "42" samples and those should have been subtracted
  // from the histogram and statistics.
  if (flagset(TestTwiceFlags::step2_rmw)) {
    std::vector<float> data(1*32*100, -10);
    for (std::size_t ii = 7*13-1; ii<data.size(); ii += 7*13)
      data[ii] = -15;
    writer->write(size3i_t{0,32,0}, size3i_t{1,32,100}, data.data());
  }

  // Overwrite the exact same area as in step 1, meaning there will now be
  // no "41" and "42" samples left.
  if (flagset(TestTwiceFlags::step2_replace)) {
    std::vector<float> data(10*32*100, -10);
    for (std::size_t ii = 7*13-1; ii<data.size(); ii += 7*13)
      data[ii] = -15;
    writer->write(size3i_t{0,32,0}, size3i_t{10,32,100}, data.data());
  }

  // Decide whether the end result will have statistics and metadata.
  if (flagset(TestTwiceFlags::step2_finalize)) {
    writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate},
                     std::ref(p),
                     flagset(TestTwiceFlags::step2_fin_incr) ? FinalizeAction::BuildIncremental : FinalizeAction::BuildFull);
    writer->close();
  }
  else {
    writer->close_incomplete();
  }

  writer.reset();

  // OPEN THE FILE AND CHECK RESULTS.
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(filename, default_context_rw());

  const int actual_nlods = reader->nlods();

  const bool actual_has_stats =
    (reader->statistics().cnt > 0 &&
     reader->histogram().samplecount > 0);

  // Read back the first sample of both regions that might or might
  // not have been written. The lod1 sample value should match lod0
  // since the simple "Decimate" algorithm was used.

  float actual_lod0first{-1};
  float actual_lod1first{-1};
  float actual_lod0second{-1};
  float actual_lod1second{-1};
  reader->read(size3i_t{0,32,0}, size3i_t{1,1,1}, &actual_lod0first, 0);
  reader->read(size3i_t{0,96,0}, size3i_t{1,1,1}, &actual_lod0second, 0);
  if (reader->nlods() > 1) {
    reader->read(size3i_t{0,32/2,0}, size3i_t{1,1,1}, &actual_lod1first, 1);
    reader->read(size3i_t{0,96/2,0}, size3i_t{1,1,1}, &actual_lod1second, 1);
  }

  const double actual_coord = reader->corners()[0][0];
  const float actual_samplerate = reader->zinc();

  // Now check all the samples, but only that the count is right,
  std::vector<float> check(10*100*100, -999.25);
  reader->read(size3i_t{0,0,0}, size3i_t{10,100,100}, check.data(), 0);
  std::int64_t actual_41 = std::count_if(check.begin(), check.end(),
                                         [](float x){return x==41;});
  std::int64_t actual_42 = std::count_if(check.begin(), check.end(),
                                         [](float x){return x==42;});
  std::int64_t actual_10 = std::count_if(check.begin(), check.end(),
                                         [](float x){return x==-10;});
  std::int64_t actual_15 = std::count_if(check.begin(), check.end(),
                                         [](float x){return x==-15;});

  SampleStatistics actual_stat(reader->statistics());
  SampleHistogram actual_hist(reader->histogram());

  // Output human readable result.
  if (verbose()) {
    std::stringstream stat_ss;
    stat_ss << "cnt = " << actual_stat.cnt
            << " min = " << actual_stat.min
            << " max = " << actual_stat.max
            << " sum = " << (long long)actual_stat.sum
            << " ssq = " << (long long)actual_stat.ssq;

    std::stringstream hist_ss;
    hist_ss << "cnt = " << actual_hist.samplecount
            << " min = " << actual_hist.minvalue
            << " max = " << actual_hist.maxvalue;
    for (std::size_t ii=0; ii<actual_hist.bins.size(); ++ii)
      if (actual_hist.bins[ii] != 0)
        hist_ss << " bins[" << ii << "]: " << actual_hist.bins[ii];

    std::cout << "\n Read back lod0 " << actual_lod0first
              << " and " << actual_lod0second
              << ", lod1 " << actual_lod1first
              << " and " << actual_lod1second
              << "\n"
              << "  Stat: " << stat_ss.str() << "\n"
              << "  Hist: " << hist_ss.str() << "\n"
              << "  Data: [-15]: " << actual_15
              << " [-10]: " << actual_10
              << " [41]: " << actual_41
              << " [42]: " << actual_42
              << "\n";
  }

    // ASSERT
  bool ok = true;
  ok = TEST_EQUAL(actual_nlods,      expect_nlods)      && ok;
  ok = TEST_EQUAL(actual_has_stats,  expect_finalized)  && ok;
  ok = TEST_EQUAL(actual_lod0first,  expect_lod0first)  && ok;
  ok = TEST_EQUAL(actual_lod1first,  expect_lod1first)  && ok;
  ok = TEST_EQUAL(actual_lod0second, expect_lod0second) && ok;
  ok = TEST_EQUAL(actual_lod1second, expect_lod1second) && ok;
  ok = TEST_EQUAL_FLOAT(actual_coord, expect_coord, 0.001) && ok;
  ok = TEST_EQUAL_FLOAT(actual_samplerate, expect_samplerate, 0.001) && ok;
  ok = TEST_EQUAL(actual_41,         expect_41)         && ok;
  ok = TEST_EQUAL(actual_42,         expect_42)         && ok;
  ok = TEST_EQUAL(actual_10,         expect_10)         && ok;
  ok = TEST_EQUAL(actual_15,         expect_15)         && ok;
  ok = TEST_EQUAL(actual_stat.cnt,   expect_stat_cnt)   && ok;
  ok = TEST_EQUAL_FLOAT(actual_stat.sum,   expect_stat_sum, 0.1)   && ok;
  ok = TEST_EQUAL_FLOAT(actual_stat.ssq,   expect_stat_ssq, 0.1)   && ok;
  ok = TEST_EQUAL(actual_stat.min,   expect_stat_min)   && ok;
  ok = TEST_EQUAL(actual_stat.max,   expect_stat_max)   && ok;
  if (expect_brickrw >= 0) {
    ok = (TEST_EQUAL(p.total(), expect_brickrw)) && ok;
  }
  else {
    // In some cases it is too difficult to manually compute the
    // expected result but it may be possible to set an upper limit.
    ok = (TEST_CHECK(p.total() >= 1)) && ok;
    ok = (TEST_CHECK(p.total() <= std::abs(expect_brickrw))) && ok;
  }
  ok = TEST_EQUAL(p.done(), p.total()) && ok;

  // Not testing histogram min/max because it is rather unclear
  // when it will be too wide (because samples were erased) or too
  // narrow (due to incremental finalize). What I can verify is
  // that because some or all samples were left un-initialized
  // the defaultvalue (i.e. zero) should be included in the range.
  // Probably. Comment in ZgyInternalBulk ctor says maybe not.
  //
  // Not checking histogram samplecount either because if the
  // histogram decided to exclude the defaultvalue then the
  // zeros won't show up in the histogram.
  //
  // Will or *should* defaultvalue be included if all samples have
  // been written? This is outside the scope of this test because
  // that kind of file is not created here.
  //
  // Not testing the actual contents of the histogram because
  // it is unlikely that testing combinations of file update and
  // histogram generation would uncover new issues.
  //
  // Not testing that uncompressed bricks are aligned, compressed
  // are usually not, and the first compressed brick may or may not
  // be aligned (this has not been decided yet).
  //
  // Not testing what happens with alignment if first writes were
  // requested (not necessarily actually) compressed and later were not.
  // Because, who cares? I expect uncompressed to always be aligned
  // but if not then this is just a minor performance issue if caching
  // is ever implemented here.
  //
  // Not testing the opposite case. Not easy to test directly, but
  // I can write a separate test that writes many compressed bricks
  // after a re-open and take a look at the resulting file size.

  return ok;
}

/**
 * Test various things that are illegal, such as trying to change
 * the size or bricksize. Make sure a good error message is produced.
 * Currently the same error is reported for all these cases.
 */
static void
test_reopen_not_if_meta()
{
  LocalFileAutoDelete lad("reopen_not_if_meta.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 100, 100)
     .datarange(-1, 1));
  writer->finalize();
  writer->close();
  writer.reset();

  must_throw("Cannot change", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .size(10, 100, 100));
  });

  must_throw("Cannot change", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .bricksize(64, 64, 64));
  });

  must_throw("Cannot change", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .datatype(SampleDataType::float32));
  });

  must_throw("Cannot change", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .datarange(-1, 1));
  });

  // This restriction might be lifted later.
  must_throw("Cannot change", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .hunit(UnitDimension::length, "m", 1));
  });
}

/**
 * Test "Compressed and cloud files cannot be updated once finalized".
 * The implementation will need an explicit up front test for that,
 * otherwise the error may come too late, after regular bricks have
 * been written. Also the eventual error message would be too cryptic.
 *
 * Beware the special case where the entire file is just one brick,
 * the code may get confused about whether finalize has been run or not.
 */
static void
test_reopen_not_if_final()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_not_if_final.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 100, 300)
     .zfp_compressor(99));
  const float fortytwo{42};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->close();
  writer.reset();

  // Disallow any write of a file that has been finalized with compressed
  // low resolution bricks. Because with current rules you won't be
  // allowed to re-finalize it. That counts as an update.
  // Current checks are incomplete; the error might not be caught
  // until the finalize and at that point the file is effectively corrupt.

  must_throw("finalized compressed file cannot be opened for update", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name()));
  });

  // Re-create the file as uncompressed.
    writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 100, 300));
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->close();
  writer.reset();

  // Disallow opening a finalized file, even if uncompressed, if the
  // second open specifies compression. Because we must assume that
  // the application is going to finalize in this session and that
  // would not be allowed. To be really pedantic I could allow the
  // compressor to be set as long as the lodcompressor is not,
  // but this is getting ridiculous.

  must_throw("finalized file cannot have compressed data appended", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .zfp_compressor(99));
  });
}

/**
 * Test "Compressed and cloud files only allow appending data".
 * In this case it might be sufficient to let the lower levels
 * report the error, unless it is *too* cryptic.
 */
static void
test_reopen_not_if_compress()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_not_if_compress.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 100, 300)
     .zfp_compressor(99));
  const float fortytwo{42}, fifteen{15}, thirtythree{33};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->close_incomplete();
  writer.reset();

  // Ok, append compressed data to a not finalized file.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .zfp_compressor(99));
  writer->write(size3i_t{0,0,64}, size3i_t{1,1,1}, &fifteen);
  writer->close_incomplete();
  writer.reset();

  // Ok but odd, append uncompressed data to a not finalized compressed file.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name()));
  writer->write(size3i_t{0,0,128}, size3i_t{1,1,1}, &fifteen);
  writer->close_incomplete();
  writer.reset();

  // Disallow updating a compressed brick with compressed data.
  // Disallow updating an uncompressed brick with compressed data.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .zfp_compressor(99));
  must_throw("is illegal", [&](){
    writer->write(size3i_t{0,0,64}, size3i_t{1,1,1}, &thirtythree);
  });
  must_throw("is illegal", [&](){
    writer->write(size3i_t{0,0,128}, size3i_t{1,1,1}, &thirtythree);
  });
  writer->close_incomplete();
  writer.reset();

  // In a mixed file, Allow updating if both old and new are uncompressed.
  // Disallow updating a compressed brick with uncompressed data.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name()));
  writer->write(size3i_t{0,0,128}, size3i_t{1,1,1}, &fortytwo);

  // TODO-Minor: This might have failed because the bricks of
  // the old lowres data, being compressed, will be leaked when
  // overwritten by the re-generated and now uncompressed bricks.
  // Mixing compressed and uncompressed this way is discouraged
  // and silly. So this behavior is a quirk, not a bug.
  //
  // Reason: The user is supposed to set a flag telling whether
  // update, r/m/w, etc. is acceptable and indirectly how much disk
  // space leak is acceptable. To simplify things the flag is not
  // exposed; it is set automatically on open. When the file is
  // on-prem and opened uncompressed, the code assumes the old parts
  // of the file are also uncompressed which means there will be no
  // leakage and the flag is set to permissive.
  /*MUST_THROW("is illegal",*/ writer->write(size3i_t{0,0,64}, size3i_t{1,1,1}, &fortytwo);/*);*/
  writer->close_incomplete();
  writer.reset();
}

/**
 * Test that compressed bricks added to an uncompressed file are
 * still written unaligned. Otherwise the compression is pointless.
 */
static void
test_reopen_alignment()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_alignment.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(64, 64, 1024));
  const float fortytwo{42};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->close_incomplete();

  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .zfp_compressor(30));
  for (int ii=64; ii<1024; ii += 64)
    writer->write(size3i_t{0,0,ii}, size3i_t{1,1,1}, &fortytwo);
  writer->close();
  writer.reset();

  std::shared_ptr<OpenZGY::IZgyReader> reader;
  reader = OpenZGY::IZgyReader::open(lad.name(), default_context_rw());

  std::shared_ptr<const FileStatistics> stats = reader->filestats();
  if (verbose())
    stats->dump(std::cout, "   ");

  // 1 MB for headers, 1 MB brick-aligned first uncompressed brick,
  // everything else heavily compressed taking almost no space.
  // Had the compressed bricks erroneously been brick-aligned,
  // maybe because the first brick was, the file would probably
  // have been ~16 MB.
  TEST_CHECK(stats->fileSize() <= 2*1024*1024 + 256*1024);
}

/**
 * Test statistics and histogram in files that contain just one brick.
 * Here, write the entire file in one go. So it isn't technically
 * part of the "reopen" tests.
 */
static void
test_reopen_create_tiny()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_create_tiny.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 20, 30));
  const float fortytwo{42}, sixteen{-16};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->write(size3i_t{0,0,1}, size3i_t{1,1,1}, &sixteen);
  writer->finalize();
  writer->close();
  writer.reset();

  // Make sure statistics were collected.
  std::shared_ptr<OpenZGY::IZgyReader> reader;
  reader = OpenZGY::IZgyReader::open(lad.name(), default_context_rw());
  SampleStatistics stat = reader->statistics();
  TEST_EQUAL(stat.cnt, 10*20*30);
  TEST_EQUAL(stat.min, -16);
  TEST_EQUAL(stat.max, 42);
  TEST_EQUAL(stat.sum, 42 - 16);
}

static float
read_one_sample(const std::string& filename, int ii, int jj, int kk, int lod)
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  std::shared_ptr<OpenZGY::IZgyReader> reader;
  float result{0};
  reader = OpenZGY::IZgyReader::open(filename, default_context_rw());
  reader->read(size3i_t{ii,jj,kk}, size3i_t{1,1,1}, &result, lod);
  return result;
}

/**
 * Test statistics and histogram in files that contain just one brick.
 * These might look like they are already finalized when they are not.
 * Or vice versa, but an extra finalize for those is harmless since
 * no actual bulk will be written.
 */
static void
test_reopen_update_tiny()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_update_tiny.zgy");

  std::shared_ptr<OpenZGY::IZgyWriter> writer;
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(10, 20, 30));
  const float fortytwo{42}, sixteen{-16};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->write(size3i_t{0,0,1}, size3i_t{1,1,1}, &sixteen);
  writer->close_incomplete();
  writer.reset();

  // If this fails then the rest of the tests are pointless.
  if (!TEST_EQUAL(read_one_sample(lad.name(), 0, 0, 0, 0), 42))
    return;

  // Make sure statistics were not collected.
  std::shared_ptr<OpenZGY::IZgyReader> reader;
  reader = OpenZGY::IZgyReader::open(lad.name(), default_context_rw());
  SampleStatistics stat = reader->statistics();
  TEST_EQUAL(stat.cnt, 0);
  reader->close(); // Otherwise it isn't allowed to open for write.

  // Open/close just to get the lowres data correct.
  // The reason this might fail: There are no lowres bricks,
  // so the library might erroneously believe that genlod
  // doesn't need to be run. Or, henlod itself might do an
  // early return. Forgetting that it is also responsible
  // for the statistics.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name()));
  writer->close();

  if (!TEST_EQUAL(read_one_sample(lad.name(), 0, 0, 0, 0), 42))
    return;

  // Make sure statistics are now good.
  reader = OpenZGY::IZgyReader::open(lad.name(), default_context_rw());
  stat = reader->statistics();
  TEST_EQUAL(stat.cnt, 10*20*30);
  TEST_EQUAL(stat.min, -16);
  TEST_EQUAL(stat.max, 42);
  TEST_EQUAL(stat.sum, 42 - 16);
  reader->close();

  // Open/close to erase the statistics.
  // Actually need to write something as well.
  writer = OpenZGY::IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name()));
  writer->write(size3i_t{0,0,2}, size3i_t{1,1,1}, &sixteen);
  writer->close_incomplete();

  // Stale statistics should have been erased.
  reader = OpenZGY::IZgyReader::open(lad.name(), default_context_rw());
  stat = reader->statistics();
  TEST_EQUAL(stat.cnt, 0);

  // Lowres data should still be good.
  TEST_EQUAL(read_one_sample(lad.name(), 0, 0, 0, 0), 42);
}

static bool
copy_file(const std::string& srcname, const std::string& dstname)
{
  std::ifstream src(srcname, std::ios_base::in|std::ios_base::binary);
  if (!TEST_CHECK(src.good()))
    return false;
  std::ofstream dst(dstname, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
  if (!TEST_CHECK(dst.good()))
    return false;
  const std::int64_t bufsize(1*1024*1024);
  std::unique_ptr<char[]> data(new char[bufsize]);
  while (src.good() && dst.good()) {
    src.read(data.get(), bufsize);
    dst.write(data.get(), src.gcount());
  }
  if (!TEST_CHECK(dst.good()))
    return false;
  dst.close();
  src.close();
  return true;
}

/**
 * Version 1 files are not allowed to be updated.
 * Or version 5 or later either, for that matter.
 */
static void
test_reopen_update_too_old()
{
  LocalFileAutoDelete lad("reopen_update_too_old.zgy");
  std::string oldname = get_testdata("Empty-v1.zgy");
  if (!TEST_CHECK(copy_file(oldname, lad.name())))
    return;
  must_throw("version is too old", [&](){
    auto writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name()));
  });
}

/**
 * There is already a check in regular open that compression
 * must have vt float. Make sure that is also checked if it
 * is only the reopen, not the original open, that was
 * uncompressed.
 */
static void
test_reopen_update_bad_vt()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_update_bad_vt.zgy");
  std::shared_ptr<OpenZGY::IZgyWriter> writer;

  // Create uncompressed int8 file.
  writer = OpenZGY::IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .datatype(SampleDataType::int8)
     .datarange(-128, 127)
     .size(10, 100, 100));
  const float fortytwo{42};
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
  writer->close_incomplete();
  writer.reset();

  // Try updating with compressed data.
  must_throw("need to be stored as float", [&](){
    writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .zfp_compressor(30));
  });
  if (writer) {
    // Should not happen. But if it did, see how far we get before it blows up.
    const float sixteen{-16};
    writer->write(size3i_t{0,0,64}, size3i_t{1,1,1}, &sixteen);
    writer->finalize();
    writer->close();
    writer.reset();
  }
}

/*
 * Test not implemented because it is covered by test_reopen.
 *
 * A finalized file will have:
 *
 *  - The single highest-lod brick not flagged as misisng.
 *    This test also works if the file is completely empty,
 *    the brick still changes from misisng to costant.
 *    but it does not work for small single-brick files.
 *
 *  - Statistical sample count != 0.
 *    This test also works if the file is completely empty, because
 *    the zeroes are included in the histogram and statistics. (the
 *    old accessor is inconsistent in that respect, so I get to
 *    choose). This is probably the simplest and most reliable test.
 *    Note that older files might have gotten this wrong, Lowres for
 *    older files will always exist, even for version 1 and compressed
 *    files, but I think stats can be missing. And there is the risk
 *    of some future change that excludes never-written samples from
 *    the statistics. Rather academic question because the empty file
 *    case isn't really that interesting.
 *
 *  - Sum of all histogram bins != 0. Similar to testing the
 *    statistics but my gut feeling is that this is less reliable.
 *    At least when looking at older files.
 *
 *  - If the file is version 1 this implies that lowres exists.
 *
 * When the application closes a file opened for write it chooses
 * whether to finalize (i.e. write statistics abd lowres) the the file
 * or not. Several strategies are possible but nor all are implemented
 * in the API.
 *
 *  1. No, and unconditionally clear any existing statistics and lowres.
 *  2. No, clear existing statistics and lowres only if bulk was changed.
 *  3. Yes, if bulk data was changed or statistics or lowres missing.
 *  4. Yes, unconditionally. Can also used to change lowres algorithms.
 *
 * For best robustness #3 should check both the statistics and the
 * low resolution bricks. It might be acceprable to only check the stats.
 *
 * If the application opens a file just to update metadata then either
 * option #2 or #3 needs to be available.
 *
 * When the application opens a file for read then it needs to know
 * how many LOD levels it can access. If finalize was not run the answer
 * is one (fullres only). Otherwise all lowres levels are avaiable.
 * If the file consists of just a single brick then the answer is 1
 * either way. So it doesn't matter that the test in that case is
 * ambiguous. It is safest to check the highest-lod brick only for
 * this opening-for-read check.
 *
static void
test_is_finalized()
{
}
*/

/*
 * Not implemented because I have not added the step1_noclose flag.
 * The complexity of do_test_reopen() is nearing critical mass.
 * To run a one-off test of this casem tempoarily comment out
 * the close and reopen, then run the main test.
 *
 * Long running processes that might abort and need to be restarted
 * probably want some checkpoint mechanism. The simplest method is to
 * just close and re-open the file. With or without finalizing it
 * depending on whether you want the checkpoint to be fully usable.
 * As of this writing finalize must  process the entire file so you
 * should carefully consider whether you really need to finalize.
 * If the file is on the cloud or if it is compressed you must not
 * finalize the checkpoint.
 *
 * After a crash the writes since the last checkpoint will be lost
 * except when modifying existing bricks in an on-prem file.
 * The data might end up as wasted space on the file.
 *
 * If a crash happens while making the checkpoint the entire file
 * might be corrupted instead of just losing the most recent data.
 * Sorry about that. It is difficult to ensure that the metadata
 * isn't changed before all new bricks have been written.
 * This is on the todo-list.
 *
 * Making a checkpoint for a cloud file that has a large segment size
 * can be considerably less expensive if done immediately after
 * flushing the delay-write buffers. Unfortunately it isn't possible
 * today for the application to know. If this really becomes a
 * problem it is possible to implement some kind of auto-checkpoint
 * mechanism.
 *
 * This tests explores a somewhat less expensive way: call finalize()
 * and just keep using the file instead of a close/re-open. This will
 * only work for uncompressed on-prem because otherwise it isn't allowed
 * to finalize more than once. For cloud access it won't work anyway
 * because you need to close the file in order to get the delayed-write
 * buffers flushed.
 *
static void
test_checkpoint_finalize()
{
  typedef TestTwiceFlags F;
  LocalFileAutoDelete lad("twice.zgy");
  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_write |
                           F::step1_finalize |
                           F::step1_noclose |
                           F::step2_write |
                           F::step2_rmw |
                           F::step2_finalize));
}
*/

/*
 * As test_checkpoint_finalize() but instead of finalizing just call
 * the (not yet implemented) sync() method to flush the metadata
 * that has pointers to the newly added bricks.
 *
 * I doubt it is worth while to implement this. On-prem the close
 * and re-open is fairly cheap and the reopened file has few
 * restrictions on what you can use it for. On the cloud I believe
 * that implementing a sync() would be a huge amount of work and
 * give very little gain.
 *
static void
test_checkpoint_nofinalize()
{
}
*/

/**
 * A file that has no written data at all will still have statistics and
 * low resolution data. There might be snags due to the contents not
 * being flagged as dirty. The test case might actually be covered by
 * one of the tests in test_reopen() but the one here is I believe
 * slightly more pedantic.
 */
static void
test_reopen_empty()
{
  LocalFileAutoDelete lad("reopen_empty.zgy");

  {
    // Create empty file with 3 lod levels
    auto writer = OpenZGY::IZgyWriter::open
        (ZgyWriterArgs()
         .iocontext(default_context_rw())
         .filename(lad.name())
         .size(10, 100, 200));
    writer->close();
  }

  // Check that lowres and statistics were generted,
  // even though the file is not "dirty" in the normal sense.
  {
    auto reader = OpenZGY::IZgyReader::open
      (lad.name(), default_context_rw());
    TEST_EQUAL(reader->nlods(), 3);
    TEST_EQUAL(reader->statistics().cnt, 10*100*200);
  }

  // Forcibly remove the lowres and statistics.
  // Also test a recent  "principle of least surprise" change.
  // After en explicit finalize saying no lowres needed,
  // close() will honor that wish. You don't need to
  // remember to call close_incomplete().
  {
    auto writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name()));
    writer->finalize(std::vector<OpenZGY::DecimationType>{}, nullptr,
                     FinalizeAction::Delete, true);
    writer->close();
  }

  // Check that lowres and statistics were deleted
  {
    auto reader = OpenZGY::IZgyReader::open
      (lad.name(), default_context_rw());
    TEST_EQUAL(reader->nlods(), 1);
    TEST_EQUAL(reader->statistics().cnt, 0);
  }

  // Open/close with no bulk data change.
  {
    auto writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name()));
    writer->close();
  }

  // Even though no bulk was written, a default close should make sure
  // the low resolution data exists.
  {
    auto reader = OpenZGY::IZgyReader::open
      (lad.name(), default_context_rw());
    TEST_EQUAL(reader->nlods(), 3);
    TEST_EQUAL(reader->statistics().cnt, 10*100*200);
  }
}

/**
 * Meant to run manually after enabling some debug logging in
 * ZgyWriter::_create_serived() but has some use even without the
 * manual verification. The test focuses on how the library keeps
 * track of dirty bricks.
 */
static void
test_reopen_track_changes()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_track_changes.zgy");

  // Survey is slightly larger than 16x16 bricks so it triggers
  // the close-to-edge logic. Vertically it is much smaller
  // than horizontally which might also cause issues.
  // Brick size is 16 x 32 x 32 just to make things interesting,
  // inline: 16 full bricks @ 16 lines/brick plus 5 lines = 261
  // xline:  17 full bricks @ 32 lines/brick plus 3 lines = 547
  // time:    6 full bricks @ 32 samples/brick plus 10 = 202
  //
  //   LOD0: 17 x 18 x 7   261 x 547 x 202
  //   LOD1:  9 x  9 x 4   131 x 274 x 101
  //   LOD2:  5 x  5 x 2    66 x 137 x  51
  //   LOD3:  3 x  3 x 1    33 x  69 x  26
  //   LOD4:  2 x  2 x 1    17 x  35 x  13
  //   LOD5:  1 x  1 x 1     9 x  18 x   7
  //
  // Sum 2142 fullres and 388 lowres = 2530 total.

  {
    // Create empty file with 6 lod levels
    SilentProgress p;
    auto writer = OpenZGY::IZgyWriter::open
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name())
       .bricksize(16,32,32)
       .size(261, 547, 202));
    if (!TEST_CHECK(writer != nullptr))
      return;
    const float fortytwo{42};
    writer->writeconst(size3i_t{0,0,0}, writer->size(), &fortytwo);
    // The code now detects all-constant and takes it as a hint to force
    // a full a build. In most cases that will be faster, and very often
    // it prevents a bad histogram. Disable that trick by forcing
    // the first brick to become allocated.
    const float sixteen{16};
    writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &sixteen);
    writer->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &fortytwo);
    writer->finalize(std::vector<OpenZGY::DecimationType>
                     {OpenZGY::DecimationType::Average}, std::ref(p));
    writer->close();
    TEST_EQUAL(p.total(), 2142 + 388);
    TEST_EQUAL(p.done(), p.total());
  }

  {
    SilentProgress p;
    auto writer = OpenZGY::IZgyWriter::reopen
      (ZgyWriterArgs()
       .iocontext(default_context_rw())
       .filename(lad.name()));
    if (!TEST_CHECK(writer != nullptr))
      return;
    const float thirteen{13}, seven{7};
    // One single brick is dirty.
    writer->writeconst(size3i_t{0,32,64}, size3i_t{16,32,32}, &thirteen);
    // Four brick-columns of 7 bricks,  because writes are not aligned
    writer->writeconst(size3i_t{241,513,0}, size3i_t{16,32,202}, &seven);
    // A total of 29 bricks are dirty, and 5 brick-columns at LOD0
    // will need to be re-read assuming the simpler genlod algorithm.
    writer->finalize(std::vector<OpenZGY::DecimationType>
                     {OpenZGY::DecimationType::Average},
                     std::ref(p), FinalizeAction::BuildIncremental);
    writer->close();
    TEST_EQUAL(p.total(), 97); // difficult to verify by hand.
    TEST_EQUAL(p.done(), p.total());
  }

  {
    auto reader = OpenZGY::IZgyReader::open
      (lad.name(), default_context_rw());
    if (!TEST_CHECK(reader != nullptr))
      return;
    TEST_EQUAL(reader->nlods(), 6);
    TEST_EQUAL(reader->statistics().cnt, 261*547*202);
    TEST_EQUAL(reader->statistics().sum,
               42*(261*547*202 - 16*32*32 - 16*32*202) +
               13*(16*32*32) + 7*(16*32*202));
    // Check lod1 to verify that at least the dirty data caused lowres
    // to be re-calculated. To also verify that the code din't do more
    // than it needed then some logging needs to be turned on.
    float data1{0};
    reader->read(size3i_t{0,16,32}, size3i_t{1,1,1}, &data1, /*lod=*/1);
    TEST_EQUAL_FLOAT(data1, 13.0, 0.001);

    float data2{0};
    // lod0:(241,513) maps to lod1:(120,256) but so does the 3 neighbor
    // brick columns which were not overwritten by seven.
    reader->read(size3i_t{120,256,0}, size3i_t{1,1,1}, &data2, /*lod=*/1);
    TEST_EQUAL_FLOAT(data2, (2*7.0 + 6*42.0)/8, 0.001);
  }
}

/**
 * This is one of the test cases run in test_reopen() but I'll also
 * run it explicitly here. Because it needs some explanation.
 *
 * Outside of unit tests this corner case is very unlikely to show up.
 *
 * After the first step the histogram range is [41..42], it excludes
 * defaultvalue (zero) and the value of the samples to be added later.
 * After the second step the file contains only -15, -10, and 0
 * samples meaning that he histogram, which cannot change its limits
 * in an incremental rebuild, ends up empty.
 *
 * Technically the empty histogram is correct. The API warns that the
 * histogram might not include all addded samples. So the verification
 * step ought to accept that the histogram might end up empty.
 *
 * The code in ZgyWriter::_create_derived() should recognizes this
 * situation and switch to a full rebuild just to avoid the confusion.
 * If that is changed then the asserts in the test need to be relaxed.
 */
static void
test_reopen_bad_histogram()
{
  typedef TestTwiceFlags F;
  LocalFileAutoDelete lad("reopen_bad_histogram.zgy");
  TEST_CHECK(do_test_reopen(lad.name(),
                            F::step1_write |
                            F::step1_finalize |
                            F::step2_replace |
                            F::step2_finalize));
}

/**
 * Extensive test of the reopen() feature.
 */
static void
test_reopen()
{
  typedef TestTwiceFlags F;
  LocalFileAutoDelete lad("reopen.zgy");

  // In theory, all possible flag combinations yield 2048 test cases.
  //
  // step1_noclose is not called from this function. Which leaves 1024.
  //
  // Nearly all test runs can turn off step1_vt_int8 and step2_nometa.
  // I will run a single test case with both those. I am reasonably
  // sure that non-float cubes should work similar to float.
  // Except there are some subtle changes with the histogram range.
  //
  // Most test cases want step2_finalize. Run just two test cases with
  // that one turned off: with and without step1_finalize.
  //
  // No negative tests here; I'll do those in separate functions.

  // Run the two most useful tests first.
  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw |
                           F::step2_finalize));

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_compress |
                           F::step2_compress |
                           F::step1_write |
                           F::step2_write |
                           F::step2_finalize));

  // Three special cases, explained above, that are not included
  // in the main loop because they have combinations of
  // step1_vt_int8, step2_nometa, and step2_finalize that I
  // don't want to test too many variations of.

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw));      // No finalize at all.

  // If the file is on seismic store, maybe just test the above three
  // cases for performance reasons. Or cherry-pick a few extra.
  // In the cloud case I might also need to test a file large enough
  // to need more than segsplit data segments.

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw |
                           F::step1_finalize)); // No step2_finalize.

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_vt_int8 |   // Only tested here.
                           F::step2_nometa |    // Only tested here.
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw |
                           F::step2_finalize));

  // The first 7 flag bits are those that most often change,
  // so I'll run all 128 combinations of those.
  // Except that turning on any of the compression options,
  // i.e. 96 of the 128 possble combinations, have limitations
  // on what the other flags can be. So there will be 32 tests
  // for uncompressed and 12 for compressed.
  int runs = 5; // for the tests above.
  for (int ii=0; ii<128; ++ii) {
    TestTwiceFlags flags = static_cast<TestTwiceFlags>(ii);
    flags = flags | TestTwiceFlags::step2_finalize;
    // Only append is allowed for compressed files. For that reason a double
    // finalize is also bad because it would need to replaces lowres bricks.
    if ((flags & (F::step1_compress | F::step2_compress)) != F::nothing &&
        (flags & (F::step2_rmw | F::step2_replace | F::step1_finalize)) != F::nothing)
      continue;
    try {
      ++runs;
      TEST_CHECK_(do_test_reopen(lad.name(), flags), "case 0x%02x", ii);
    }
    catch(const std::exception& ex) {
      TEST_CHECK_(false, "case 0x%02x raised %s", ii, ex.what());
    }
  }
  // Just in case I messed up the "valid" test.
  TEST_EQUAL(runs, 5+32+12);
}

/**
 * Extensive test of the reopen() feature.
 * Requesting incremental build of derived data.
 */
static void
test_reopen_incr()
{
  typedef TestTwiceFlags F;
  LocalFileAutoDelete lad("reopen_incr.zgy");

  // Run roughly the same tests as in test_reopen().
  // Skip those that don't have both finalize step1 and step2 turned on,
  // and skip compressed files because incremental is not allowed there.
  // Paranoia: Ideally I'd like to test a few int16_t cubes here as well
  // because histograms for those work slightly differently.

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw |
                           F::step2_fin_incr | F::step2_finalize));

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_compress |
                           F::step2_compress |
                           F::step1_write |
                           F::step2_write |
                           F::step2_fin_incr |F::step2_finalize));

  TEST_CHECK(do_test_reopen(lad.name(),
                           F::step1_vt_int8 |   // Only tested here.
                           F::step2_nometa |    // Only tested here.
                           F::step1_write |
                           F::step2_write |
                           F::step2_rmw |
                           F::step2_fin_incr | F::step2_finalize));

  // here will be 32 tests for uncompressed and 12 for compressed.
  int runs = 3; // for the tests above.
  for (int ii=0; ii<128; ++ii) {
    TestTwiceFlags flags = static_cast<TestTwiceFlags>(ii);
    flags = flags | TestTwiceFlags::step2_fin_incr | TestTwiceFlags::step2_finalize;
    // Only append is allowed for compressed files. For that reason a double
    // finalize is also bad because it would need to replaces lowres bricks.
    if ((flags & (F::step1_compress | F::step2_compress)) != F::nothing)
      continue;
    try {
      ++runs;
      TEST_CHECK_(do_test_reopen(lad.name(), flags), "case 0x%02x", ii);
    }
    catch(const std::exception& ex) {
      TEST_CHECK_(false, "case 0x%02x raised %s", ii, ex.what());
    }
  }
  // Just in case I messed up the "valid" test.
  TEST_EQUAL(runs, 3+32);
}

/**
 * Similar to the "reopen" case but instead if close and re-open
 * the code just does a finalize and then continues writing.
 * This may be useful for applications that keep a file open
 * for a long time, updating it frequently and assuming that
 * other parts of the application can read low resolution
 * data and assume it is kept up to date.
 *
 * This set of tests turn on the keepopen flag which is only useful
 * when both calls to finalize are turned on and an incremental
 * second finalize is requested. And might as well only test the
 * case with the most "interesting" writes i.e. everything except
 * step2_replace. The nometa flag is implied because changing the
 * metadata was to be done by thereopen. Compression is not allowed
 * (because the test does multiple finalize) which means there aren't
 * really that many combinations left to test. In fact, just one.
 */
static void
test_reopen_keepopen()
{
  typedef TestTwiceFlags F;
  LocalFileAutoDelete lad("reopen_keepopen.zgy");

  TEST_CHECK(do_test_reopen(lad.name(),
                            F::step1_write     |
                            F::step1_finalize  |
                            F::step1_keepopen  |
                            F::step2_write     |
                            F::step2_rmw       |
                            F::step2_finalize  |
                            F::step2_fin_incr));
  // do_test_reopen() verifies that the second finalize needs less I/O
  // than the 88 operations a full one needs, but it doesn't check the
  // exact number which means it doesn't check that with and without
  // keepopen needs the exact same bulk I/O. Maybe not important.
}

/**
 * Similar to test_reopen_keepopen() but here the focus is on reading
 * low resolution data while the file is still open for write.
 * That is probably the only reason the application would call
 * finalize and continue reading.
 *
 * Ideally test both the case with a single open/create and the case
 * where the file is created, written, closed, reopened, written more,
 * and than read lodres data. Currently do just the simple one.
 *
 * TODO-Future-@@@: The library might change to detect that the
 * application is asking for stale data and recompute and store
 * the data transparently to the applicaton. If this gets implemented
 * then that feature would also be tested here.
 *
 * Minor: Also test that with two calls to finalize in a row there
 * should be no bulk data I/O at all in the second finalize.
 */
static void
test_reopen_rwlod()
{
  LocalFileAutoDelete lad("reopen_rwlod.zgy");
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  using namespace OpenZGY;
  SilentProgress p;

  std::shared_ptr<IZgyWriter> writer =
    IZgyWriter::open
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name())
     .size(66, 100, 100)
     .bricksize(8, 32, 64)
     .datarange(-100,+155));

  std::vector<float> data80(10*32*100, 80);
  writer->write(size3i_t{0,32,0}, size3i_t{10,32,100}, data80.data());
  TEST_EQUAL(writer->nlods(), 1); // Not yet finalized.
  writer->finalize(std::vector<DecimationType>{DecimationType::Average},
                   std::ref(p), FinalizeAction::BuildFull);
  TEST_EQUAL(p.total(), 88);
  p.reset();
  writer->close();
  writer = IZgyWriter::reopen
    (ZgyWriterArgs()
     .iocontext(default_context_rw())
     .filename(lad.name()));

  const float ten{10};
  writer->write(size3i_t{0,32,1}, size3i_t{1,1,1}, &ten);

  writer->finalize(std::vector<DecimationType>{DecimationType::Average},
                   std::ref(p), FinalizeAction::BuildIncremental);
  TEST_EQUAL(p.total(), 14);
  p.reset();
  writer->finalize(std::vector<DecimationType>{DecimationType::Average},
                   std::ref(p), FinalizeAction::BuildIncremental);
  TEST_EQUAL(p.total(), 0);
  TEST_EQUAL(writer->nlods(), 5);

  // Check that the low resolution sample reflects the change done
  // after the first finalize. The API currently won't allow this.
  // read() and readconst() don't have a "lod" parameters.
  // That can easily be added but due to all the caveats, such
  // as taking care of periodic finalize, it is better to keep
  // the feature hidden until it is likely that it will be used.
#if 0 // TODO-@@@: Allow reading lowres data from file open for write.
  float check{-1};
  writer->read(size3i_t{0,32/2,0}, size3i_t{1,1,1}, &check, 1);
  TEST_EQUAL(check, float(71.25));
#endif
}

/**
 * The short version: When OpenZGY saves a file that isn't finalized
 * then it is flagged as v4 because the old reader cannot handle that.
 * Just like compressed data will do. Also there will be no
 * incremental finalize on the next close since there is no starting
 * point. On the positive side, compressed files that are not
 * finalized can still be re-opened and appended to.
 *
 * The long version: The above is complicated by the fact that
 * finalize does more than one thing so there exists a "halfway
 * finalized" state where low resolution data is present and
 * statistics are not. Or vice versa. Or low resolution is missing
 * because it is not needed for tiny files.
 *
 * The version number (3 or 4) should reflect the current state of the
 * file, not its history. So e.g. a file that was v4 only due to
 * missing low resolution data should be changed to v3 when finalized.
 * The tests that determine the version number are similar to the tests
 * that decide whether incremental builds are allowed and whether a
 * compressed file may be re-opened.
 *
 * Why does this matter?
 *
 * Petrel creates an empty file that will be written to later. This
 * file must not be finalized because that would prevent it from being
 * later written with compressed data. And if the finalize was done
 * with compression it would prevent even ubcompressed writes. If no
 * compression at all is involved then a finalize at this point "only"
 * a serious performance issue.
 *
 * So:
 * - Create empty file not finalized. It will be version 4.
 *   Except the special single-brick case which will be v3
 *   because there isn't any lowres to be generated.
 *   It doesn't matter whether compression is requested;
 *   the file will be uncompressed because it has no data.
 *
 * Features:
 * (ZFP) Has compressed data.
 * (LOD) Has low resolution data (N/A for single brick file).
 * (S/H) Has statistics and histogram (empty histogram due to wrong range fails)
 *
 * Behavior depending on the above:
 * - Mark as V4 if ZFP || !LOD.
 * - Not allowed to reopen if ZFP && LOD.
 * - Track changes for BuildIncremental if LOD && S/H
 *
 * Behavior for single-brick file which cannot have LOD:
 * - Mark as V4 if ZFP
 * - Not allowed to reopen if ZFP.
 * - Never track changes.
 *
 * Notes:
 * - Applications using the old readers are assumed to handle missing
 *   statistics and histogram correctly, so there is no need to flag
 *   as v4 only for that purpose.
 *
 * - A single brick file cannot have lowres, this lack will not trigger
 *   the v4 flag either.
 *
 * - Allowing BuildIncremental has a stricter check. Err on the size
 *   of forcing a BuildFull where an incremental might have sufficed.
 *   Here, missing S/H makes us assume the file is not finalized yet.
 *   A single brick file always fails this test because incremental
 *   builds would be pointless. (test: See if the statistics can shrink).
 *
 * - Reopen of a compressed file is not allowed if finalized because the
 *   low resolution bricks (just like the full resolution) can only be
 *   written once. Reopen of a compressed single brick file is allowed
 *   for consistency even though it would be pointless. Not due to the
 *   double finalize but because if the single data brick is already
 *   written then there are no more fullres bricks to write. And if it
 *   is not written then the file cannot be compressed.
 */
static void
test_reopen_setversion()
{
  LocalFileAutoDelete lad("reopen_setversion.zgy");
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  using namespace OpenZGY;
  constexpr double inf{std::numeric_limits<double>::infinity()};

  {
    // Empty file finalized, request compression
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                                .filename(lad.name())
                                .zfp_compressor(99)
                                .size(33, 128, 92));
    writer->close();
  }

  {
    // In spite of compression being requested there wasn't actually
    // any compressed data written. So the version should still be 3.
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 2);
    TEST_EQUAL(reader->filestats()->fileVersion(), 3);
    TEST_EQUAL(reader->statistics().cnt, 33*128*92);
    TEST_EQUAL(reader->statistics().min, 0);
    TEST_EQUAL(reader->statistics().max, 0);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
    // Histogram range is unspecified when there is just a single value.
    //TEST_EQUAL(reader->histogram().minvalue,    reader->statistics().min);
    //TEST_EQUAL(reader->histogram().maxvalue,    reader->statistics().max);
  }

  {
    // NEW -> Empty file not finalized
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                                .filename(lad.name())
                                .size(33, 128, 92));
    writer->close_incomplete();
  }

  {
    // Check 1 lod i.e. only fullres, and v4 for that reason.
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 1);
    TEST_EQUAL(reader->filestats()->fileVersion(), 4);
    TEST_EQUAL(reader->statistics().cnt, 0);
    TEST_EQUAL(reader->statistics().min, inf);
    TEST_EQUAL(reader->statistics().max, -inf);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
  }

  {
    // REOPEN -> Empty file finalized
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(lad.name()));
    writer->close();
  }

  {
    // Should now be ZGY-Public compatible. Stats show only zeros.
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 2);
    TEST_EQUAL(reader->filestats()->fileVersion(), 3);
    TEST_EQUAL(reader->statistics().cnt, 33*128*92);
    TEST_EQUAL(reader->statistics().min, 0);
    TEST_EQUAL(reader->statistics().max, 0);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
  }

  {
    // NEW -> Start over, make a file with compressed data not finalized.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                                .filename(lad.name())
                                .zfp_compressor(99)
                                .size(33, 128, 92));
    std::vector<float> data(64*64*64, 1000.0f);
    data[0] = 42.0f;
    writer->write(size3i_t{0,0,0}, size3i_t{64,64,64}, data.data());
    writer->close_incomplete();
  }

  {
    // Check 1 lod i.e. only fullres, and v4 for two reasons.
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 1);
    TEST_EQUAL(reader->filestats()->fileVersion(), 4);
    TEST_EQUAL(reader->statistics().cnt, 0);
    TEST_EQUAL(reader->statistics().min, inf);
    TEST_EQUAL(reader->statistics().max, -inf);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
    // Histogram range might be the value ranges seen till now,
    // but this is an implementation detail since the histogram
    // is still empty.
    TEST_EQUAL(reader->histogram().minvalue, 42.0);
    TEST_EQUAL(reader->histogram().maxvalue, 1000.0);
  }

  {
    // REOPEN -> Write more compressed data and then finalize.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(lad.name())
                                  .zfp_compressor(99));
    const float fifteen{1500.0f};
    writer->write(size3i_t{0,64,0}, size3i_t{1,1,1}, &fifteen);
    writer->close();
  }

  {
    // Check all lods, still v4 because of the compression.
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 2);
    TEST_EQUAL(reader->filestats()->fileVersion(), 4);
    TEST_EQUAL(reader->statistics().cnt, 33*128*92);
    TEST_EQUAL(reader->statistics().min, 0.0);
    TEST_EQUAL(reader->statistics().max, 1500.0);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
    TEST_EQUAL(reader->histogram().minvalue,    reader->statistics().min);
    TEST_EQUAL(reader->histogram().maxvalue,    reader->statistics().max);
  }

  {
    // NEW -> Start over, tiny file (one brick), written not finalized.
    // Not finalized means the statistics and histogram are unset.
    // Low resolution data is N/A. For the v4 test the lowres is treated
    // as if present i.e. not v4.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                                .filename(lad.name())
                                .size(33, 22, 11));
    const float fortytwo{42.0f};
    writer->write(size3i_t{1,2,3}, size3i_t{1,1,1}, &fortytwo);
    writer->close_incomplete();
  }

  {
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 1);
    TEST_EQUAL(reader->filestats()->fileVersion(), 3);
    TEST_EQUAL(reader->statistics().cnt, 0);
    TEST_EQUAL(reader->statistics().min, inf);
    TEST_EQUAL(reader->statistics().max, -inf);
    TEST_EQUAL(reader->histogram().samplecount, 0);
    // Histogram range might be the value ranges seen till now,
    // but this is an implementation detail since the histogram
    // is still empty.
    TEST_EQUAL(reader->histogram().minvalue, 0.0);
    TEST_EQUAL(reader->histogram().maxvalue, 42.0);
  }

  {
    // REOPEN -> Finalize the tiny file.
    // One twist: Explicitly request an incremental build.
    // The code should switch to a full build since there
    // was no previous finalize but that might not be
    // obvious at this point.
    // Another thing is that incremental rebuild makes
    // absolutely no sense in this case.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(lad.name()));
    writer->finalize(std::vector<OpenZGY::DecimationType>
                     {OpenZGY::DecimationType::Average},
                     nullptr, FinalizeAction::BuildIncremental);
    writer->close();
  }

  {
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    TEST_EQUAL(reader->nlods(), 1);
    TEST_EQUAL(reader->filestats()->fileVersion(), 3);
    // Off topic: This also verifies that statistics and histogram were
    // created even though genlod() didn't actually have anything to do.
    TEST_EQUAL(reader->statistics().cnt, 33*22*11);
    TEST_EQUAL(reader->statistics().min, 0);
    TEST_EQUAL(reader->statistics().max, 42);
    TEST_EQUAL(reader->histogram().samplecount, reader->statistics().cnt);
    TEST_EQUAL(reader->histogram().minvalue,    reader->statistics().min);
    TEST_EQUAL(reader->histogram().maxvalue,    reader->statistics().max);

  }

  {
    // NEW -> Start over, as above but compressed data. tiny file (one brick),
    // written compressed not finalized. Statistics and histogram are unset.
    // Low resolution data is N/A. For the v4 test the lowres is treated
    // as if present i.e. not v4. For compressed re-open lowres is
    // treated as not present i.e. reopen ok.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                                .filename(lad.name())
                                .zfp_compressor(99)
                                .size(33, 22, 11));
    const float fortytwo{42.0f};
    writer->write(size3i_t{1,2,3}, size3i_t{1,1,1}, &fortytwo);
    writer->close_incomplete();
  }

  {
    // REOPEN -> Finalize the compressed tiny file with compressed lowres.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                  .zfp_compressor(99)
                                  .filename(lad.name()));
    writer->finalize(std::vector<OpenZGY::DecimationType>
                     {OpenZGY::DecimationType::Average},
                     nullptr, FinalizeAction::BuildIncremental);
    writer->close();
  }

#if 0 // Don't really care either way. For consistency maybe throw,
  {
    // REOPEN -> Finalize the compressed tiny file one more time.
    // Normally this is not allowed, but since the finalize does not
    // actually store any data in this case it should be ok.
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                  .zfp_compressor(99)
                                  .filename(lad.name()));
    writer->finalize(std::vector<OpenZGY::DecimationType>
                     {OpenZGY::DecimationType::Average},
                     nullptr, FinalizeAction::BuildIncremental);
    writer->close();
  }
#endif
}

/**
 * Opening an empty file created by the old ZGY accessor has some
 * challenges with respect to alignment.
 */
static void
test_reopen_zgypublic()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  LocalFileAutoDelete lad("reopen_zgypublic.zgy");
  std::string oldname = get_testdata("EmptyOldFile.zgy");
  if (!TEST_CHECK(copy_file(oldname, lad.name())))
    return;

  ZgyWriterArgs secondargs = ZgyWriterArgs().filename(lad.name());

  {
    if (verbose())
      std::cout << "Append, first brick all 1000\n";
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(secondargs);
    std::vector<float> data(64*64*64, 1000.0f);
    std::fill(data.data() + 64*64*64 - 64, data.data() + 64*64*64, 0.0f);
    writer->write(size3i_t{0,0,0}, size3i_t{64,64,64}, data.data());
    writer->close();
  }

  {
    if (verbose())
      std::cout << "Append, write second brick all 2000\n";
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(secondargs);
    // Will not be a round number if the code failed to add padding.
    TEST_EQUAL(writer->filestats()->fileSize(), 64*64*64*16);
    std::vector<float> data(64*64*64, 2000);
    writer->write(size3i_t{64,64,64}, size3i_t{64,64,64}, data.data());
    writer->close();
  }

  {
    if (verbose())
      std::cout << "Overlap, write 3000 in some samples\n";
    std::shared_ptr<OpenZGY::IZgyWriter> writer =
      OpenZGY::IZgyWriter::reopen(secondargs);
    std::vector<float> data(65*65*65, 3000);
    writer->write(size3i_t{96,96,96}, size3i_t{65,65,65}, data.data());
    writer->close();
  }

  static const auto inside = [](int ii, int jj, int kk, const size3i_t& start, const size3i_t& size) {
    return (ii >= start[0] && ii < start[0] + size[0] &&
            jj >= start[1] && jj < start[1] + size[1] &&
            kk >= start[2] && kk < start[2] + size[2]);
  };
  static const auto expect = [](int ii, int jj, int kk) -> float {
    return (inside(ii, jj, kk, size3i_t{96,96,96}, size3i_t{65,65,65}) ?
            3000.0f :
            inside(ii, jj, kk, size3i_t{64,64,64}, size3i_t{64,64,64}) ?
            2000.0f :
            inside(ii, jj, kk, size3i_t{63,63,0}, size3i_t{1,1,64}) ?
            0.0f :
            inside(ii, jj, kk, size3i_t{0,0,0}, size3i_t{64,64,64}) ?
            1000.0f :
            0.0f);
    };

  {
    if (verbose())
      std::cout << "Read back and check\n";
    std::shared_ptr<OpenZGY::IZgyReader> reader =
      OpenZGY::IZgyReader::open(lad.name());
    std::vector<float> check(192*192*192, -999);
    reader->read(size3i_t{0,0,0}, size3i_t{192,192,192}, check.data(), 0);
    reader->close();
    for (int ii=0; ii<192; ++ii) {
      for (int jj=0; jj<192; ++jj) {
        for (int kk=0; kk<192; ++kk) {
          float value_expect = expect(ii, jj, kk);
          float value_actual = check[(std::int64_t)ii*192*192 + (std::int64_t)jj*192 + kk];
          if (std::abs(value_actual - value_expect) > 0.001) {
            if (!TEST_EQUAL_FLOAT(value_actual, value_expect, 0.001)) {
              std::cout << "FAIL at (" << ii << "," << jj << "," << kk << ")\n";
              return;
            }
          }
        }
      }
    }
  }
}

#ifdef HAVE_SD

/**
 * Opening an empty file created by the old ZGY accessor and cloud plug-in
 * is not allowed.
 */
void
test_reopen_zgycloud()
{
  Test_Utils::CloudFileAutoDelete cad
    ("reopen_zgycloud.zgy", default_sd_context_rw());
  {
    // Expected to fail because the header segment is not padded.
    // Otherwise this file, containing no bricks, would have worked.
    const std::string oldname = get_sdtestdata("EmptyOldFile.zgy");
    Test_Utils::copy_sd_to_sd(oldname, cad.name());
    std::shared_ptr<IZgyWriter> writer;
    must_throw("Only files uploaded by OpenZGY", [&](){
      writer = IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(cad.name())
                                  .iocontext(default_sd_context_rw()));
    });
  }
  {
    // Expected to fail because the first segment has both header and bulk.
    const std::string oldname = get_sdtestdata("Synt2.zgy");
    Test_Utils::copy_sd_to_sd(oldname, cad.name());
    std::shared_ptr<IZgyWriter> writer;
    must_throw("Only files uploaded by OpenZGY", [&](){
      writer = IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(cad.name())
                                  .iocontext(default_sd_context_rw()));
    });
  }
  {
    // Expected to fail because only the current version (3/4) can be updated,
    // and that check should be done first because that error is more clear.
    // Without the check it would still fail the "uploaded by OpenZGY" test.
    const std::string oldname = get_sdtestdata("Synt2-v1.zgy");
    Test_Utils::copy_sd_to_sd(oldname, cad.name());
    std::shared_ptr<IZgyWriter> writer;
    must_throw("version is too old for this library", [&](){
      writer = IZgyWriter::reopen(ZgyWriterArgs()
                                  .filename(cad.name())
                                  .iocontext(default_sd_context_rw()));
    });
  }
}

#endif

static std::vector<std::int16_t>
make_sequence(std::int64_t start, std::int64_t count, std::int64_t step)
{
  std::vector<std::int16_t> result(count);
  std::int16_t value{static_cast<std::int16_t>(start % 32768)};
  for (std::int64_t ii = 0; ii < (std::int64_t)result.size(); ii += step) {
    // Assuming no rollover.
    result[ii] = value++;
  }
  return result;
}

template<typename T>
static std::string
array_to_string(const std::vector<T>& vec, double scale)
{
  if (vec.empty())
    return "()";
  std::stringstream ss;
  for (T it : vec)
    ss << ", " << (it / scale);
  return "(" + ss.str().substr(2) + ")";
}

#ifdef HAVE_SD

/**
 * Extract the segment information from FileStatistics -> xx_segments()
 * as the abbreviated 3-element vector of (first, all_middle, last)
 * and expand that to a list of all the segment sizes.
 * Yes this is way too roundabout. Also it obviously cannot check
 * that the segments that are expected to be the same size actually are.
 *
 * Note that if the FileStatistics come from a file opened for write
 * then the result depends on the segsplit setting. Data still in the
 * write buffer is assigned to the last segment; it won't be split
 * into the correctly sized segments until xx_close().
 */
static std::vector<std::int64_t>
list_segments(std::shared_ptr<const OpenZGY::FileStatistics> fs, int verbose)
{
  const std::vector<std::int64_t>& segs = fs->segmentSizes();
  const std::int64_t fsize = fs->fileSize();
  const std::int64_t numsegs =
    (fsize <= segs[0] && segs.size() == 1 ? 1 :
     fsize <= segs[0] + segs[1] && segs.size() == 2 ? 2 :
     segs.size() != 3 ? 0 :
     2 + (fsize - segs[0] - segs[2] + (segs[1] - 1)) / segs[1]);
  std::vector<std::int64_t> result;
  if (segs.size() < 3) {
    result = segs;
  }
  else {
    result.push_back(segs[0]);
    for (std::int64_t seg = 1; seg < numsegs - 1; ++seg)
      result.push_back(segs[1]);
    result.push_back(segs[2]);
  }
  if (verbose)
    std::cout << "segments " << array_to_string(result, 512*1024)
              << " bytes " << fsize
              << std::endl;
  return result;
}

static std::vector<std::int64_t>
list_segments(const std::string& name, int verbose)
{
#if 0
  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(name, default_sd_context_rw());
  std::shared_ptr<const OpenZGY::FileStatistics> fs = reader->filestats();
  reader->close();
  return list_segments(fs, verbose);
#else
  // Bypass OpenZGY completely. This will get the actual size of each
  // segment instead of assuming that all but the first and last have
  // the same size.
  std::vector<std::int64_t> result = Test_Utils::get_segsizes(name);
  if (verbose)
    std::cout << "segments " << array_to_string(result, 512*1024)
              << std::endl;
  return result;
#endif
}

static void
test_reopen_sd()
{
  typedef OpenZGY::IZgyWriter::size3i_t size3i_t;
  Test_Utils::CloudFileAutoDelete cad("reopen_sd.zgy", default_sd_context_rw());
  OpenZGY::SeismicStoreIOContext context =
    OpenZGY::SeismicStoreIOContext(*default_sd_context_rw())
    .segsplit(8)
    .segsize(2);
  if (verbose())
    std::cout << std::endl;

  std::shared_ptr<OpenZGY::IZgyWriter> writer =
    OpenZGY::IZgyWriter::open(ZgyWriterArgs()
                              .iocontext(&context)
                              .filename(cad.name())
                              .size(2*64, 3*64, 15*64) // 45 MB
                              .datatype(SampleDataType::int16)
                              .datarange(-32768,+32767));

  // Write 11 bricks, 5.5 MB (2.75 segs), plus 12 lowres => 5.75 data segs
  // Lowres: lod1 6 bricks, lod2 3 bricks, lod3 2 bricks, lod4 1 brick
  std::vector<std::int16_t> data1 = make_sequence(1, 11*64, 64);
  writer->write(size3i_t{0,0,0}, size3i_t{1,1,11*64}, data1.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate}, nullptr);

  // Check segment size before closing; note that the last segment
  // is the data buffer which can (if segsplit>1) cover more than
  // one real segment.
#ifndef TEST_WRITETHREADS // If writethreads forced on, this will fail.
  {
    std::vector<std::int64_t> check =  list_segments(writer->filestats(), verbose());
    if (TEST_EQUAL(check.size(), (std::size_t)2)) {
      TEST_EQUAL(check[0], 512*1024);
      TEST_EQUAL(check[1/**/], 23*512*1024);
    }
  }
#endif

  writer->close();

  // This re-opens the file for read to check the segment sizes.
  {
    std::vector<std::int64_t> check =  list_segments(cad.name(), verbose());
    if (TEST_EQUAL(check.size(), (std::size_t)7)) {
      TEST_EQUAL(check[0], 512*1024);
      TEST_EQUAL(check[1], 4*512*1024);
      TEST_EQUAL(check[2], 4*512*1024);
      TEST_EQUAL(check[3], 4*512*1024);
      TEST_EQUAL(check[4], 4*512*1024);
      TEST_EQUAL(check[5], 4*512*1024);
      TEST_EQUAL(check[6], 3*512*1024);
    }
  }

  // Write 11 bricks in addition to the 23 already present. Total now 34.
  // One header segment and 9 data segments.
  // Also the same 12 lowres bricks need to be updated. The three last
  // ones are found in the re-opened last segment and can be updated there.
  // Provided that segsplit is large enough to hold everything we want
  // to update this time around. 8 is sufficient, i.e. 32 bricks.
  // The remaining 9 cause more space to be allocated so we end up with
  // 43 bricks or 10.75 data bricks total.
  writer =
    OpenZGY::IZgyWriter::reopen(ZgyWriterArgs()
                                .iocontext(&context)
                                .filename(cad.name()));
  std::vector<std::int16_t> data2 = make_sequence(11*64+1, 11*64, 64);
  writer->write(size3i_t{64,0,0}, size3i_t{1,1,11*64}, data2.data());
  writer->finalize(std::vector<OpenZGY::DecimationType>{OpenZGY::DecimationType::Decimate}, nullptr);

  // Check segment size before closing; see notes above.
  // All but the last segments already written will remain.
#ifndef TEST_WRITETHREADS // If writethreads forced on, this will fail.
  {
    std::vector<std::int64_t> check =  list_segments(writer->filestats(), verbose());
    if (TEST_EQUAL(check.size(), (std::size_t)7)) {
      // Previous write.
      TEST_EQUAL(check[0], 512*1024);
      TEST_EQUAL(check[1], 4*512*1024);
      TEST_EQUAL(check[2], 4*512*1024);
      TEST_EQUAL(check[3], 4*512*1024);
      TEST_EQUAL(check[4], 4*512*1024);
      TEST_EQUAL(check[5], 4*512*1024);
      // Reopened 3 bricks, this and the new 15 bricks still not flushed.
      TEST_EQUAL(check[6], 23*512*1024);
    }
  }
#endif

  writer->close();

  std::shared_ptr<OpenZGY::IZgyReader> reader =
    OpenZGY::IZgyReader::open(cad.name(), &context);

#ifndef TEST_WRITETHREADS // If writethreads forced on, this will fail.
  {
    std::vector<std::int64_t> check =
      list_segments(reader->filestats(), verbose());
    if (TEST_EQUAL(check.size(), (std::size_t)12)) {
      TEST_EQUAL(check[0], 512*1024);
      TEST_EQUAL(check[1], 4*512*1024);
      TEST_EQUAL(check[2], 4*512*1024);
      TEST_EQUAL(check[3], 4*512*1024);
      TEST_EQUAL(check[4], 4*512*1024);
      TEST_EQUAL(check[5], 4*512*1024);
      TEST_EQUAL(check[6], 4*512*1024);
      TEST_EQUAL(check[7], 4*512*1024);
      TEST_EQUAL(check[8], 4*512*1024);
      TEST_EQUAL(check[9], 4*512*1024);
      TEST_EQUAL(check[10], 4*512*1024);
      TEST_EQUAL(check[11], 3*512*1024);
    }
  }

  {
    // Same as above but more thorough as it bypasses OpenZGY.
    std::vector<std::int64_t> check =
      list_segments(cad.name(), verbose());
    if (TEST_EQUAL(check.size(), (std::size_t)12)) {
      TEST_EQUAL(check[0], 512*1024);
      TEST_EQUAL(check[1], 4*512*1024);
      TEST_EQUAL(check[2], 4*512*1024);
      TEST_EQUAL(check[3], 4*512*1024);
      TEST_EQUAL(check[4], 4*512*1024);
      TEST_EQUAL(check[5], 4*512*1024);
      TEST_EQUAL(check[6], 4*512*1024);
      TEST_EQUAL(check[7], 4*512*1024);
      TEST_EQUAL(check[8], 4*512*1024);
      TEST_EQUAL(check[9], 4*512*1024);
      TEST_EQUAL(check[10], 4*512*1024);
      TEST_EQUAL(check[11], 3*512*1024);
    }
  }
#endif

  std::vector<std::int16_t> check1(11*64);
  std::vector<std::int16_t> check2(11*64);
  reader->read(size3i_t{0,0,0}, size3i_t{1,1,11*64}, check1.data(), 0);
  reader->read(size3i_t{64,0,0}, size3i_t{1,1,11*64}, check2.data(), 0);

  for (std::size_t ii = 0; ii < check1.size(); ++ii)
    if (check1[ii] != data1[ii])
      if (!TEST_EQUAL(check1[ii], data1[ii]))
        break;

  for (std::size_t ii = 0; ii < check2.size(); ++ii)
    if (check2[ii] != data2[ii])
      if (!TEST_EQUAL(check2[ii], data2[ii]))
        break;
}

#endif // HAVE_SD

class Register_reopen
{
public:
  Register_reopen()
  {
    register_test("reopen.plain",            test_reopen_plain);
    register_test("reopen.not_if_meta",      test_reopen_not_if_meta);
    register_test("reopen.not_if_final",     test_reopen_not_if_final);
    register_test("reopen.not_if_compress",  test_reopen_not_if_compress);
    register_test("reopen.alignment",        test_reopen_alignment);
    register_test("reopen.create_tiny",      test_reopen_create_tiny);
    register_test("reopen.update_tiny",      test_reopen_update_tiny);
    register_test("reopen.update_too_old",   test_reopen_update_too_old);
    register_test("reopen.update_bad_vt",    test_reopen_update_bad_vt);
    register_test("reopen.empty",            test_reopen_empty);
    register_test("reopen.bad_histogram",    test_reopen_bad_histogram);
    register_test("reopen.reopen",           test_reopen);
    register_test("reopen.incr",             test_reopen_incr);
    register_test("reopen.keepopen",         test_reopen_keepopen);
    register_test("reopen.rwlod",            test_reopen_rwlod);
    register_test("reopen.setversion",       test_reopen_setversion);
    register_test("reopen.zgypublic",        test_reopen_zgypublic);
    register_test("reopen.track_changes",    test_reopen_track_changes);
#ifdef HAVE_SD
    register_sd_test("reopen.plain_sd",         test_reopen_plain_sd);
    register_sd_test("reopen.sd",               test_reopen_sd);
    register_sd_test("reopen.zgycloud",         test_reopen_zgycloud);
#endif
  }
} dummy_reopen;

} // namespace
