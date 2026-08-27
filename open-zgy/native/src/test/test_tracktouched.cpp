// Copyright 2017-2023, SLB
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
#include "../impl/tracktouched.h"
#include "../impl/logger.h"
#include "../impl/enum.h"

#include <tuple>
#include <cstdint>
#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace InternalZGY;

/**
 * \file test_tracktouched.cpp
 *
 * Test plan for class TrackTouched
 *
 * When on the fly generation is active, there are a few invariants:
 *
 * - Any write of a single brick-column (or single brick) should never
 *   trigger generating more than one low resolution brick-column
 *   (or single brick) at each level.
 *
 * - If a write causes a brick to be generated for level N, the same
 *   write will have caused at least one brick in each of levels
 *   1..N-1 to have been generated. Otherwise there would not have
 *   been any reason for that brick to suddenly become eligible.
 *
 * Consider this survey seen from above.
 * The brick-columns with no number inside
 * are not part of the survey; they are just
 * shown to better visualize what is going on.
 *
 *      <--- inline numbers -->
 *      0      32      64      96      128
 *  0   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |   |   |   |
 * 16   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      |10 |11*|12 |13*|14 |15*|16 |17*|18*|   |   |   |
 * 32   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      |19 |20 |21 |22 |23 |24 |25 |26 |27 |   |   |   |
 * 48   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      |28 |29*|30 |31#|32 |33*|34 |35#|36#|   |   |   |
 * 64   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      |37 |38*|39 |40#|41 |42*|43 |44#|45#|   |   |   |
 * 80   +---+---+---+---+---+---+---+---+---+---+---+---+
 *      |   |   |   |   |   |   |   |   |   |   |   |   |
 *
 *
 * surveysize in bricks: 9, 5, 3
 * 5 lod levels:
 *    lod1  5*3*2
 *    lod2  3*2*1
 *    lod3  2*1*1
 *    lod4  1*1*1
 *
 * Test 1: Write one brick-column at a time, in the order shown.
 * (*) means that lod1 can now be generated. (#) means lod1, lod2,
 * and possibly higher can be generated.
 *
 * Some examples:
 *  - 11: The 4 brick columns needed for the first lod1 column has been written.
 *  - 18: Only 2 brick columns needed here for the lod1 column.
 *  - 31: The 16 colums needed for the 4 columns needed for the first lod2.
 *  - 35: One more lod1, and second lod2 is ready.
 *  - 36: Only 4 and not 16 input columns neded for a lod2 here.
 *  - 44: Ditto
 *  - 45: The entire survey is now ready. Generate lod 3 and lod4.
 *
 * Alternative test: The very first brick (not the whole brick column)
 * is the last brick written, not the first.
 *
 * - 11: Will not trigger building a lod1.
 * - 31: Will now trigger a lod1 but not a lod2
 *
 * Test 2: Make the vertical size 51 bricks, so we get 7 lod levels.
 *   lod1 26,
 *   lod2 13,
 *   lod3  7,
 *   lod4  4,
 *   lod5  2,
 *   lod6  1
 *
 * The order of lod generation should be the same, except for step 45
 * which will now have a few extra levels to generate.
 *
 * Test 3: The very first brick-columnn will be flagged partially
 * written, which will cause several low resolutiom bricks to be delayed.
 *
 * Test 4: Write three brick-columns at a time, using the same order.
 *
 * Other tests are documented only above the test function itself.
 */

namespace {
#if 0
}
#endif

/**
 * Add a constructor to InteralLodMode for testing convenience,
 * taking arguments (level, force) IN THAT ORDER.
 * The main point of InteralLodMode is to avoid positional
 * arguments with a non-obvious ordering. In the unit test
 * it should be safe enough.
 */
struct MyLevelAndForce : public InternalZGY::LevelAndForce
{
  MyLevelAndForce(int the_level, int the_force)
  {
    this->level = the_level;
    this->force = the_force;
  }
};

// Moved from tracktouched, as the overloads are now only used in tests.

static std::int64_t
countIO(const TrackTouched::tasklist_t& tasklist)
{
  std::int64_t result{0};
  for(const TrackTouched::task_t& task : tasklist)
    result += TrackTouched::countIO(task);
  return result;
}

static std::int64_t
countIO(
     const TrackTouched& touched,
     const InternalZGY::LevelAndForce& lodmode,
     const TrackTouched::indexlist_t& roi)
{
  return countIO(touched.getWork(lodmode, roi));
}

/**
 * Duplicated from impl/meta.cpp
 */
static std::vector<std::array<std::int64_t,3>>
calcLodSizes(const std::array<std::int64_t,3>& size_in_bricks)
{
  std::array<std::int64_t,3> size = size_in_bricks;
  std::vector<std::array<std::int64_t,3>> result;
  if (size[0] < 1 || size[1] < 1 || size[2] < 1) {
    result.push_back(std::array<std::int64_t,3>{0,0,0});
  }
  else {
    result.push_back(size);
    while (size[0] > 1 || size[1] > 1 || size[2] > 1) {
      for (int ii=0; ii<3; ++ii)
        size[ii] = (size[ii] + 1) / 2;
      result.push_back(size);
    }
  }
  return result;
}

/**
 * Duplicated from impl/meta.cpp
 */
static std::vector<std::int64_t>
calcLutOffsets(const std::vector<std::array<std::int64_t,3>>& lods_in,
               bool isalpha)
{
  std::vector<std::array<std::int64_t,3>> lods(lods_in);
  std::vector<std::int64_t> result;
  std::reverse(lods.begin(), lods.end());
  std::int64_t pos = 0;
  for (const auto& e : lods) {
    result.push_back(pos);
    pos += (e[0] * e[1] * (isalpha ? 1 : e[2]));
  }
  std::reverse(result.begin(), result.end());
  result.push_back(pos);
  return result;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Create a TrackedTouched instance with higher level arguments:
 * size instead of lodsizes and brickoffsets.
 * Also supply a custom logger, because this test is so low level
 * that it cannot use e.g. the default one set up in ZgyUnternalBulk.
 * Note that the standard logger writes to stderr by default.
 */
static std::shared_ptr<TrackTouched>
makeTracker(
     const TrackTouched::index3_t& size_in_bricks)
{
  const TrackTouched::lodsizes_t lodsizes = calcLodSizes(size_in_bricks);
  const TrackTouched::offsets_t brickoffsets = calcLutOffsets(lodsizes, false);
  auto logger = LoggerBase::standardCallback
    (verbose() ? 9 : LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"),
     "testlogger: ", "");
  return std::make_shared<TrackTouched>
    (lodsizes, brickoffsets, logger);
}

static std::string
formatWork(const TrackTouched::tasklist_t& work, int sequence)
{
  std::stringstream ss;
  ss << std::setw(4) << sequence << ":";
  for (const auto& it : work)
    ss << " [" << it.lod
       << "](" << it.pos[0]
       << ", " << it.pos[1]
       << ", " << it.pos[2]
       << ")";
  return ss.str();
}

static void
checkWork(const TrackTouched::tasklist_t& work, int sequence, const std::vector<int>& expect)
{
  if (verbose() /* && !work.empty()*/)
    std::cerr << formatWork(work, sequence) << std::endl;
  if (TEST_CHECK(sequence >= 1 && sequence <= (int)expect.size())) {
    if (TEST_EQUAL((int)work.size(), expect[sequence - 1])) {
      TEST_CHECK((int)work.size() <= 0 || work[0].lod == 1);
      TEST_CHECK((int)work.size() <= 1 || work[1].lod == 2);
      TEST_CHECK((int)work.size() <= 2 || work[2].lod == 3);
      TEST_CHECK((int)work.size() <= 3 || work[3].lod == 4);
    }
  }
}

/**
 * Remove duplicates from a vector of ijk positions.
 */
static TrackTouched::tasklist_t
removeDuplicates(const TrackTouched::tasklist_t& list)
{
  if (list.size() <= 1)
    return list;
  static auto lessthan = [](const TrackTouched::task_t& a, const TrackTouched::task_t& b)
    {
      if (a.lod != b.lod) return a.lod < b.lod;
      return a.pos < b.pos;
    };
  static auto equals = [](const TrackTouched::task_t& a, const TrackTouched::task_t& b)
  {
    return a.lod == b.lod && a.pos == b.pos;
  };
  TrackTouched::tasklist_t copy(list);
  std::sort(copy.begin(), copy.end(), lessthan);
  copy.erase(std::unique(copy.begin(), copy.end(), equals), copy.end());
  return copy;
}

/*
 * Run thru an example where the expected results have been calculated by hand.
 * Mostly without peeking at the output from the test run.
 *
 * The "tall" parameter makes a very slight change: The vertical size becomes
 * much larger than the horizontal, so it is the vertical extent that
 * determines the number of lod layers. Since the compute currently only deals
 * with brick columns, the lowres count is the only observable change.
 * And that shows only in the very last write.
 */
static void
do_test_1(bool tall)
{
  if (verbose())
    std::cerr << std::endl;
  const TrackTouched::index3_t size_in_bricks{9, 5, tall ? 35 : 3};
  std::shared_ptr<TrackTouched> tt =  makeTracker(size_in_bricks);

  // Refer to the figure. "*" means one brick at level 1
  // should be ready. "#" expects that, and more.
  std::vector<int> expect_levels{
    0,0,0,0,0,0,0,0,0,
    0,1,0,1,0,1,0,1,1,
    0,0,0,0,0,0,0,0,0,
    0,1,0,2,0,1,0,2,2,
    0,1,0,2,0,1,0,3, (tall? 6 : 4)
  };

  std::array<std::int64_t,3>brickpos{0,0,0};
  int sequence{0};
  // Iterate with inlines faster than crosslines to match the figure.
  for (brickpos[1] = 0; brickpos[1] < size_in_bricks[1]; ++brickpos[1]) {
    for (brickpos[0] = 0; brickpos[0] < size_in_bricks[0]; ++brickpos[0]) {
      ++sequence;
      tt->setCWritten(TrackTouched::indexlist_t{brickpos}, /*lod=*/0, /*rmw=*/false);
      TrackTouched::tasklist_t work =
        tt->getWorkAndClear(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
      checkWork(work, sequence, expect_levels);
    }
  }
  //tt->showDirty(true, 1);
  //tt->showDirty(false, 1);
}

static void
test_1()
{
  do_test_1(false);
}

static void
test_2()
{
  do_test_1(true);
}

/**
 * Slight variation of test 1: The very first brick-column is marked as
 * partly written, and won't be fully written until after the main loop.
 * Several "expected" values will be one less than they were in test 1.
 */
static void
test_3()
{
  if (verbose())
    std::cerr << std::endl;
  const TrackTouched::index3_t size_in_bricks{9, 5, 3};
  std::shared_ptr<TrackTouched> tt = makeTracker(size_in_bricks);

  // Refer to the figure. "*" means one brick at level 1
  // should be ready. "#" expects that, and more.
  std::vector<int> expect_levels{
    0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,1,0,1,1,
    0,0,0,0,0,0,0,0,0,
    0,1,0,1,0,1,0,2,2,
    0,1,0,2,0,1,0,2,3,4
  };

  std::array<std::int64_t, 3>brickpos{ 0,0,0 };
  int sequence{ 0 };
  // Iterate with inlines faster than crosslines to match the figure.
  for (brickpos[1] = 0; brickpos[1] < size_in_bricks[1]; ++brickpos[1]) {
    for (brickpos[0] = 0; brickpos[0] < size_in_bricks[0]; ++brickpos[0]) {
      ++sequence;
      const TrackTouched::indexlist_t this_roi{ brickpos };
      const TrackTouched::indexlist_t no_roi{};
      tt->setCWritten(TrackTouched::indexlist_t{brickpos}, /*lod=*/0, /*rmw=*/(sequence == 1));
      // Without roi hint
      TrackTouched::tasklist_t work =
        tt->getWork(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
      // With roi hint, and reset the bricks we are asked to process
      TrackTouched::tasklist_t work_roi =
        tt->getWorkAndClear(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{ brickpos });
      // After clearing, nothing left behind, even without the roi.
      TrackTouched::tasklist_t work_not =
        tt->getWork(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
      checkWork(work, sequence, expect_levels);
      checkWork(work_roi, sequence, expect_levels);
      TEST_EQUAL(work_not.size(), 0);
    }
  }
  ++sequence;
  // To finish up, either properly write the partial first column,
  // or use mode=2 telling the system this is the last chance.
  tt->setCWritten(TrackTouched::indexlist_t{{0,0,0}}, /*lod=*/0, /*rmw=*/false);
  TrackTouched::tasklist_t work =
    tt->getWorkAndClear(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
  checkWork(work, sequence, expect_levels);
  //tt->showDirty(true, 1);
  //tt->showDirty(false, 1);
}

static void
test_4()
{
  if (verbose())
    std::cerr << std::endl;
  const TrackTouched::index3_t size_in_bricks{9, 5, 3};
  std::shared_ptr<TrackTouched> tt = makeTracker(size_in_bricks);

  std::array<std::int64_t, 3>brickpos{ 0,0,0 };
  int sequence{ 0 };
  TrackTouched::tasklist_t allwork;
  // Iterate with inlines faster than crosslines to match the figure.
  for (brickpos[1] = 0; brickpos[1] < size_in_bricks[1]; ++brickpos[1]) {
    for (brickpos[0] = 0; brickpos[0] < size_in_bricks[0]; ++brickpos[0]) {
      ++sequence;
      tt->setCWritten(TrackTouched::indexlist_t{brickpos}, /*lod=*/0, /*rmw=*/false);
      if ((sequence % 3) == 0 || sequence == 45) {
        TrackTouched::tasklist_t work =
          tt->getWorkAndClear(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
        if (verbose())
          std::cerr << formatWork(work, sequence) << std::endl;
        allwork.insert(allwork.end(), work.begin(), work.end());
      }
    }
  }

  // Verify that each low resolution brick was generated exactly
  // once. The test may have a false negative if some weird
  // position outside the survey is included. I am not that paranoid.
  TEST_EQUAL(allwork.size(), (std::size_t)24);
  allwork = removeDuplicates(allwork);
  TEST_EQUAL(allwork.size(), (std::size_t)24);
  TrackTouched::tasklist_t firstandlast{ allwork.front(), allwork.back() };
  TEST_EQUAL(formatWork(firstandlast, 0), std::string("   0: [1](0, 0, 0) [4](0, 0, 0)"));

  //tt->showDirty(true, 1);
  //tt->showDirty(false, 1);
}

static void
test_force()
{
  if (verbose())
    std::cerr << std::endl;
  const TrackTouched::index3_t size_in_bricks{9, 5, 3};
  std::shared_ptr<TrackTouched> tt = makeTracker(size_in_bricks);

  std::vector<int> expect_levels{ 0,0,4,0,4 };

  tt->setCWritten(TrackTouched::indexlist_t{{0,0,0}}, /*lod=*/0, /*rmw=*/false);
  TrackTouched::tasklist_t work;

  // No-op with force 0
  work = tt->getWorkAndClear(MyLevelAndForce(-1, 0), TrackTouched::indexlist_t{});
  checkWork(work, 1, expect_levels);

  // No fully ready bricks yet, so nothing to do for incremental
  work = tt->getWorkAndClear(MyLevelAndForce(-1, 1), TrackTouched::indexlist_t{});
  checkWork(work, 2, expect_levels);

  // Generate also for partial bricks with force 2
  work = tt->getWorkAndClear(MyLevelAndForce(-1, 2), TrackTouched::indexlist_t{});
  checkWork(work, 3, expect_levels);

  // Already cleared, so there shouldn't be any work left....
  work = tt->getWorkAndClear(MyLevelAndForce(-1, 2), TrackTouched::indexlist_t{});
  checkWork(work, 4, expect_levels);

  // Unless a full rebuild is requested with force 3,
  // causing all 24 low resolution brick columns to be generated.
  work = tt->getWorkAndClear(MyLevelAndForce(-1, 3), TrackTouched::indexlist_t{});
  TEST_EQUAL(work.size(), (std::size_t)24);
}

static void
test_roi()
{
  if (verbose())
    std::cerr << std::endl;
  const TrackTouched::index3_t size_in_bricks{9, 5, 3};
  std::shared_ptr<TrackTouched> tt = makeTracker(size_in_bricks);

  std::array<std::int64_t, 3>brickpos{ 0,0,0 };
  for (brickpos[0] = 0; brickpos[0] < size_in_bricks[0]; ++brickpos[0])
    for (brickpos[1] = 0; brickpos[1] < size_in_bricks[1]; ++brickpos[1])
      for (brickpos[2] = 0; brickpos[2] < size_in_bricks[2]; ++brickpos[2])
        tt->set1Written(brickpos, /*lod=*/0, /*rmw=*/false);

  // TEST: A region fully outside the survey
  // should never contain any work.
  TrackTouched::indexlist_t roi_outside{
    std::array<std::int64_t, 3>{-100,-100,-100},
    std::array<std::int64_t, 3>{100,200,300}
  };
  TrackTouched::tasklist_t work_outside = tt->getWork(MyLevelAndForce(-1, 1), roi_outside);
  TEST_EQUAL(work_outside.size(), (std::size_t)0);

  // TEST: ROI consisting of two bricks.
  // The roi hint normally means that nothing outside the region has changed,
  // so it is pointless to check. In this test I am lying. Passing a too
  // small region. Limited to the first and last brick. Not all the lod1 will
  // be generated. Which in turn means that few (in this case, one) lod2 will
  // be generated. The reason why we get even one, is that the very last lod2
  // brick column only requires a single lod1 column input. The other three
  // are outside the survey.
  TrackTouched::indexlist_t roi{
    std::array<std::int64_t, 3>{8,4,2},
    std::array<std::int64_t, 3>{0,0,0}
  };
  TrackTouched::tasklist_t work = tt->getWorkAndClear(MyLevelAndForce(-1, 1), roi);

  std::string actual(formatWork(work, 0));
  std::string expect("   0: [1](0, 0, 0) [1](4, 2, 0) [2](2, 1, 0)");
  if (verbose())
    std::cerr << actual << std::endl;
  TEST_EQUAL(actual, expect);

  TrackTouched::tasklist_t work_no_more = tt->getWork(MyLevelAndForce(-1, 1), roi);
  // We already did this, right?
  TEST_EQUAL(work_outside.size(), (std::size_t)0);

  //tt->showDirty(true, 1);
  //tt->showDirty(false, 1);
}

/**
 * A survey with just a single brick and thus no lod levels can very
 * easily trigger corner cases.
 */
static void
test_micro()
{
  typedef TrackTouched::index3_t index3_t;
  typedef TrackTouched::indexlist_t indexlist_t;

  if (verbose())
    std::cerr << std::endl;
  std::shared_ptr<TrackTouched> tt = makeTracker(index3_t{1,1,1});

  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 0), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 1), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 2), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 3), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(10, 2), indexlist_t{}).size(), 0);

  //tt->showDirty(true, 0);
  //tt->showDirty(false, 0);

  tt->setWritten(indexlist_t{{0,0,0}}, /*lod=*/0, /*rmw=*/false);

  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 0), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 1), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 2), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 3), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(10, 2), indexlist_t{}).size(), 0);

  TEST_EQUAL(countIO(*tt, MyLevelAndForce(-1, 3), indexlist_t{}), 0);

  //tt->showDirty(true, 0);
  //tt->showDirty(false, 0);
}

/**
 * A survey with just four bricks and one lod.
 */
static void
test_small()
{
  typedef TrackTouched::index3_t index3_t;
  typedef TrackTouched::indexlist_t indexlist_t;

  if (verbose())
    std::cerr << std::endl;
  std::shared_ptr<TrackTouched> tt = makeTracker(index3_t{2,2,1});

  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 0), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 1), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 2), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 3), indexlist_t{}).size(), 1);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(10, 2), indexlist_t{}).size(), 0);

  //tt->showDirty(true, 0);
  //tt->showDirty(false, 0);

  tt->setWritten(indexlist_t{{0,0,0}}, /*lod=*/0, /*rmw=*/false);

  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 0), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 1), indexlist_t{}).size(), 0);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 2), indexlist_t{}).size(), 1);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(-1, 3), indexlist_t{}).size(), 1);
  TEST_EQUAL(tt->getWork(MyLevelAndForce(10, 2), indexlist_t{}).size(), 1);

  // This counts also the 4 bricks completely outside the survey.
  TEST_EQUAL(countIO(*tt, MyLevelAndForce(-1, 3), indexlist_t{}), 9);

  //tt->showDirty(true, 0);
  //tt->showDirty(false, 0);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("touched.test1",           test_1);
      register_test("touched.test2",           test_2);
      register_test("touched.test3",           test_3);
      register_test("touched.test4",           test_4);
      register_test("touched.force",           test_force);
      register_test("touched.roi",             test_roi);
      register_test("touched.micro",           test_micro);
      register_test("touched.small",           test_small);
    }
  } dummy;
} // namespace for registration
