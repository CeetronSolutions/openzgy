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

#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>

// Define this if you want to look at the generated assembly code.
//#define LIBRARY_ONLY

#ifndef LIBRARY_ONLY
#include "test_all.h"
#include "test_utils.h"
#include "../impl/timer.h"
#include <chrono>
#include <random>
#include <iostream>
#include <functional>

using InternalZGY::PrintingTimer;

#endif

// Choose to run tests with float or double.
// Yes I could template this, but really I only need the float results.
// And templates would make looking at the assembler output more difficult.

#define FLOAT float

/**
 * \file range.cpp
 *
 * Experiment with different versions of the range() and
 * clip() functions. There are deliberately simplified.
 * E.g. there is no templating. Function names are declared
 * with "C" linkage to make them easier to spot in the
 * generated assembly code.
 *
 * What OpenZGY needs is a range() function that ignores
 * both NaN and Inf input, i.e. the functionality in range3()
 * below. This means that if range1() or range2() is used
 * instead then the result needs to be checked for bad values
 * and if found then fall back to using range3() instead.
 * The fallback will rarely be needed which means that if
 * range1() or range2() are notably faster then it might
 * make sense to try those first.
 *
 * What OpenZGY needs is a clip() function that replaces NaN
 * with a user supplied value and clips other input values,
 * including +/- Inf, to the value range of the target type.
 * As clip1() does. When the function is invoked then the result
 * of range() will be known. This means that clip() might know
 * a priori if no actual clipping or NaN handling will be needed.
 * If so, one of the simpler version of clip() might be used.
 *
 * All of the functions should allow parallelization with
 * OpenMP. There would be some cutoff on input size where the
 * overhead of OpenMP makes this pointless. If the cutoff ends
 * up smaller than 256 KB (probably the lowest practical
 * OpenZGY brick size) then parallelization would make sense
 * in most cases.
 *
 * NOTE: The tests in this file are somewhat pointless unless
 * the local timers are enabled, because the time measured by the
 * test framework will include the setup time which is much more
 * than the actual test. I can fix this by re-using a static
 * test data buffer but then the tests would depend on each other.
 */

/**
 * Compute the min/max range of all members in a buffer.
 * If the input contains NaN then the resulting range
 * will be (NaN,NaN). Else, if the input contains
 * +/-Inf then those will be included in the result.
 * If the input is empty the result will have min>max.
 *
 * With g++ -O3, double: 14 instructions, one of which is an
 * almost-always-true conditional jump.
 */
extern "C"
std::pair<FLOAT, FLOAT>
range1(const FLOAT *ptr, const FLOAT* end)
{
  FLOAT min = std::numeric_limits<FLOAT>::infinity();
  FLOAT max = -std::numeric_limits<FLOAT>::infinity();
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (!(min <= value))
      min = value;
    if (!(max >= value))
      max = value;
  }
  return std::make_pair(min, max);
}

/**
 * Compute the min/max range of all members in a buffer.
 * NaN in the input is ignored. If the input contains
 * +/-Inf then those will be included in the result.
 * If the input has no non-NaN values the result will
 * have min>max.
 *
 * With g++ -O3, double: 9 instructions, one of which
 * is an almost-always-true conditional jump.
 */
extern "C"
std::pair<FLOAT, FLOAT>
range2(const FLOAT *ptr, const FLOAT* end)
{
  FLOAT min = std::numeric_limits<FLOAT>::infinity();
  FLOAT max = -std::numeric_limits<FLOAT>::infinity();
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (min > value)
      min = value;
    if (max < value)
      max = value;
  }
  return std::make_pair(min, max);
}

/**
 * Compute the min/max range of all members in a buffer.
 * NaN and +/- inf in the input is ignored.
 * If the input has no finite values the result will
 * have min>max.
 *
 * With g++ -O3, double: 14 instructions, two of which
 * are almost-always-true conditional jumps. The feature
 * of returning the allfinite variable apparently
 * only costs one instruction.
 */
extern "C"
std::pair<FLOAT, FLOAT>
range3(const FLOAT *ptr, const FLOAT* end, bool *allfinite)
{
  FLOAT min = std::numeric_limits<FLOAT>::infinity();
  FLOAT max = -std::numeric_limits<FLOAT>::infinity();
  std::size_t bad = end - ptr;
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (std::isfinite(value)) {
      --bad;
      if (min > value)
        min = value;
      if (max < value)
        max = value;
    }
  }
  *allfinite = (bad == 0);
  return std::make_pair(min, max);
}

/**
 * Hand-optimized version of range2, i.e. NaN ignored
 * and Inf included. The inner loop is partly rolled out
 * so it handle 4 input values per iteration.
 *
 * With g++ -O3, double: 24 instructions, one of which is
 * an almost-always-true conditional jump. So there are
 * 6 instructions and 1/4 jump per input value. It is
 * also possible that the CPU might be slightly better
 * at prefetching in this case. I haven't benchmarked
 * this though.
 *
 * NOTE: Logic to handle size not a multiple of 4 is
 * missing; it complicates the code but won't contribute
 * noticeably to the execution time.
 *
 * Note: If I add code to handle both the first 0-3
 * values and the last 0-3 values in a way that makes
 * the pointers 256-bit aligned, and the compiler use
 * more efficient sse instructions?
 */
extern "C"
std::pair<FLOAT, FLOAT>
range4(const FLOAT *ptr, const FLOAT* end)
{
  FLOAT min = std::numeric_limits<FLOAT>::infinity();
  FLOAT max = -std::numeric_limits<FLOAT>::infinity();
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (min > value)
      min = value;
    if (max < value)
      max = value;
    value = *ptr++;
    if (min > value)
      min = value;
    if (max < value)
      max = value;
    value = *ptr++;
    if (min > value)
      min = value;
    if (max < value)
      max = value;
    value = *ptr++;
    if (min > value)
      min = value;
    if (max < value)
      max = value;
  }
  return std::make_pair(min, max);
}

/**
 * Here is my hand-optimized version of range3.
 * Ouch! What was I thinking?
 *
 * With g++ -O3 and few non-finite values, 3 jump instructions per
 * input value and a very real risk of branch prediction failure.
 *
 * As with the other variants it is possible to parallelize this with
 * OpenMP. Note that the first loop will almost always stop after the
 * first iteration. So only the second loop should run multi-threaded.
 */
extern "C"
std::pair<FLOAT, FLOAT>
range5(const FLOAT *ptr, const FLOAT* end)
{
  FLOAT min = std::numeric_limits<FLOAT>::infinity();
  FLOAT max = -std::numeric_limits<FLOAT>::infinity();
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (std::isfinite(value)) {
      min = max = value;
      break;
    }
  }
  // Due to the loop above I can skip the check for min > max
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (!std::isfinite(value)) continue;
    else if (min > value) min = value;
    else if (max < value) max = value;
  }
  return std::make_pair(min, max);
}

/**
 * Round and clip, correctly handling both NaN, Inf, and value outside
 * range. The assembly code is not trivial to follow but it looks like
 * it has about 4 jumps and the last one, for the +/- 0.5, will not be
 * predictable. If there are many Nan/Inf/Crop values then there may be
 * other branch prediction failures.
 *
 * Precisely -0.5 will map to -1, precisely +0.5 will map to +1.
 */
extern "C"
void
clip1(const FLOAT *ptr, const FLOAT *end, std::int16_t *out, FLOAT nan)
{
  while (ptr < end) {
    FLOAT value = *ptr++;
    if (std::isnan(value))
      *out++ = nan;
    else if (value <= (FLOAT)std::numeric_limits<std::int16_t>::min())
      *out++ = std::numeric_limits<std::int16_t>::min();
    else if (value >= (FLOAT)std::numeric_limits<std::int16_t>::max())
      *out++ = std::numeric_limits<std::int16_t>::max();
    else
      *out++ = static_cast<std::int16_t>(value < 0 ? value - 0.5 : value + 0.5);
  }
}

/**
 * Round and clip, only works correctly of the data is already known
 * to have no NaN/Inf/Clip values. A previous call to range(), which
 * we need anyway, might provide this information. Depending on which
 * range() variant was used.
 *
 * With g++ -O3 this looks like it gets two jumps, one of them not
 * predictable. If this is correct then using this version might be
 * pointless, as the non-predictable branch will probably be
 * responsible for most of the time spent. And clip1() also has just
 * one of those, I think.
 */
extern "C"
void
clip2(const FLOAT *ptr, const FLOAT *end, std::int16_t *out, FLOAT/*nan*/)
{
  while (ptr < end) {
    const FLOAT value = *ptr++;
    *out++ = static_cast<std::int16_t>(value < 0 ? value - 0.5 : value + 0.5);
  }
}

/**
 * As clip2() but try another way of doint the rounding. The assembly
 * code is not trivial to follow but it looks like there will still be
 * one non-predictable branch.
 */
extern "C"
void
clip3(const FLOAT *ptr, const FLOAT *end, std::int16_t *out, FLOAT/*nan*/)
{
  while (ptr < end) {
    const FLOAT value = *ptr++;
    *out++ = static_cast<std::int16_t>(std::floor(value + 0.5));
  }
}

/**
 * As clip2() but try yet another way of doing the rounding.
 * With g++ -O3 the call to std::round() is not intrinsic
 * and the code won't even ue sse since ::round() doesn't.
 * This will likely run way slower.
 */
extern "C"
void
clip4(const FLOAT *ptr, const FLOAT *end, std::int16_t *out, FLOAT/*nan*/)
{
  while (ptr < end) {
    const FLOAT value = *ptr++;
    *out++ = static_cast<std::int16_t>(std::round(value));
  }
}

#ifndef LIBRARY_ONLY

static void getdata(std::vector<float> *data, std::int64_t size)
{
  *data = Test_Utils::random_vector(size);
}

#if 0 // Enable if testing both float and double
static void getdata(std::vector<double> *data, std::int64_t size)
{
  *data = Test_Utils::random_double_vector(size);
}
#endif

typedef std::function<std::pair<FLOAT,FLOAT>(const FLOAT*, const FLOAT*)> range_fn;
typedef std::function<void(const FLOAT*, const FLOAT*, std::int16_t*, FLOAT)> clip_fn;

static void test_optimize_range(const char *name, const range_fn& fn)
{
  std::vector<FLOAT> data;
  getdata(&data, 10*1024*1024);
  PrintingTimer t1(name, 1, verbose());
  /*auto result = */fn(data.data(), data.data() + data.size());
  //std::cout << "Range is "
  //          << result.first << " to " << result.second << std::endl;
}

static void test_optimize_clip(const char *name, const clip_fn& fn)
{
  std::vector<FLOAT> data;
  getdata(&data, 10*1024*1024);
  std::vector<std::int16_t> result(data.size());
  PrintingTimer t1(name, 1, verbose());
  fn(data.data(), data.data() + data.size(), result.data(), -999.25);
}

static void test_optimize_range1()
{
  test_optimize_range("range1", range1);
}

static void test_optimize_range2()
{
  test_optimize_range("range2", range2);
}

static void test_optimize_range3()
{
  bool allfinite = false;
  test_optimize_range("range3",
                      [&allfinite](const FLOAT *beg, const FLOAT *end)
                      {
                        return range3(beg, end, &allfinite);
                      });
}

static void test_optimize_range4()
{
  test_optimize_range("range4", range4);
}

static void test_optimize_range5()
{
  test_optimize_range("range5", range5);
}

static void test_optimize_clip1()
{
  test_optimize_clip("clip1", clip1);
}

static void test_optimize_clip2()
{
  test_optimize_clip("clip2", clip2);
}

static void test_optimize_clip3()
{
  test_optimize_clip("clip3", clip3);
}

static void test_optimize_clip4()
{
  test_optimize_clip("clip4", clip4);
}

// Conclusions for range() with g++ -O3 on an Intel x86_64:
// range1 is more expensive than the others, so ignore that one.
// Unrolling the range2 loop, in range4, had too little effect.
// range5 wasn't as bad as I feared; possibly the more modern CPUs
// are better at handling unpredictable brances. But given that the
// results for range3 and range5 are similar I'll drop the latter.
// That leaves range2 (20 ms/10^7 values) and range3 (30 ms/10^7 values).
//
// There is a 30% speedup by trying range2 first and then using range3
// as a fallback in case an infinite range was returned. Given that
// the results are compiler- and hardware dependent an that the speed
// difference isn't that dramatic I will probably just stick with range3.
//
// Conclusions for range() with g++ -O3 on an Intel x86_64:
// Forget clip4; round() is slower. Not terribly so but it still
// performs worse. clip1 takes ~90 ms, clip2 ~70 ms, clip3 ~50 ms
// so it is possible to almost double the speed if it is reasonably
// easy to keep track of whether range() encountered any NaN or
// outside-range values.
//
// Note that in the real code the inner loop of clip() will also
// be doing a multiply-add which would cost the same in all cases.
// So the percentage increase in execution time goes down. If I find
// I need to parallelize the operation it is an open question whether
// using clip3 with clip1 only as a fallback is worth the effort.

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("optimize.range1",         test_optimize_range1);
      register_test("optimize.range2",         test_optimize_range2);
      register_test("optimize.range3",         test_optimize_range3);
      register_test("optimize.range4",         test_optimize_range4);
      register_test("optimize.range5",         test_optimize_range5);
      register_test("optimize.clip1",          test_optimize_clip1);
      register_test("optimize.clip2",          test_optimize_clip2);
      register_test("optimize.clip3",          test_optimize_clip3);
      register_test("optimize.clip4",          test_optimize_clip4);
    }
  } dummy;
}

#endif
