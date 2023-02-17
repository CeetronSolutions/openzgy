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

//Adapted from: Zgy/PlatformLib/RoundAndClip.h and others.
#pragma once

#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>
#include <memory.h>

/**
 * \file roundandclip.h
 * \brief Conversion between scalar types.
 *
 * The original code in the ZGY library and in Salmon went a bit overboard
 * with optimizations. Up to and including inline assembly code. It is
 * possible to bring this back, but not until it can be verified that
 * those hacks are still relevant with current compilers.
 *
 * Note: Use the more portable std::isnan() in \<cmath\> for all platforms.
 * The bare "isnan" might not exist on all platforms, and/or it might be
 * a macro which gets undefined if \<cmath\> happens to be included elsewhere.
 * For similar reasons use std::isfinite unless it turns out to be proveably
 * slower than the alternatives.
 *
 * Update: It has been proven. In particular with msvc, the methods in the
 * standard library are *s*l*o*w. Hand coding some functions using bit-fiddling
 * is ugly but useful. I have run a single test case where the standard
 * functions had 1.4 times slowdown in g++. And more that 3 times slowdown
 * in msvc / visual studio.
 *
 * Risks:
 * - Algorithm correctness. (Verify with unit tests).
 *
 * - Correctness for odd types, e.g. const volatile long double. (Add tests).
 *
 * - SNaN handling or other subtle rules might differ. (We shouldn't care).
 *
 * - This kind of micro-optimizing is problematic because the effect is very
 *   compiler dependent. What is faster today might end up slower in a newer
 *   compiler. Even reorganizing the calling function (which is obviously
 *   calling these tests in a tight loop) might help the optimizer do a better
 *   job.
 *
 * - The detailed performance measurements were done just for isfinite(float),
 *   and only in one specific tight loop where the test was inlined.
 *   Applying the fix also to double and also to isnan might not be a good idea.
 *
 * - Would isnan() work better as just (x!=x) ? Would that be safe?
 *
 * - Looking at the library as a whole, the time spend in this particular
 *   function might not be that noticeable. So the change might not have been
 *   worth the effort. Especially with some of the callers already switched
 *   to sse2 code.
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * IsFiniteT() and IsNanT() use template magic to avoid having to test integral
 * types. They also avoid std::isfinite() and std::isnan() because those are
 * inefficient in g++ (use sse2 where regular instructions would be faster)
 * and criminally inefficient in msvc (use a function call, even in the non
 * standard _isfinite() and friends).
 *
 * iec559/ieee754 float can be interpreted as integral sign/value. Mask out the
 * sign bit. If what remains is exactly 0x7F800000 (for floats) then the number
 * is infinite. Anything larger is a NaN. No need to distinguish between
 * positive and negative infinity, nor between quiet and signalling NaN.
 * So the sign bit is not needed. The difference between Inf and NaN checks
 * is just a less vs. less-or-equal test.
 *
 * Implementation note about the cast inside the non-specialized templates.
 * With g++ this is not needed because the compiler sees that isfinite
 * and isnan will never be called with an integer, so it doesn't matter
 * that there are no overloads for integral types. On windows, the compiler
 * reports an error. The cast to long double removes the error.
 * To avoid an unnecessary cast, explicit specialization for float and
 * double should be made even if they just forward to the std:: version.
 *
 * Caveat: Even though isnan *can* be implemented very similar to isfinite,
 * it might not need to be. Testing (x!=x) or std::isnan() might be better.
 */
template <typename T> inline bool IsFiniteT(T value)
{
  if (!std::numeric_limits<T>::is_specialized ||
      std::numeric_limits<T>::is_integer) {
    return true;
  }
  /*
  else if (std::is_same<T, float>::value || std::is_same<T, const float>::value) {
    return std::isfinite(static_cast<float>(value));
  }
  else if (std::is_same<T, double>::value || std::is_same<T, const double>::value) {
    return std::isfinite(static_cast<double>(value));
  }*/
  else {
    return std::isfinite(static_cast<long double>(value));
  }
}

template <> inline bool IsFiniteT<float>(float value)
{
#if 1 // Measured to be faster.
  if (std::numeric_limits<float>::is_iec559) {
    std::uint32_t punnedVal;
    memcpy(&punnedVal, &value, sizeof(std::uint32_t));
    return ((punnedVal & 0x7FFFFFFF) < 0x7F800000);
  }
  else {
    return std::isfinite(value);
  }
#else
  return std::isfinite(value);
#endif
}

template <> inline bool IsFiniteT<double>(double value)
{
#if 0 // Not sure yet whether this is faster.
  if (std::numeric_limits<double>::is_iec559) {
    std::uint64_t punnedVal;
    memcpy(&punnedVal, &value, sizeof(std::uint64_t));
    return ((punnedVal & 0x7FFFFFFFFFFFFFFF) < 0x7FF0000000000000);
  }
  else {
    return std::isfinite(value);
  }
#else
  return std::isfinite(value);
#endif
}

template <typename T> inline bool IsNanT(T value)
{
  if (!std::numeric_limits<T>::is_specialized ||
      std::numeric_limits<T>::is_integer) {
    return false;
  }
  else {
    return std::isnan(static_cast<long double>(value));
  }
}

template <> inline bool IsNanT<float>(float value)
{
#if 0 // Only slightly faster on Windows, maybe.
  if (std::numeric_limits<float>::is_iec559) {
    std::uint32_t punnedVal;
    memcpy(&punnedVal, &value, sizeof(std::uint32_t));
    return ((punnedVal & 0x7FFFFFFF) > 0x7F800000);
  }
  else {
    return std::isnan(value);
  }
#else
  return std::isnan(value);
#endif
}

template <> inline bool IsNanT<double>(double value)
{
#if 0 // Not sure yet whether this is faster.
  if (std::numeric_limits<double>::is_iec559) {
    std::uint64_t punnedVal;
    memcpy(&punnedVal, &value, sizeof(std::uint64_t));
    return ((punnedVal & 0x7FFFFFFFFFFFFFFF) > 0x7FF0000000000000);
  }
  else {
    return std::isnan(value);
  }
#else
  return std::isnan(value);
#endif
}

/**
 * Simple round-to-nearest from a floating point number to 32-bit int.
 * Result is unspecified if input is not finite or does not fit in an int.
 */
static inline std::int32_t RoundD2I(double d)
{
  return d>0 ? (int)(d+0.5) : (int)(d-0.5);
}

/**
 * Cast a double to another type.
 * If the target type is floating point, this is a simple cast.
 * If the target type is integral the result is rounded to nearest
 * and clipped to the target's value range. If you supply a second
 * argument, NaN in the input is replaced with this value. plus or
 * minus infinite are still clipped to the target range.
 * If you do not supply a second argument, NaN input will result
 * in an undefined return value. It probably ends up being
 * 0x80000000 for signed/unsigned int and 0 for narrowed integral
 * types, but the caller should not depend on this.
 */
template <typename T>
T RoundAndClip(double value)
{
  // Can skip this (and avoid an annoying warning on windows)
  // if we are sure that all floating types are explicitly
  // specialized below. Or replace it with an assert.
  // But an assert could be even worse due to the performance hit.
  if (!std::numeric_limits<T>::is_integer)
    return static_cast<T>(value);

  // Be careful changing this code, it is more subtle than it looks.
  // If T is a 64-bit integer, this has greater precision than a double.
  // E.g. for T == unsigned long long and an input value of, say,
  // 2^64+1 this is clearly outside the valid range (max 2^64-1).
  // But both those values end up being rounded to 2^64 when assigned
  // to a double. So if the test used the more intuitive > instead of
  // >= it would fail. And the assignment at the end would overflow,
  // resulting in a hardware dependant result. 0x8000000000000000 on
  // Pentium. NaNs are supposed to fail both tests; in this case I want
  // the 0x8xx result as a fallback in case the caller didn't remove NaNs.

  if (value <= (double)std::numeric_limits<T>::min())
    return std::numeric_limits<T>::min();
  if (value >= (double)std::numeric_limits<T>::max())
    return std::numeric_limits<T>::max();

  // Since C++ by default rounds toward zero, we need special handling.
  // In the general case we cannot use RoundD2I(), and round() is slow.
  if (value < 0)
    return static_cast<T>(value - 0.5);
  return static_cast<T>(value + 0.5);
}

/**
 * Cast a double to another type, including substituting NaN
 * with a value supplied by the caller if and only if the target is integral.
 */
template <typename T>
T RoundAndClip(double value, T nan)
{
  if (!std::numeric_limits<T>::is_integer)
    return static_cast<T>(value);
  else if (IsNanT(value))
    return nan;
  else
    return RoundAndClip<T>(value);
}

} // namespace
