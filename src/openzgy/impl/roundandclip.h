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
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * IsFiniteT uses template magic to avoid having to test integral types.
 */
template <typename T> inline bool IsFiniteT(T /*value*/)
{
  return true;
}

template <> inline bool IsFiniteT<float>(float value)
{
  return std::isfinite(value);
}

template <> inline bool IsFiniteT<double>(double value)
{
  return std::isfinite(value);
}

template <> inline bool IsFiniteT<long double>(long double value)
{
  return std::isfinite(static_cast<double>(value));
}

/**
 * IsNanT uses template magic to avoid having to test integral types.
 */
template <typename T> inline bool IsNanT(T /*value*/)
{
  return false;
}

template <> inline bool IsNanT<float>(float value)
{
  return std::isnan(value);
}

template <> inline bool IsNanT<double>(double value)
{
  return std::isnan(value);
}

template <> inline bool IsNanT<long double>(long double value)
{
  return std::isnan(static_cast<double>(value));
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
