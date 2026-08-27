// Copyright 2017-2020, Schlumberger
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

//Adapted from: Zgy/UtilityLib/StatisticData
#pragma once

#include "../declspec.h"

#include <cstdint>
#include <math.h>
#include <cmath>
#include <string>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file statisticdata.h
 * \brief provides class \ref InternalZGY::StatisticData.
 */

/**
 * \brief Holds the result of computing statistics. Not thread
 * safe. Also has subtle issues with expandable histograms.
 *
 * \details
 *
 * Sorry for writing an entire essay in the class header.
 *
 * Thread safety:
 *
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 *
 * There are subtle issues mostly related to expandable histograms.
 *
 * A StatisticData with min>max is defined as empty, meaning that no
 * finite samples as been added.
 *
 * - An empty StatisticData must have zero count, sum, and sum of squares.
 * - An empty instance can still have a non-zero count of non-finite samples.
 * - A zero count doesn't necessarily mean that the instance is empty.
 * - Subtracting samples from an instance cannot narrow the min/max range.
 * - Subtracting samples cannot convert an instance from non-empty to empty.
 * - Sample counts may be negative, to allow removing samples.
 *
 * The reason why a zero count does not imply an empty instance is
 * that in the expression "stats+=(new-old)" will often not change the
 * count. So, the the temporary (new-old) will often have a zero
 * count. The sum, ssq, etc. in the temporary are still valid and
 * should be added to stats. And the min/max range range is
 * potentially widened.
 *
 * For consistency, executing "stats-=stats" will not result in an
 * empty instance. The range will be kept.
 *
 * Executing "stats+=some_empty_instance" is not a no-op because the
 * "empty" instance can represent a buffer full of NaN or Inf samples.
 * So stats might have its "inf" count adjusted. And only that. For
 * the same reason, "some_empty_instance+=stats" is not the same as
 * "some_empty_instance=stats". Furthermore, Adding a StatistisData
 * with all zero counts and min<=max might extend the range of the
 * target.
 *
 * While a StatisticData cannot shrink its min/max range itself,
 * higher level code that considers both the histogram and the
 * statistics might be able to do so. See trimRange().
 *
 * For defensive coding, I suggest using "stats-=old;stats+=new"
 * instead of "stats+=(new-old)" And use "stats=StatisticData()"
 * if the intent is to empty the instance.
 *
 * If invariants are violated (e.g. empty && count>0) then the result
 * of any operation is unspecified.
 *
 * NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM removes the following issue.
 * Special handling might be needed for single-value instances.
 * A single-value instance has had "count" values added, all of which
 * are the same. Statistics will have min==max. The histogram will
 * have a small bogus range to avoid division by zero. That bogus
 * range should preferably be reset once a second value is seen.
 *
 * NOTE: THE FOLLOWING DESCRIBES THE PREVIOUS, BUGGY VERSION.
 * ONLY INCLUDED TO HELP DEBUG ANY REGRESSIONS.
 *
 * The code used to ignore all fields in "other" when other.cnt
 * is zero, and assume this->min_ and this->max_ are bogus if this->cnt_
 * is zero. Both those cases lead to some very subtle and unexpected
 * behavior. Consider an incremental change:
 *
 *    - stats -= added; stats += removed; // Safe if count >= 0.
 *    - added -= removed; stats += added; // WRONG (2)
 *    - stats -= removed; stats += added; // Slightly wrong (3)
 *
 * In (2) the net count of changed objects will typically be zero and
 * operator+= will incorrectly ignore changes to sum, ssq, etc. (3)
 * can also fail in obscure cases if the caller removes and adds more
 * data than actually exists. E.g. subtracting also the samples in the
 * padding area and adding them back later, in the belief that this
 * gives the same result. Because the first -= might give a temporary
 * result with count zero. Another problem is that the code relies on
 * testing count==0 instead of the conventional min <= max to
 * determine whether the range is valid or not. Which is why (3) might
 * fail, but is also a problem in its own right.
 *
 * After the change, a StatisticData with all zeros is no longer a no-op
 * when used either as lhs or rhs of operator+= and operator-=.
 * This will now expand the range of lhs to include zero. I have changed
 * the default constructor to set min>max to bring back the nop-op
 * behavior. There is some risk that code creates an all-zero instance
 * by hand. Or worse, an instance with count zero and sum etc. garbage.
 *
 * END DESCRIBING THE OLD AND BUGGY HANDLING.
 */
class OPENZGY_TEST_API StatisticData
{
public:
  typedef std::int64_t count_type;
  count_type  getcnt() const { return cnt_; }  /**< \brief Net number of added samples. */
  count_type  getinf() const { return inf_; }  /**< \brief Net number of not added (infinite) samples. */
  double      getsum() const { return sum_; }  /**< \brief Sum of added samples. */
  double      getssq() const { return ssq_; }  /**< \brief Sum-of-squares of added samples. */
  double      getmin() const { return min_; }  /**< \brief Minimum added sample value. */
  double      getmax() const { return max_; }  /**< \brief Maximum added sample value. */
  bool        empty()  const { return min_ > max_; }  /**< \brief true if no values added. */
  bool        single() const { return min_ == max_; } /**< \brief true if all values are the same. */

  StatisticData();
  StatisticData(count_type cnt, count_type inf, double sum, double ssq, double min, double max);
  StatisticData(const count_type *bins, int nbins, double range_min, double range_max, bool is_int8);

  inline void add(double value);

  StatisticData& operator+=(const StatisticData& other);
  StatisticData& operator-=(const StatisticData& other);
  StatisticData& operator*=(count_type factor);

  inline void scale(double oldmin, double oldmax, double newmin, double newmax);
  void trimRange(const count_type *bins, int nbins, double range_min, double range_max);
  std::string toString() const;

private:
  count_type  cnt_;  /**< Number of added samples. */
  count_type  inf_;  /**< Number of not added (infinite) samples. */
  double      sum_;  /**< Sum of added samples. */
  double      ssq_;  /**< Sum-of-squares of added samples. */
  double      min_;  /**< Minimum added sample value. */
  double      max_;  /**< Maximum added sample value. */
};

/**
 * Add a single sample to the statistics.
 */
inline void
StatisticData::add(double value)
{
  if (std::isfinite(value)) {
    sum_ += value;
    ssq_ += value*value;
    if (cnt_ == 0)
      min_ = max_ = value;
    if (value < min_)
      min_ = value;
    if (value > max_)
      max_ = value;
    ++cnt_;
  }
  else
    ++inf_;
}

/**
 * Calculate the linear transform needed to convert from one range
 * (typically the natural data range of the integral storage type)
 * to the data range that the application wants to see.
 * Then update the statistics so they look like the transform had
 * been done on every single data point before adding it.
 *
 * Note: If caler has (slope, intercept) instead of the 4 factors,
 * pass oldmin=0, oldmax=1, newmin=intercept, newmax=intercept+slope.
 */
inline void StatisticData::scale(double oldmin, double oldmax, double newmin, double newmax)
{
  /*
  The decoded value Y is given by a linear transform of the coded value X:

    Y = a + b*X

  where a and b are given by the coding range and the value range of type T (see
  below). The statistics of Y are then:

    SUM_Y = SUM(a + b*x)
          = n*a + b*SUM(x) = n*a + b*SUM_X

    SSQ_Y = SUM((a + b*x)^2)
          = SUM(a^2 + 2*a*b*x + b^2*x^2)
          = n*a^2 + 2*a*b*SUM(x) + b^2*SUM(x^2)
          = n*a^2 + 2*a*b*SUM_X + b^2*SSQ_X

    MIN_Y = MIN(a + b*x)
          = a + b*MIN(x)
          = a + b*MIN_X

  and similar for MAX_Y.
  */

  const double b = (newmax - newmin) / (oldmax - oldmin);
  const double a = newmin - oldmin*b;

  ssq_ = cnt_*a*a + 2*a*b*sum_ + b*b*ssq_;
  sum_ = cnt_*a + b*sum_;
  min_ = a + b*min_;
  max_ = a + b*max_;
}

} // end namespace
