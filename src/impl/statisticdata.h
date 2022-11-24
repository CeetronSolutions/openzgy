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
 * \brief Holds the result of computing statistics.
 *
 * \details Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_TEST_API StatisticData
{
public:
  typedef std::int64_t count_type;
  count_type  getcnt() const { return cnt_; }  /**< \brief Number of added samples. */
  count_type  getinf() const { return inf_; }  /**< \brief Number of not added (infinite) samples. */
  double      getsum() const { return sum_; }  /**< \brief Sum of added samples. */
  double      getssq() const { return ssq_; }  /**< \brief Sum-of-squares of added samples. */
  double      getmin() const { return min_; }  /**< \brief Minimum added sample value. */
  double      getmax() const { return max_; }  /**< \brief Maximum added sample value. */

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
