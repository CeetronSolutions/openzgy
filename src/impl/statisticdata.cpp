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

#include "statisticdata.h"

#include <cstdint>
#include <limits>
#include <assert.h>
#include <algorithm>
#include <sstream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Initialize a StatisticData to empty.
 * min/max are by definition irrelevant when count is zero,
 * but they will be set to 0 here just to have a consistent value.
 */
StatisticData::StatisticData()
  : cnt_(0), inf_(0), sum_(0), ssq_(0)
  , min_(std::numeric_limits<double>::infinity())
  , max_(-std::numeric_limits<double>::infinity())
{
}

/**
 * Initialize a StatisticData from discrete values.
 */
StatisticData::StatisticData(count_type cnt, count_type inf, double sum, double ssq, double min, double max)
  : cnt_(cnt), inf_(inf), sum_(sum), ssq_(ssq), min_(min), max_(max)
{
}

/**
 * Create a StatisticData from histogram information.
 *
 * If is_int8 is true, the histogram was populated from 8-bit data
 * and the histogram has 256 bins. In this case the result will be
 * precise. The method needs to be told about this. Otherwise it
 * would have to add 1/2 bin width slop, not knowing that there can
 * only be a single value in each bin. That value might not be
 * integral any longer due to scaling, buit there can still only be one.

 * If is_int8 is false, this isn't as accurate as the actual statistics
 * since each sample will be rounded to the center of the bin it falls in.
 * Also, Inf will always be 0 since the histogram doesn't count outliers.
 * The result is still useful as a consistency check. And also to update
 * the min/max range in the real statistics after subtracting samples.
 *
 * Used by HistogramBuilder::fastAddSignedChar256() which provides a
 * not insignificant performance boost in int8 files.
 *
 * Used by HistogramBuilder::gethiststats(), but that method itself is no
 * longer invoked because it gives incorrect results when the histogram
 * range cannot grow.
 */
StatisticData::StatisticData(const count_type *bins, int nbins, double range_min, double range_max, bool is_int8)
: cnt_(0), inf_(0), sum_(0), ssq_(0)
, min_(std::numeric_limits<double>::infinity())
, max_(-std::numeric_limits<double>::infinity())
{
  const double begin = range_min;
  const double width = (range_max - range_min) / (nbins - 1);
  const double slop = is_int8 ? 0 : width/2;

  for (int binno = 0; binno < nbins; ++binno) {
    count_type count = bins[binno];
    double value = begin + binno*width;
    if (count != 0) {
      if (cnt_ == 0)
        min_ = value - slop; // first non-empty bin found
      max_ = value + slop; // last (so far) non-empty bin found.
    }
    cnt_ += count;
    sum_ += count * value;
    ssq_ += count * value * value;
  }
}

/**
 * Add more samples to an existing StatisticData.
 * other is also allowed to hold negative counts,
 * this will cause samples to be removed. But the
 * min/max range will still be expanded.
 *
 * CAVEAT: The code used to ignore all fields in "other" when other.cnt
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
 * bu hand. Or worse, an instance with count zero and sun etc. garbage.
 */
StatisticData& StatisticData::operator+=(const StatisticData& other)
{
  if (other.min_ <= other.max_) {
    if (min_ <= max_) {
      min_ = std::min(min_, other.min_);
      max_ = std::max(max_, other.max_);
    }
    else {
      min_ = other.min_;
      max_ = other.max_;
    }
  }
  cnt_ += other.cnt_;
  sum_ += other.sum_;
  ssq_ += other.ssq_;
  inf_ += other.inf_;
  return *this;
}

/**
 * Remove samples from an existing StatisticData.
 * The resulting min/max range will be the union of the two classes,
 * not the intersection. Could also have chosen to leave the
 * min/max uncanged. But it is more consistent to do a union,
 * because -a+b will then work like b-a.  Either way, you can call
 * trimRange() afterwards if you know what the range ought to be.
 */
StatisticData& StatisticData::operator-=(const StatisticData& other)
{
  if (other.min_ <= other.max_) {
    if (min_ <= max_) {
      min_ = std::min(min_, other.min_);
      max_ = std::max(max_, other.max_);
    }
    else {
      min_ = other.min_;
      max_ = other.max_;
    }
  }
  cnt_ -= other.cnt_;
  sum_ -= other.sum_;
  ssq_ -= other.ssq_;
  inf_ -= other.inf_;
  return *this;
}

/**
 * Multiply StatisticData with a constant N, equivalent to creating
 * a new instance and adding the old one to it N times. N can also be
 * negative. The min/max range is not affected.
 */
StatisticData& StatisticData::operator*=(count_type factor)
{
  cnt_ *= factor;
  inf_ *= factor;
  sum_ *= factor;
  ssq_ *= factor;
  return *this;
}

/**
 * Try to narrow the range found in the statistics based on the
 * information we have in the histogram. This is particularly
 * useful after subtracting something from the histogram.
 *
 * The range found in the histogram is normally wider than or equal to
 * the statistics range. If this is not the case, we know the min/max
 * in the statistics are incorrect.  We cannot know the precise new
 * range without rescanning all the data, but we can approximate it
 * using the histogram.  There is no point in doing any of this if the
 * histogram is empty. In that case min/max are irrelevant, and will
 * be reset once data has been added.  If any of the overflow bins are
 * in use, these are ignored.  That should not happen in normal
 * circumstances, but could occur if stats are calculated from
 * application provided data, but maintained in a fixed histogram
 * because that is what the file format dictates.
 */
void StatisticData::trimRange(const count_type *bins, int nbins, double range_min, double range_max)
{
  if (cnt_ != 0) {
    const double begin = range_min;
    const double width = (range_max - range_min) / (nbins - 1);
    const double slop = width/2;

    int lobin = 0;
    while (bins[lobin] == 0 && lobin < nbins)
      ++lobin;
    int hibin = nbins - 1;
    while (bins[hibin] == 0 && hibin >= 0)
      --hibin;
    if (lobin <= hibin) {
      // test is needed because all bins might have been empty,
      // due to an inconsistency or because all the values were
      // in the overflow bins or NaN (currently not counted).
      // Not much we can do if that is the case.
      // Applying the transform gives the center of the bin,
      // the actual range of the data can be 1/2 bin more at either end.
      double histmin = begin + width*lobin - slop;
      double histmax = begin + width*hibin + slop;
      if (min_ < histmin)
        min_ = histmin;
      if (max_ > histmax)
        max_ = histmax;
    }
  }
}

std::string
StatisticData::toString() const
{
  if (getcnt() == 0 && getsum() == 0 && getssq() == 0 && getinf() == 0 && getmin() > getmax())
    return "empty";
  std::stringstream ss;
  ss << "min " << getmin() << " max " << getmax() << " cnt " << getcnt()
     << " sum " << getsum() << " ssq " << getssq()
     << " avg " << (getcnt() ? getsum() / getcnt() :
                    std::numeric_limits<double>::quiet_NaN());
  return ss.str();
}

} // end namespace
