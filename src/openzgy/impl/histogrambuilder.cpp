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

//Adapted from: Zgy/UtilityLib/HistogramBuilder

#include "histogrambuilder.h"

#include <limits>
#include <assert.h>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Create a new HistogramBuilder.
 */
HistogramBuilder::HistogramBuilder(size_type nbins, double min, double max)
: stats_(), histo_(nbins, min, max)
{
}

/**
 * Copy constructor from a histogram passed as discrete information.
 * Useful for converting from some other histogram type.
 * Always create a fixed-range histogram with the same range
 * and number of bins as the source. If this is not what you need,
 * convert it using operator+=().
 * The above, below, and infinite counts are left at zero.
 */
HistogramBuilder::HistogramBuilder(const count_type* bins, int nbins, double min, double max, count_type scnt, double ssum, double sssq, double smin, double smax)
: stats_(scnt, 0, ssum, sssq, smin, smax), histo_(bins, nbins, min, max)
{
}

/**
 * Add the samples found in another histogram,
 * just as if the samples had been added one at a time using Add().
 */
HistogramBuilder& HistogramBuilder::operator+=(const HistogramBuilder& other)
{
  // Update both the histogram and the statistics in parallel.
  histo_ += other.gethisto();
  stats_ += other.getstats();

  // Try to narrow the statistics range based on the now updated histogram.
  // This is also needed with operator+= because we can add a "negative" statistic.
  stats_.trimRange(gethisto().getbins(), gethisto().getsize(), gethisto().getmin(), gethisto().getmax());
  return *this;
}

/**
 * Subtract the samples found in another histogram,
 * more or less undoing the effect of Add.
 * The min/max range in the statistics might be left
 * showing a too wide range.
 */
HistogramBuilder& HistogramBuilder::operator-=(const HistogramBuilder& other)
{
  // Update both the histogram and the statistics in parallel.
  histo_ -= other.gethisto();
  stats_ -= other.getstats();

  // Try to narrow the statistics range based on the now updated histogram.
  stats_.trimRange(gethisto().getbins(), gethisto().getsize(), gethisto().getmin(), gethisto().getmax());
  return *this;
}

/**
 * Multiply HistogramBuilder with a constant N, equivalent to creating
 * a new instance and adding the old one to it N times. N can be negative.
 */
HistogramBuilder& HistogramBuilder::operator*=(count_type factor)
{
  stats_ *= factor;
  histo_ *= factor;
  return *this;
}

/**
 * Two histograms are considered equal if they return the same bin count
 * for any input value. In practice there is some slop, as only the center
 * of each bin of either histogram is checked. The statistics information
 * need not match. Nor is there any check that the histograms have the
 * same range or even the same number of bins.
 */
bool HistogramBuilder::operator==(const HistogramBuilder& other) const
{
  return this->histo_ == other.histo_;
}

/**
 * See operator== for a description.
 */
bool HistogramBuilder::operator!=(const HistogramBuilder& other) const
{
  return this->histo_ != other.histo_;
}

/**
 * Calculate statistics directly from the histogram.
 * See the overloaded StatisticData constructor for details.
 * Here we don't know whether the histogram was built from
 * 8-bit data, so we play it safe and assume it wasn't.
 *
 * TODO-Low: YAGNI: no longer used because it gives incorrect
 * results when the histogram range cannot grow. Remove,
 * along with the other remnants of automatically expanding
 * histogram range.
 */
StatisticData
HistogramBuilder::gethiststats() const
{
  return StatisticData(gethisto().getbins(), gethisto().getsize(), gethisto().getmin(), gethisto().getmax(), false);
}

/**
 * Calculate the linear transform needed to convert from one range
 * (typically the natural data range of the integral storage type)
 * to the data range that the application wants to see.
 * Then update the histogram and associated statistics so they look
 * like the transform had been done on every single data point
 * before adding it.
 */
void
HistogramBuilder::scale(double oldmin, double oldmax, double newmin, double newmax)
{
  stats_.scale(oldmin, oldmax, newmin, newmax);
  histo_.scale(oldmin, oldmax, newmin, newmax);
}

/**
 * Erase all values from the histogram.
 * The range remains unchanged.
 */
void
HistogramBuilder::clear()
{
  histo_.clear();
  stats_ = StatisticData();
}

} // end namespace
