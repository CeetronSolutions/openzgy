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

//Adapted from: Zgy/UtilityLib/HistogramData

#include "histogramdata.h"

#include <limits>
#include <assert.h>
#include <memory.h>
#include <cmath>
#include <sstream>

#define XXX_DISABLE_NAN_CHECK 0

namespace InternalZGY {
#if 0
}
#endif

/**
 * Create a new HistogramData.
 */
HistogramData::HistogramData(size_type nbins, double min, double max)
: min_(min), max_(max), nbins_(nbins<2 ? 2 : nbins), bins_(new count_type[nbins<2 ? 2 : nbins])
{
  clear();
}

/**
 * Copy constructor from a histogram passed as discrete information.
 * Useful for converting from some other histogram type.
 * Always create a fixed-range histogram with the same range
 * and number of bins as the source. If this is not what you need,
 * convert it using operator+=().
 */
HistogramData::HistogramData(const count_type* bins, int nbins, double min, double max)
: min_(min), max_(max), nbins_(nbins<2 ? 2 : nbins), bins_(new count_type[nbins<2 ? 2 : nbins])
{
  clear();
  for (size_type binno = 0; binno < nbins && binno < nbins_; ++binno)
    bins_[binno] = bins[binno];
}

/**
 * Standard copy constructor.
 * Needs to be explicit because of the allocated data for bins.
 */
HistogramData::HistogramData(const HistogramData& other)
: bins_(nullptr)
{
  *this = other;
}

/**
 * Standard assignment operator.
 * Needs to be explicit because of the allocated data for bins.
 */
HistogramData& HistogramData::operator=(const HistogramData& other)
{
  if (this != &other) {
    // Handle allocated data.
    nbins_  = other.nbins_;
    bins_.reset(new count_type[nbins_]);
    for (int ii=0; ii<nbins_; ++ii)
      bins_[ii] = other.bins_[ii];

    // The rest is simple copying
    min_    = other.min_;
    max_    = other.max_;
  }

  return *this;
}

/**
 * Add the samples found in another histogram,
 * just as if the samples had been added one at a time using Add().
 */
HistogramData& HistogramData::operator+=(const HistogramData& other)
{
  // Add the new data to our bins. The range doesn't need to match.
  addBins(other.getbins(), other.getsize(), other.getmin(), other.getmax(), true/*add*/);
  return *this;
}

/**
 * Subtract the samples found in another histogram,
 * more or less undoing the effect of add.
 * The min/max range in the statistics might be left
 * showing a too wide range.
 */
HistogramData& HistogramData::operator-=(const HistogramData& other)
{
  addBins(other.getbins(), other.getsize(), other.getmin(), other.getmax(), false/*subtract*/);
  return *this;
}

/**
 * Multiply HistogramBuilder with a constant N, equivalent to creating
 * a new instance and adding the old one to it N times. N can be negative.
 */
HistogramData& HistogramData::operator*=(count_type factor)
{
  for (int ii=0; ii<nbins_; ++ii)
    bins_[ii] *= factor;
  return *this;
}

/**
 * Two histograms are considered equal if they return the same bin count
 * for any input value. In practice there is some slop, as only the center
 * of each bin of either histogram is checked. The statistics information
 * need not match. Nor is there any check that the histograms have the
 * same range or even the same number of bins.
 */
bool HistogramData::operator==(const HistogramData& other) const
{
  return this->compare(other) && other.compare(*this);
}

/**
 * See operator== for a description.
 */
bool HistogramData::operator!=(const HistogramData& other) const
{
  return !(this->compare(other) && other.compare(*this));
}

/**
 * Given a sample value, return how many times this sample was found.
 * This is almost the inverse of AddOne, but the function does no
 * scaling of the value.
 */
HistogramData::count_type HistogramData::get(double value) const
{
  if (XXX_DISABLE_NAN_CHECK || std::isfinite(value)) {

    // Get the linear transform from sample value to bin number,
    double A = 0, B = 0;
    calculateConversionFactors(&A, &B);

    // calculate bin index.
    int n(RoundD2I(A + B*value));
    if (n >= 0 && n < nbins_)
      return bins_[n];
    else
      return 0;
  } else {
    return 0;
  }
}

/**
 * Return the total number of samples added to any bin.
 * Currently this is returned by summing all the bins,
 * if this is a performance problem that it is possible
 * to maintain a separate count member.
 * Called from compare().
 */
HistogramData::count_type HistogramData::getcount() const
{
  count_type result = 0;
  for (size_type ii = 0; ii < this->getsize(); ++ii)
    result += this->getbins()[ii];
  return result;
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
HistogramData::scale(double oldmin, double oldmax, double newmin, double newmax)
{
  // Same logic as in StatisticData::scale() to handle our min/max.
  const double b = (newmax - newmin) / (oldmax - oldmin);
  const double a = newmin - oldmin*b;
  min_ = a + b*min_;
  max_ = a + b*max_;
}

/**
 * Erase all values from the histogram.
 * The range remains unchanged.
 */
void HistogramData::clear()
{
  memset(bins_.get(), 0, nbins_*sizeof(bins_[0]));
}

/**
 * Get the linear transform (offset, scale) that will map values from the range
 * oldmin..oldmax to the range newmin..newmax. If either range is empty,
 * the identity transform will be returned.
 */
void
HistogramData::getLinearTransform(double *offset, double *scale, double oldmin, double oldmax, double newmin, double newmax)
{
  // Calculate a linear transform that will convert oldmin to newmin and oldmax to newmax
  if ((newmax <= newmin) || (oldmax <= oldmin)) {
    *scale = 1;
    *offset = 0;
  } else {
    *scale = (newmax-newmin) / (oldmax-oldmin);
    *offset = newmin - *scale * oldmin;
  }
}

/**
 * Get the linear transform for mapping from application values to bin numbers.
 * The resulting bin number is meant to be rounded to nearest integral value.
 * If the histogram range is not initialized yet (can happen if histogram is empty)
 * then a transform is returned that will map all values to a non-existant bin.
 */
void HistogramData::calculateConversionFactors(double* A, double* B) const
{
  // Get linear transform for mapping input range [min_ .. max_] to bin numbers [0 .. nbins_-1]
  // return getLinearTransform(A, B, min_, max_, 0, nbins_-1);
  // However it is better to inline the code because we want different error handling.
  if (max_ <= min_) {
    // Invalid ramge. All input values map to an illegal bin (-1).
    *B = 0;
    *A = -1;
  } else {
    *B = (nbins_ - 1)/(max_ - min_);
    *A = -min_ * *B;
  }
}

/**
 * Add or subtract data from another set of histogram bins,
 * instead of adding data from individual samples.
 * Caller needs to ensure that the range is big enough.
 */
void HistogramData::addBins(const count_type* bins, int nbins, double min, double max, bool add)
{
  // We aren't going to change the range again in the loop below.
  // So the conversion factors only need to be calculated once.
  double A=0, B=0;
  calculateConversionFactors(&A, &B);

  // Conversion needed in other, to go in the opposite direction
  // i.e. from bin number to actual value.
  double begin = min;
  double width = (max - min) / (nbins - 1);

  for (size_type binno = 0; binno < nbins; ++binno) {
    double value = begin + binno*width;
    count_type count = bins[binno];
    addOne(value, add ? count : -count, A, B);
  }
}

/**
 * Check all our bins, to see that the other histogram has the exact
 * same count when looking up the value that is the center point of
 * our own bins. Note that the other histogram could still have more
 * samples than we do, if it has a wider range or more densly sampled
 * bins. When comparing two histograms you probably want to compare
 * in both directions. Or at least do an extra check on the total
 * sample count - but I prefer the double check as it is guaranteed
 * to be symmetric.
 */
bool HistogramData::compare(const HistogramData& other) const
{
  // special handling needed for empty histograms.
  // If histogram is not empty, we are guaranteed
  // that width is not zero. For empty one it might.

  count_type self_count = this->getcount();
  count_type other_count = other.getcount();

  if (self_count == 0 || other_count == 0)
    return self_count == other_count; // true only if both are empty.

  double width = (this->getmax() - this->getmin()) / (this->getsize() - 1);
  double begin = this->getmin();
  for (int ii=0; ii<this->getsize(); ++ii) {
    double value = begin + ii*width;
    if (this->getbins()[ii] != other.get(value))
      return false;
  }
  return true;
}

std::string
HistogramData::toString(bool verbose) const
{
  std::int64_t count{0};
  for (const std::int64_t* ptr = getbins(), *end = getbins() + getsize(); ptr < end; ++ptr)
    count += *ptr;
  std::stringstream ss;
  ss << "min "  << getmin() << " max " << getmax() << " cnt " << getcount()
     << " (" << count << ")";
  if (verbose)
    for (const std::int64_t* ptr = getbins(), *end = getbins() + getsize(); ptr < end; ++ptr)
      if (*ptr != 0)
        ss << " [" << (ptr - getbins()) << "]:" << *ptr;
  return ss.str();
}

} // end namespace
