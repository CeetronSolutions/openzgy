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
#include <iomanip>

#define XXX_DISABLE_NAN_CHECK 0

namespace InternalZGY {
#if 0
}
#endif

/**
 * Create a new HistogramData.
 */
HistogramData::HistogramData(size_type nbins, double min, double max)
  : min_(min)
  , max_(max)
  , nbins_(nbins<2 ? 2 : nbins)
  , bins_(new count_type[nbins<2 ? 2 : nbins])
  , fixed_(true)
{
  clear();
  if (fixed_ && min_ >= max_)
    throw std::runtime_error("Bad histogram range in fixed histogram");
}

/**
 * Create a new HistogramData.
 */
HistogramData::HistogramData(size_type nbins, double min, double max, bool fixed)
  : min_(min)
  , max_(max)
  , nbins_(nbins<2 ? 2 : nbins)
  , bins_(new count_type[nbins<2 ? 2 : nbins])
  , fixed_(fixed)
{
  clear();
  if (fixed_ && min_ >= max_)
    throw std::runtime_error("Bad histogram range in fixed histogram");
}

/**
 * Copy constructor from a histogram passed as discrete information.
 * Useful for converting from some other histogram type.
 * Always create a fixed-range histogram with the same range
 * and number of bins as the source. If this is not what you need,
 * convert it using operator+=().
 */
HistogramData::HistogramData(const count_type* bins, int nbins, double min, double max)
  : min_(min)
  , max_(max)
  , nbins_(nbins<2 ? 2 : nbins)
  , bins_(new count_type[nbins<2 ? 2 : nbins])
  , fixed_(true)
{
  clear();
  for (size_type binno = 0; binno < nbins && binno < nbins_; ++binno)
    bins_[binno] = bins[binno];
  if (fixed_ && min_ >= max_)
    throw std::runtime_error("Bad histogram range in copied fixed histogram");
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
    fixed_  = other.fixed_;
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
 * Given a sample value, return how many times this sample was found.
 */
HistogramData::count_type HistogramData::get(double value) const
{
  if (XXX_DISABLE_NAN_CHECK || std::isfinite(value)) {

    // Get the linear transform from sample value to bin number,
    double A = 0, B = 0;
    calculateConversionFactors(&A, &B);

    // calculate bin index.
    const int n(RoundD2I(A + B*value));
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
  const count_type* ptr = this->getbins();
  const count_type* end = ptr + this->getsize();
  while (ptr < end)
    result += *ptr++;
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
 * Given two linear transforms as y = A = B*x, return a single transform
 * for applying the first and then the second transform.
 */
void
HistogramData::combineLinearTransform(
     double *A, double *B,
     double A1, double B1,
     double A2, double B2)
{
  *A = A2 + A1*B2;
  *B = B1 * B2;
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
  if (max_ > min_) {
    *B = (nbins_ - 1)/(max_ - min_);
    *A = -min_ * *B;
  } else {
    // Invalid range. All input values map to an illegal bin (-1).
    *B = 0;
    *A = -1;
  }
}

/**
 * Add or subtract data from another set of histogram bins, instead of
 * adding data from individual samples. Caller needs to ensure that
 * the range is big enough.
 *
 * Danger, technically the NaN/Inf handling triggers undefined behavior.
 * Because those might get cast to integer.
 *
 * NaN samples should be skipped. Handling of samples outside range,
 * including +/- inf, is to skipped those as well.
 *
 * TODO-WIP-BrickedAPI: Revisit the clip or skip decision.
 *
 * Skip, the current choice, probably looks better and is marginally
 * more efficient. ExpandableHistogram never clips or skips because it
 * gets expanded to make sure all samples fit between the *center*
 * vaues of the first and last bin.
 *
 * Clip would have worked better if the histogram was only expanded
 * enough to fit between the extreme values of the histogram.
 * Otherwise, rounding issues would become even trickier.
 *
 * If it turns out that the caller needs to know the number of clipped
 * or skipped values, the caller should request a histogram that is
 * two bins wider. Have the loop below do clipping. The caller then
 * discards the first and last bin from the actual histogram to skip
 * the outside samples. Optionally add those bins to their neighbors
 * first, if clipping is desired instead of skipping. Caveat, this
 * trick probably won't work for fixed range histograms because there
 * isn't any obvious place to do the extra trimming at the end.
 *
 * See ExpandableBuilder::antiZcBinWidthFromRange().
 */
void HistogramData::addBins(
     const count_type* more_bins,
     int more_nbins,
     double more_min,
     double more_max,
     bool add)
{
  // Transform bin number in "more" to real values to bin number in "this".
  double A1{0}, A2{0}, B1{0}, B2{0}, A{0}, B{0};
  getLinearTransform(&A1, &B1, 0, more_nbins - 1, more_min, more_max);
  getLinearTransform(&A2, &B2, this->min_, this->max_, 0, this->nbins_ - 1);
  combineLinearTransform(&A, &B, A1, B1, A2, B2);

  // Paranoia. Cheap as long as it is done outside the loop.
  if (!IsFiniteT(A) || !IsFiniteT(B) || more_min > more_max || min_ > max_)
    return;

  if (more_nbins == this->nbins_ &&
      std::fabs(A) < 0.01 &&
      std::fabs(B*more_nbins - 1) < 0.01)
  {
    // Short cut when bin numbers map 1:1
    if (add) {
      for (size_type binno = 0; binno < this->nbins_; ++binno)
        this->bins_[binno] += more_bins[binno];
    }
    else {
      for (size_type binno = 0; binno < this->nbins_; ++binno)
        this->bins_[binno] -= more_bins[binno];
    }
  }
  else {
    const int sign {add ? 1 : -1};
    for (size_type more_binno = 0; more_binno < more_nbins; ++more_binno) {
      const int this_binno = (RoundD2I(A + B*more_binno));
      if (this_binno >= 0 && this_binno < this->nbins_)
        this->bins_[this_binno] += sign * more_bins[more_binno];
    }
  }
}

std::string
HistogramData::toString(bool verbose) const
{
  std::int64_t count{0};
  for (const std::int64_t* ptr = getbins(), *end = getbins() + getsize(); ptr < end; ++ptr)
    count += *ptr;
  std::stringstream ss;
  ss << "cnt " << getcount()
     << " (" << count << ")"
     << " min "  << float(getmin()) << " max " << float(getmax())
     << " in " << getsize() << " bins"
     << " width " << float(getbinwidth());
  if (verbose) {
    ss << "\n";
    const double width = getbinwidth();
    const double begin = getmin();
    const count_type* bins = getbins();
    const int binno_width = getsize() < 1000 ? 3 : getsize() < 10000 ? 4 : 5;
    for (int ii = 0; ii < getsize(); ++ii) {
      const double center = begin + ii * width;
      const double lo = begin + (ii - 0.5) * width;
      const double hi = begin + (ii + 0.5) * width;
      // Show non-empty bins, first, last, and the bin for real zero.
      if (bins[ii] != 0 || (lo<=0 && hi>=0) || ii == 0 || ii == getsize()-1) {
        ss << "  bin[" << std::setw(binno_width) << ii << "]: "
           << std::setw(5) << bins[ii]
           << " value " << center
           << std::showpos << " (" << lo << ".." << hi << std::noshowpos << ")\n";
      }
    }
  }
  return ss.str();
}

/**
 * Output the histogram in agr format, to be viewed by grace or xmgrace.
 *
 * This is very rudimentary and is only intended for debugging.
 * An obvious improvement would be to scale the X axis to the correct values
 * instead of just showing the bin number.
 */
void
HistogramData::toGraph(std::ostream& os) const
{
  std::int64_t xmax = std::max(10, this->getsize());
  std::int64_t ymax{10};
  for (const std::int64_t* ptr = getbins(), *end = getbins() + getsize(); ptr < end; ++ptr)
    ymax = std::max(ymax, std::abs(*ptr));
  // Filter out the highest, typically the zero bin
  std::int64_t ymax2{10};
  for (const std::int64_t* ptr = getbins(), *end = getbins() + getsize(); ptr < end; ++ptr)
    if (std::abs(*ptr) != ymax)
      ymax2 = std::max(ymax2, std::abs(*ptr));
  // The bin with the highest frequency is not allowed to be higher
  // than twice the second-highest frequency.
  std::int64_t real_ymax = ymax;
  if (ymax > 2*ymax2)
    ymax = 2*ymax2;

  os << "# Grace project file\n"
     << "#\n"
     << "# Suggestions for viewing the result:\n"
     << "#   grace -hardcopy -hdevice PNG mychart.agr && display mychart.png\n"
     << "#   grace -hardcopy -hdevice SVG mychart.agr && inkscape mychart.svg\n"
     << "#   xmgrace mychart.agr\n"
     << "#\n"
     << "# Real ymax " << real_ymax << " real second " << ymax2 << "\n"
     << "# Using ymax " << ymax << "\n"
     << "# Total sample count " << this->getcount() << "\n"
     << "@version 50002\n"
     << "@page size 1600 800\n"
     << "\n"
     << "@g0 on\n"
     << "@g0 type Chart\n"
     << "@with g0\n"
     << "@    world xmin 0\n"
     << "@    world xmax " << xmax + 1<< "\n"
     << "@    world ymin " << -ymax/2 << "\n" // in case some negative
     << "@    world ymax " << ymax + 1<< "\n"
    //<< "@    world ymax 100 \n"
     << "@    view xmin 0.1\n"
     << "@    view xmax 1.8\n"
     << "@    view ymin 0.1\n"
     << "@    view ymax 0.9\n"
     << "\n"
     << "@    xaxis  on\n"
     << "@    xaxis  type zero false\n"
     << "@    xaxis  tick on\n"
     << "@    xaxis  tick major " << xmax / 8 << "\n"
     << "@    xaxis  tick out\n"
     << "@    xaxis  tick op left\n"
     //<< "@    xaxis  tick type auto\n"
     //<< "@    xaxis  ticklabel type auto\n"
     << "@    xaxis  tick spec type both\n"
     << "\n";
  os << "@    xaxis  tick spec 17\n";
  for (int ii=0; ii<17; ++ii) {
    int step = getsize() / 16;
    int binno = ii * step;
    double real = getmin() + binno * getbinwidth();
    if ((ii%2) == 0)
      os << "@    xaxis  tick major " << ii << ", " << binno << "\n"
         << std::setprecision(4)
         << "@    xaxis  ticklabel " << ii << ", \"" << real << "\"\n";
    else
      os << "@    xaxis  tick minor " << ii << ", " << binno << "\n";
  }
  os << "@    yaxis  on\n"
     << "@    yaxis  type zero false\n"
     << "@    yaxis  tick on\n"
     << "@    yaxis  tick major " << ymax / 5 << "\n"
     << "@    yaxis  tick out\n"
     << "@    yaxis  tick op left\n"
     << "@    yaxis  ticklabel type auto\n"
     << "@    yaxis  tick type auto\n"
     << "\n"
     << "@    s0 type bar\n"
     << "@    s0 symbol 0\n"
     << "@    s0 symbol color 2\n"
     << "@    s0 symbol pattern 1\n"
     << "@    s0 symbol fill color 2\n"
     << "@    s0 symbol fill pattern 1\n"
     << "@    s0 symbol linestyle 0\n"
     << "@    s0 symbol linewidth 0\n"
     << "@    s0 symbol size 0.2\n"
     << "@    s0 line type 0\n"
     << "\n"
     << "@target G0.S0\n"
     << "@type bar\n";
  for (int ii = 0; ii < this->getsize(); ++ii) {
    std::int64_t value = getbins()[ii];
    //os << ii << " " << value << "\n";
    os << ii << " " << (std::min(value, ymax)) << "\n";
    //os << ii << " " << (std::min(value, ymax)) / (ymax / 100.0) << "\n";
  }
  os << "&\n";
}

} // end namespace
