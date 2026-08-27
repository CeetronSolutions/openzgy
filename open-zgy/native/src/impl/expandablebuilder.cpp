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


#include "expandablebuilder.h"

#include <limits>
#include <cmath>
#include <sstream>
#include <iostream>

namespace InternalZGY {
#if 0
}
#endif

// NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM does not allow explicitly
// specifying the histogram range. There are still unit tests
// doing this. Eventually this constructor should be removed.
ExpandableBuilder::ExpandableBuilder(
     size_type nbins,
     double min,
     double max,
     bool fixed)
  : HistogramBuilder(nbins, min, max, fixed)
{
  // TODO-WIP-BrickedAPI: There are unit tests that create builders that
  // violate the assumptions. Should be removed or fixed.
  //if (!fixed && (max-min <= 0 || min != -max))
  //  throw std::runtime_error("Wrong limits for dynamic HistogramData");
}

ExpandableBuilder::ExpandableBuilder(size_type nbins)
  : HistogramBuilder(nbins,
                     smallestAntiZcRange(nbins, /*center=*/true).first,
                     smallestAntiZcRange(nbins, /*center=*/true).second,
                     /*fixed=*/false)
{
}

/**
 * Create a fixed-width builder, like a plain HistogramBuilder.
 */
ExpandableBuilder::ExpandableBuilder(size_type nbins, double min, double max)
  : HistogramBuilder(nbins, min, max)
{
}

/**
 * Copy a builder, keeping the histogram size, limits, and isfixed.
 * Optionally copy the bin contents and the statistics.
 *
 * TODO-WIP-BrickedAPI: Is this dangerous? If creating a dynamic
 * histogram this way, there is no check for it being symmetric
 * and anti-zero-centric.
 */
ExpandableBuilder::ExpandableBuilder(
     const HistogramData& histo,
     const StatisticData& stats,
     bool copy)
  : HistogramBuilder(histo.getsize(),
                     histo.getmin(),
                     histo.getmax(),
                     histo.isfixed())
{
  if (copy) {
    this->histo_ = histo; // Minor: bins allocated twice.
    this->stats_ = stats;
  }
}

/**
 * \brief Static version of operator+= and operator-=.
 *
 * \details
 * To be used by application code that only holds stats and histo pointers
 * making it impractical and expensive to instanciate a builder and later
 * copy back the stats and histo.
 *
 * This method violates encapsulation since it includes one line of code
 * from HistogramBuilder::operator+=. Also, it would perhaps be better
 * if our operator+= and operator-= just called this function.
 */
void /*static*/
ExpandableBuilder::staticAddOrSub(
       StatisticData* stats,
       HistogramData* histo,
       const StatisticData& other_stats,
       const HistogramData& other_histo,
       bool add)
{
  SimpleTimerEx tt(HistoBuilderAdHocTimers::instance().addorsub);
  // Code from ExpandableBuilder::operator+= and HistogramBuilder::operator+=
  // The logic here MUST BE KEPT IN SYNC with those two.
  widenRange(stats, histo, other_stats.getmin(), other_stats.getmax());
  if (add) {
    *stats += other_stats;
    *histo += other_histo;
  }
  else {
    *stats -= other_stats;
    *histo -= other_histo;
  }
  // Might want to remove this feature for more reproducible results.
  // In fact, trackChanges() already bypassed this step because
  // it used += on the lower level types. Which is probably better,
  // as trimRange() will give wrong results if saples were clipped
  // from a fixed range histogram.
  //StatisticData prev(*stats);
  //stats->trimRange(histo->getbins(),
  //                 histo->getsize(), histo->getmin(), histo->getmax());
  //if (prev.getmin() != stats->getmin() || prev.getmax() != stats->getmax())
  //  std::cerr << "trimRange in " << prev.getmin() << " to " << prev.getmax()
  //            << " out " << stats->getmin() << " to " << stats->getmax()
  //            << std::endl;
  //else
  //  std::cerr << "trimRange no change "
  //            << stats->getmin() << " to " << stats->getmax() << std::endl;
}

/**
 * NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM removes the following issue.
 * Special handling might be needed for dynamic histograms if *this is
 * currently empty or single. The current histogram might in that case
 * have a too wide range. See the very long comment in the
 * StatisticData class.
 *
 * If other is empty or this is empty that this is almost a no-op
 * and a copy respectively. But stats_.inf_ can still be != 0.
 *
 * Note that "other" might also be the simpler, fixed-size type.
 */
ExpandableBuilder&
ExpandableBuilder::operator+=(const HistogramBuilder& other)
{
  widenRange(other.getstats().getmin(), other.getstats().getmax());
  HistogramBuilder::operator+=(other);
  return *this;
}

/**
 * See operator+= for details.
 * The special case foo -= foo might technically have cleared the
 * statistics range, but who would use that obscure feature?
 */
ExpandableBuilder&
ExpandableBuilder::operator-=(const HistogramBuilder& other)
{
  widenRange(other.getstats().getmin(), other.getstats().getmax());
  HistogramBuilder::operator-=(other);
  return *this;
}

ExpandableBuilder&
ExpandableBuilder::operator*=(count_type factor)
{
  HistogramBuilder::operator*=(factor); return *this;
}

double /*static*/
ExpandableBuilder::getsmallestbinwidth()
{
  // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM assumed here.
  // For expandable histograms, bin width is constrained to be 2^N.
  // It is also constrained by not being a denormal number, even when
  // stored as a float. And just to be safe, use a larger value.
  // Approximately 1.0e-30 should still be way less that any real value.
  //return (float)pow(2.0, -125.0); // technically this would work.
  //return (float)pow(2.0, -149.0); // even this, with denormals.
  return (float)pow(2.0, -100.0); // defensive coding.
}

std::pair<double,double> /*static*/
ExpandableBuilder::antiZcRangeFromBinWidth(size_type nbins, double binwidth, bool center)
{
  // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM assumed here. Also even bin count.
  // Technically I could probably handle odd number of bins, but why
  // bother adding that extra complexity? Also, doubleRange() currently
  // demands a multiple of 4.
  if ((nbins % 4) != 0 || nbins < 2)
    throw std::runtime_error("Internal error: Odd histogram range");

  // +/- zero maps to a boundary between two bins.
  // TODO-WIP-BrickedAPI: there are choices here.
  // * Use N/2-1 and N/2-1 to make a completely symmetric range.
  // * Do not use N/2-1 and N/2 because this would make -0 map to an
  //   odd bin number which might mess with the final collapse for making
  //   the range zero centric.
  // * Use N/2 and N/2+1 for slightly unsymmetric that might help keep
  //   the special -128..+127 look better. Maybe. That causes a lot of
  //   ripple effects, though.
  // TODO-WIP-BrickedAPI: Another choice, made by caller via "center" arg:
  // * The histogram only needs to be wide enough to have values min and max
  //   fit inside the extreme range of the histogram.
  // * The histogram needs to be wide enough to allow values min and max fit
  //   between the center of the first bin and the center of the last bin.
  // The choice might affect how numerical instability is handled. Which in
  // turn depends on whether a very slight tweak gets added.
  // But note: Won't the second option make sense anyway because all it
  // does is to maybe double the range once more than needed? If the
  // intermediate size is 4x of final sise, or more, this should be taken
  // care of when finalizing.
  //
  // There is also a choice about tweaking the parameters slightly
  // to avoid numerical instability, or rather, have the corner cases
  // not be "all integer samples".
  static const double fakemul   =  1.0 + pow(2.0, -30);
  static const double fakeshift = -0.5 * pow(2.0, -30) * binwidth;
  binwidth *= fakemul;
  double min = -(nbins/2+0) * binwidth + (center ? binwidth/2 : 0);
  double max =  (nbins/2-0) * binwidth - (center ? binwidth/2 : 0);
  return std::make_pair(min+fakeshift, max+fakeshift);
}

std::pair<double,double> /*static*/
ExpandableBuilder::smallestAntiZcRange(size_type nbins, bool center)
{
  return antiZcRangeFromBinWidth(nbins, getsmallestbinwidth(), center);
}

/**
 * Compute the smallest possible bin width that allows all finite samples
 * to fit inside the histogram.
 *
 * NOTE: Check whether expanding is needed before calling this function.
 * And/or check that the histogram range isn't inadvertently narrowed.
 * That would mess up existing samples.
 */
double /*static*/
ExpandableBuilder::antiZcBinWidthFromRange(size_type nbins, double min, double max)
{
  if (min > max || (min == 0 && max == 0))
    return getsmallestbinwidth();

  // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM assumed here. Also even bin count.

  double factor = std::max(std::fabs(min), std::fabs(max));
  double width = factor / (nbins/2);
  // Round up to a power of two. There is some risk of numerical
  // instability. And there are subtle issues with the range being not
  // precisely symmetrical. If the result is wrong, make sure it fails
  // in the "too small" direcction.
  //std::cout << "Width before round up: " << width << std::endl;
  width = std::pow(2.0, std::ceil(std::log(width*0.9)/std::log(2.0)));
  //std::cout << "Width after round up: " << width << std::endl;

  // Check the result. If too small, doubling the width should fix it.
  // This is also where we decide whether it is enough that min,max
  // falls inside the extreme range of the histogram, or whether it
  // needs to be between the centers of the first and last bin.
  // The first case might cause samples to be dropped due to numerical
  // instability. Unless there is an option to clip instead of drop.
  // The second case can, also due to numerical instability, cause
  // trouble for the -128..+127 special case.
  std::pair<double, double> range =
    antiZcRangeFromBinWidth(nbins, width, /*center=*/true/*?*/);
  if (min < range.first || max > range.second) {
    //std::cout << "checking ... need double\n";
    width *= 2;
  }
  else {
    //std::cout << "checking ... ok\n";
  }

  //std::cout << "WidthFromRange(" << nbins
  //          << ", range " << min << " to " << max
  //          << ") = " << width << "\n";

  // If later using antiZcRangeFromBinWidth() to create a new histogram,
  // remember that the constructor expects center=true regardless
  // of what was used here.
  return width;

  // TODO-WIP-BrickedAPI: Should I return e.g. width * tiny_number
  // to reduce issues with numerical inaccuracy?
  // Problematic for large histogram where the accumulated
  // "error" can get close to one bin width.
  //
  // Or consider a kludge in the sample-value-to-bin code
  // which rounds toward zero if the found bin is very
  // close to a boundary between bins.
}

/**
 * Return true if WidenRange will change anything, i.e. if the histogram
 * is not fixed and if the current range is too small to hold the data
 * in the supplied statistics.
 */
bool /*static*/
ExpandableBuilder::widenRangeIsNeeded(
       const StatisticData* stats,
       const HistogramData* histo,
       double new_min, double new_max)
{
  // A fixed-size histogram with a useless range is bad news.
  if (histo->isfixed() && histo->getmin() >= histo->getmax())
    throw std::runtime_error("Internal error: badhistogram range");
  if (histo->isfixed())
    return false; // Can never widen a fixed histogram.
  if (new_min > new_max)
    return false; // other builder has nothing to contribute.
  // TODO-WIP-BrickedAPI: Questionable. Histogram might already
  // have a range even when it and the statistics is empty.
  // Mitigate by adding  code in widenRange() to avoid inadvertently
  // narrowing the range.
  // Can I simply remove this instead, since a "not established"
  // range in the histogram is always -tiny..+tiny? If so, this is
  // dependant on SYMMETRIC_EXPANDABLE_HISTOGRAM.
  //if (stats->empty())
  //  return true;  // range has not been established yet.
  if (new_min < histo->getmin() || new_max > histo->getmax())
    return true;  // Current range too small
  return false;   // Current range is adequate.
}

bool
ExpandableBuilder::widenRangeIsNeeded(double new_min, double new_max) const
{
  return widenRangeIsNeeded(&stats_, &histo_, new_min, new_max);
}

/**
 * Make sure the histogram is wide enough to include the point(s) to add.
 * This is a no-op if the histogram was created as fixed.
 *
 * Choose the new range so that new_max fits between the center value
 * of the new histogram range. (center=true). Not just fit between the
 * extreme values. This simplifies issues with clipping and rounding,
 * and the only effect (I believe) is to sometimes set a range twice
 * as big as strictly needed. This ought to be fixed in finalize().
 *
 * NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM is an important assumption.
 */
void /*static*/
ExpandableBuilder::widenRange(
       StatisticData* stats,
       HistogramData* histo,
       double new_min, double new_max)
{
  if (widenRangeIsNeeded(stats, histo, new_min, new_max)) {
    //const double old_min{histo->getmin()}; // debugging
    //const double old_max{histo->getmax()}; // debugging
    const size_type nbins = histo->getsize();
    const double width = antiZcBinWidthFromRange(nbins, new_min, new_max);
    const std::pair<double,double> range =
      antiZcRangeFromBinWidth(nbins, width, /*center=*/true);
    HistogramData h(nbins, range.first, range.second, /*fixed=*/false);
    h += *histo;
    *histo = h;
    //std::cout << "Expand " << old_min << " " << old_max
    //          << " to fit " << new_min << " " << new_max
    //          << " ==> " << h.getmin() << " " << h.getmax() << "\n";
  }
  else {
    //std::cout << "NO expand "
    //          << histo->getmin() << " " << histo->getmax()
    //          << " to fit " << new_min << " " << new_max << "\n";
  }
}

void ExpandableBuilder::widenRange(double new_min, double new_max)
{
  return widenRange(&stats_, &histo_, new_min, new_max);
}

#if 0 // Used in commented-out ad-hoc debugging.
namespace {
  static std::string fmt(const HistogramBase& h)
  {
    const double width = (h.getmax() - h.getmin()) / (h.getsize() - 1);
    std::stringstream ss;
    ss << "HistogramBase size = " << h.getsize()
       << " range " << h.getmin() << " to " << h.getmax()
       << " extreme " << h.getmin() - width/2 << " to " << h.getmax() + width/2
       << " width " << width;
    return ss.str();
  }
}
#endif

/**
 * First step in preparing the histogram for use: Reduce the number of bins
 * by dropping empty buckets. Input is expandable, symmetric, anti-zero-centric.
 * Output is fixed, possibly non-symmetric, anti-zero-centric.
 * The output will be missing the actual bin counts.
 */
HistogramBase /*static*/
ExpandableBuilder::finalTrimEnds(const HistogramData& in)
{
  const size_type nbins = in.getsize();
  const double binwidth = (in.getmax() - in.getmin()) / (nbins - 1);
  size_type skipleft = 0;
  size_type skipright = 0;
  while (skipleft < nbins && in.getbins()[skipleft] == 0)
    ++skipleft;
  while (skipright < nbins && in.getbins()[nbins-1-skipright] == 0)
    ++skipright;
  if (skipleft + skipright >= nbins) // Empty. Leave size and range as is.
    skipleft = skipright = 0;
  else if (skipleft + skipright == nbins-1) { // Require histo size >= 2.
    if (skipleft != 0)
      --skipleft;
    else
      --skipright;
  }
  return HistogramBase(nbins - (skipleft + skipright),
                       in.getmin() + (skipleft * binwidth),
                       in.getmax() - (skipright * binwidth));
}

/**
 * Called before compressing the histogram by combining every "factor"
 * bins into one. The input must be anti-zero-centric. Slightly widen
 * the input histogram so both the number of bins with a negative
 * center value and the total number of bins are a multiple of factor.
 * This is done so the histogram after compression will still be
 * anti-zero-centric. Also make sure the effective "factor" doesn't
 * end up slightly not integral. If all sample values are positive
 * then the zero-centric property is not important. It might or might
 * not be kept.
 */
HistogramBase /*static*/
ExpandableBuilder::finalTweakFactor(const HistogramBase& in, size_type factor)
{
  const size_type left_bins  = (size_type)std::max(0.0, std::ceil(-(in.getmin() / in.getbinwidth())));
  const size_type right_bins = in.getsize() - left_bins;
  const size_type left_pad   = factor - (left_bins % factor);
  const size_type right_pad  = factor - (right_bins % factor);

  //std::cout << "factor " << factor
  //          << " left "  << left_bins << "+" << left_pad
  //          << " right " << right_bins << "+" << right_pad
  //          << "\n";

  return HistogramBase(in.getsize() + left_pad + right_pad,
                       in.getmin() - left_pad * in.getbinwidth(),
                       in.getmax() + right_pad * in.getbinwidth());
}

/**
 * Called before compressing the histogram by a factor 2 and
 * converting it from anti-zero-centric to zero-centric. This is more
 * or less the opposite of finalTweakFactor(). The factor is implied
 * to be 2. The number of negative and positive bins should now both
 * be odd, i.e. *not* a multiple of factor. That will convert the
 * histogram to zero-centric.
 */
HistogramBase /*static*/
ExpandableBuilder::finalTweakMakeZC(const HistogramBase& in)
{
  const size_type left_bins  = (size_type)std::max(0.0, std::ceil(-(in.getmin() / in.getbinwidth())));
  const size_type right_bins = in.getsize() - left_bins;
  const size_type left_pad   = (left_bins%2)  == 0 ? 1 : 0;
  const size_type right_pad  = (right_bins%2) == 0 ? 1 : 0;

  //std::cout << "zc factor 2"
  //          << " left "  << left_bins << "+" << left_pad
  //          << " right " << right_bins << "+" << right_pad
  //          << "\n";

  return HistogramBase(in.getsize() + left_pad + right_pad,
                       in.getmin() - left_pad * in.getbinwidth(),
                       in.getmax() + right_pad * in.getbinwidth());
}

/**
 * Second step in preparing the histogram for use: Reduce the number of bins
 * by combining buckets, leaving at most new_nbins. Typically, new_nbins
 * will be twice the bin count of the desired output histogram size.
 * Input and output are both anti-zero-centric.
 * The input bin width is constrained to 2^n, while the output bin width
 * is only constrained to m*2^n. m and n are integers.
 */
HistogramBase /*static*/
ExpandableBuilder::finalSqueeze(const HistogramBase& in, size_type want_nbins, size_type extra)
{
  const size_type factor = ((in.getsize()+want_nbins-1) / want_nbins) + extra;
  const HistogramBase tweaked = finalTweakFactor(in, factor);
  const size_type new_nbins = (tweaked.getsize() + factor - 1) / factor;

  //std::cout << "Tweaked:  " << fmt(tweaked) << "\n";

  // The squeezed histogram should have the same values for extreme
  // min and max range. The min and max center value will differ.
  const double old_width =
    (tweaked.getmax() - tweaked.getmin()) / (tweaked.getsize() - 1);
  const double new_width = old_width * factor;
  const double old_edge_min = tweaked.getmin() - (old_width / 2);
  const double old_edge_max = tweaked.getmax() + (old_width / 2);
  const double new_center_min = old_edge_min + (new_width / 2);
  const double new_center_max = old_edge_max - (new_width / 2);
  return HistogramBase(new_nbins, new_center_min, new_center_max);
}

/**
 * Convert a histogram that is known to be anti-zero-centric to a
 * zero-centric histogram that is half the size. Or one bin more.
 */
HistogramBase /*static*/
ExpandableBuilder::finalAntiToZeroCentric(const HistogramBase& in)
{
  const HistogramBase tweak = finalTweakMakeZC(in);

  // TODO-WIP-BrickedAPI: Avoid copy/paste.
  // Maybe have a constructor (nbins, histogram-to-take-range-from)
  // The squeezed histogram should have the same values for extreme
  // min and max range. The min and max center value will differ.
  const double old_width =
    (tweak.getmax() - tweak.getmin()) / (tweak.getsize() - 1);
  const double new_width = old_width * 2;
  const double old_edge_min = tweak.getmin() - (old_width / 2);
  const double old_edge_max = tweak.getmax() + (old_width / 2);
  const double new_center_min = old_edge_min + (new_width / 2);
  const double new_center_max = old_edge_max - (new_width / 2);
  return HistogramBase(tweak.getsize() / 2, new_center_min, new_center_max);
}

/**
 * Add empty bins to make the histogram size exactly as requested.
 * Works both for anti-zero-centric and zero-centric. For cosmetic
 * reasons the same number of bins are added on each end.
 */
HistogramBase /*static*/
ExpandableBuilder::finalInflate(const HistogramBase& in, size_type want_nbins)
{
  if (want_nbins > in.getsize()) {
    const size_type padding = want_nbins - in.getsize();
    const double    binwidth = (in.getmax() - in.getmin()) / (in.getsize() - 1);
    const size_type padding_left = padding / 2;
    const size_type padding_right = (padding + 1) / 2;
    return HistogramBase(want_nbins,
                         in.getmin() - (padding_left * binwidth),
                         in.getmax() + (padding_right * binwidth));
  }
  else {
    return in;
  }
}

/**
 * Convert the internal, symmetric anti-zero-centric histogram to a
 * zero-centric histogram with fewer bins.
 */
HistogramData
ExpandableBuilder::finalize(size_type want_nbins) const
{
  if (this->histo_.isfixed()) {
    if (want_nbins == this->histo_.getsize())
      return this->histo_;
    else
      throw std::runtime_error("Cannot resize a fixed histogram");
  }

  if (this->histo_.getsize() < 2*want_nbins || (this->histo_.getsize()%2) != 0)
    throw std::runtime_error("Internal error: Bad histogram range");

  // The finalXXX() methods only compute the required size and range,
  // they don't copy any bin counts.

  HistogramBase h(2, -1, +1); // Bogus range gets fixed immediately.
  for (int extra = 0; extra < 3; ++extra) {
    h = finalTrimEnds(this->histo_);
    //std::cout << "Trimmed:  " << fmt(h) << "\n";
    h = finalSqueeze(h, 2*want_nbins, extra);
    //std::cout << "Squeezed: " << fmt(h) << "\n";
    if (h.getsize() > 2*want_nbins) {
      // Usually we get larger bin width (and smaller size) than requested,
      // but protect against rounding errors etc. Try again if needed.
      //std::cout << "Overshoot! Asked finalSqueeze() for "
      //          << 2*want_nbins
      //          << " got " << h.getsize()
      //          << ". Try again.\n";
      continue;
    }
    h = finalAntiToZeroCentric(h);         // Now it is zero-centric.
    if (h.getsize() > want_nbins) {
      //std::cout << "Overshoot! Asked finalAntiToZeroCentric() for "
      //          << 2*want_nbins
      //          << " got " << h.getsize()
      //          << ". Try again.\n";
      continue;
    }
    //std::cout << "Zcentric: " << fmt(h) << "\n" << std::flush;
    h = finalInflate(h, want_nbins);
    //std::cout << "Inflated: " << fmt(h) << "\n";
    break;
  }

  // Create a fixed-size histogram to replace our temporary one.
  HistogramData histogram(h.getsize(), h.getmin(), h.getmax());
  histogram += this->histo_;
  //this->histo_ = histogram;
  return histogram;
}

} // end namespace
