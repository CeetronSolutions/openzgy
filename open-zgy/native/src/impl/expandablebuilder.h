// Copyright 2017-2023, Schlumberger
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

/*

/home/paal/git/Salmon-master/src/Salmon/Zgy/UtilityLib/HistogramBuilder.h
/home/paal/git/Salmon-master/src/Salmon/Zgy/UtilityLib/HistogramBuilder.cpp

*/

#pragma once

#include "../declspec.h"
#include "histogrambuilder.h"
#include "statisticdata.h"
#include "histogramdata.h"
#include "roundandclip.h"

#include <cmath>
#include <limits>
#include <iterator>
#include <assert.h>
#include <stdexcept>
#include <atomic>
#ifdef LINUX
#include <stdint.h> // for intptr_t
#endif

#define XXX_DISABLE_NAN_CHECK 0

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file expandablebuilder.h
 * \brief Collect statistics and histogram for bulk data.
 */

/**
 * \brief Collect statistics and histogram for bulk data.
 * The histogram doesn't need to be known a priori.
 *
 * \details
 *
 * Histograms may be created as fixed (min/max will never change), or
 * dynamic (the histogram may be extended when samples are seen
 * outside the current range).
 *
 * If fixed is true, min and max must be specified and give the sample
 * range. The result of attempting to add samples outside this range
 * is not yet decided; the choice might be significant in relation to
 * numerical instability. And for that issue, the choice also affects
 * dynamic histogram. TODO-WIP-BrickedAPI: revisit the current choice.
 *   - (current) completely ignored.
 *   - clipped to highest or lowest bin.
 *   - placed in a "below" or "above" bins.
 *
 * If fixed is false, min and max cannot be specified.
 */
class OPENZGY_TEST_API ExpandableBuilder : public HistogramBuilder
{
public:
  // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM does not allow explicitly
  // specifying the histogram range. There are still unit tests
  // doing this. Eventually this constructor should be removed.
  ExpandableBuilder(size_type nbins, double min, double max, bool fixed);
  explicit ExpandableBuilder(size_type nbins);
  ExpandableBuilder(size_type nbins, double min, double max);
  ExpandableBuilder(const HistogramData& histo, const StatisticData& stats, bool copy);
  ExpandableBuilder(const ExpandableBuilder&) = default;
  ExpandableBuilder(ExpandableBuilder&&) = default;
  ExpandableBuilder& operator+=(const HistogramBuilder& other);
  ExpandableBuilder& operator-=(const HistogramBuilder& other);
  ExpandableBuilder& operator*=(count_type factor);

  // Data collection. Shadows the method in the base class.
  template <typename It> void add(It begin, It end);

  // Convert to a zero-centric histogram with fewer bins.
  HistogramData finalize(size_type want_nbins) const;

  // Static version of operator+= and operator-=
  static void staticAddOrSub(
       StatisticData* stats,
       HistogramData* histo,
       const StatisticData& other_stats,
       const HistogramData& other_histo,
       bool add);

private:
  bool widenRangeIsNeeded(double new_min, double new_max) const;
  static bool widenRangeIsNeeded(
       const StatisticData* stats,
       const HistogramData* histo,
       double new_min, double new_max);
  void widenRange(double new_min, double new_max);
  static void widenRange(
       StatisticData* stats,
       HistogramData* histo,
       double new_min, double new_max);

  static double getsmallestbinwidth();
  static std::pair<double,double> antiZcRangeFromBinWidth(
       size_type N, double binwidth, bool center);
  static std::pair<double,double> smallestAntiZcRange(size_type N, bool center);
  static double antiZcBinWidthFromRange(size_type nbins, double min, double max);
  static HistogramBase finalTrimEnds(
       const HistogramData& in);
  static HistogramBase finalTweakFactor(
       const HistogramBase& in, size_type factor);
  static HistogramBase finalTweakMakeZC(const HistogramBase& in);
  static HistogramBase finalSqueeze(
       const HistogramBase& in, size_type want_nbins, size_type extra);
  static HistogramBase finalAntiToZeroCentric(
       const HistogramBase& in);
  static HistogramBase finalInflate(
       const HistogramBase& in, size_type want_nbins);
};

/**
 * Add samples from an iterator to the statistics and histogram.
 */
template <typename It>
void ExpandableBuilder::add(It begin, It end)
{
  // verify that we need to do anything
  if (begin == end) {
    return;
  }

  // Short cuts for simple cases that can be handled more efficiently.
  if (isSignedChar256<It>()) {
    // TODO-WIP-BrickedAPI: The overhead here might skew the total time.
    //auto& timer = HistoBuilderAdHocTimers::instance().histogram_addbytes;
    //SimpleTimerEx tt(timer);
    fastAddSignedChar256<It, 1>(begin, end);
    //timer.addBytesRead(end - begin);
    return;
  }

  // Calculate statistics and histogram information in a single pass.
  // It might be possible to skip tempstats and just pass in &stats_
  // instead. (tryAdd would needs to use += to add stats instead of
  // replacing them. And the calls to widenRange() would operate on
  // total stats seen.) Since StatisticData::operator+= is so cheap,
  // it probably doesn't make much difference.
  StatisticData tempstats;
#if 0
  tryAdd<It, 1>(begin, end, &tempstats);
  //std::cout << "tryAdd #1 got " + tempstats.toString() + " HISTO " + histo_.toString() + "\n" << std::flush;
#else
  if ((end - begin) < 32768)
  {
    // TODO-WIP-BrickedAPI: The overhead here might skew the total time.
    auto& timer = HistoBuilderAdHocTimers::instance().expandable_addsmall;
    SimpleTimerEx tt(timer);
    tryAdd<It, 1>(begin, end, &tempstats);
    timer.addBytesRead(end - begin);
  }
  else {
    auto& timer = HistoBuilderAdHocTimers::instance().expandable_addbig;
    SimpleTimerEx tt(timer);
    tryAdd<It, 1>(begin, end, &tempstats);
    timer.addBytesRead(end - begin);
  }
#endif
  stats_ += tempstats;

  // Specific to HAVE_EXPANDABLE_BUILDER. Since we didn't check the
  // range beforehand, some or all our samples might have ended up
  // outside the bins. When this happens, we need needs to back out
  // the histogram counts that were just added, fix the histogram
  // range, and try again. The statistics are already ok, so don't
  // touch those.
  //
  // TODO-WIP-BrickedAPI: If caller is multi-threaded then "this" is
  // probably a temp instance. In which case it might be ok to leave
  // the retry logic to the caller. If add() is called with too little
  // data at a time, this might improve performance

  if (tempstats.getcnt() != 0) { // not all NaN
    if (widenRangeIsNeeded(tempstats.getmin(), tempstats.getmax())) {
      SimpleTimerEx tt(HistoBuilderAdHocTimers::instance().expandable_expand);
#if 1 // 1 for normal, 0 for debugging.
      tryAdd<It, -1>(begin, end, nullptr);
      widenRange(tempstats.getmin(), tempstats.getmax());
      tryAdd<It, 1>(begin, end, nullptr);
#else
      StatisticData t2, t3;
      tryAdd<It, -1>(begin, end, &t2);
      std::cout << "tryAdd #2 sub " + t2.toString() + " HISTO " + histo_.toString() + "\n" << std::flush;
      widenRange(tempstats.getmin(), tempstats.getmax());
      tryAdd<It, 1>(begin, end, &t3);
      std::cout << "tryAdd #3 add " + t2.toString() + " HISTO " + histo_.toString() + "\n" << std::flush;
#endif
    }
  }
}

} // end namespace
