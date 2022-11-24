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
#pragma once

#include "../declspec.h"
#include "statisticdata.h" // needed due to data members.
#include "histogramdata.h" // needed due to data members.

#include <cmath> // needed for std::isfinite

#include <limits>
#include <iterator> // needed for std::iterator_traits
#include <assert.h>
#ifdef LINUX
#include <stdint.h> // for intptr_t
#endif

#define XXX_DISABLE_NAN_CHECK 0

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file histogrambuilder.h
 * \brief Collect statistics and histogram for bulk data.
 */

/**
 * \brief Collect statistics and histogram for bulk data.
 *
 * Note the following caveat with histogram data. The histogram is described by
 * number of bins and a pair of min/max values. It is not obvious whether those
 * values describe the center value of the first and last bin, or whether they
 * are the open ended range representing the whole histogram. I.e. the first value
 * (inclusive) of the first bin and the last value (exclusive) of the last bin.
 * HistogramBuilder uses the former definition.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are meant to be fully populated prior to being
 * made available to others.
 */
class OPENZGY_TEST_API HistogramBuilder
{
public:

  typedef int                        size_type;
  typedef StatisticData::count_type  count_type;

  HistogramBuilder(size_type _nbins, double _min, double _max);
  HistogramBuilder(const count_type* bins, int nbins, double min, double max, count_type  scnt, double ssum, double sssq, double smin, double smax);

  HistogramBuilder& operator+=(const HistogramBuilder& other);
  HistogramBuilder& operator-=(const HistogramBuilder& other);
  HistogramBuilder& operator*=(count_type factor);
  bool operator==(const HistogramBuilder& other) const;
  bool operator!=(const HistogramBuilder& other) const;

  // Queries
  const StatisticData& getstats() const { return stats_; } /**< StatisticBuilder on added data */
  const HistogramData& gethisto() const { return histo_; } /**< Histogram from added data */
  StatisticData        gethiststats() const;               /**< Statistics calculated from bins */

  // Data collection
  template <typename It> void add(It begin, It end);

  // Other operations
  void scale(double oldmin, double oldmax, double newmin, double newmax);

private:
  void clear();
  template <typename It> void tryAdd(It begin, It end, StatisticData *localstats);
  template <typename It> void fastAddSignedChar256(It begin, It end);

private:
  StatisticData     stats_;     /**< Statistics from added data */
  HistogramData     histo_;     /**< Histogram from added data */
};

/**
 * Short cut version of Add, only works for "char" type, no scaling,
 * histogram range fixed at entire range of the "char" type, number of bins==256.
 */
template <typename It>
void HistogramBuilder::fastAddSignedChar256(It begin, It end)
{
  assert(this->gethisto().getsize() == 256 &&
         std::numeric_limits<typename std::iterator_traits<It>::value_type>::is_integer &&
         std::numeric_limits<typename std::iterator_traits<It>::value_type>::min()+128 == 0 &&
         std::numeric_limits<typename std::iterator_traits<It>::value_type>::max() == +127 &&
         std::numeric_limits<typename std::iterator_traits<It>::value_type>::min() == this->gethisto().getmin() &&
         std::numeric_limits<typename std::iterator_traits<It>::value_type>::max() == this->gethisto().getmax());

  // For char input, there is no risk of overflow and no scaling - just an offset.
  for (It i = begin; i != end; ++i)
    histo_.addToBinNumberUNSAFE(static_cast<char>(*i) + 128);

  // stats_ needs to be re-calculated. Assuming that the typical call to add()
  // contains 256k samples, it is more efficient to do it this way.
  // For char data (unlike wider types), we will get a precise result.
  // It would have been even faster to defer this until the stats were
  // actually requested, but this way is safer and hopefully fast enough.
  // Note: if I ever rewrite the calling code to parallelize the histogram
  // generation, the chunks will likely be smaller and in that case something
  // will need to be done to avoid calculating stats too often.
  stats_ = StatisticData(gethisto().getbins(), gethisto().getsize(), gethisto().getmin(), gethisto().getmax(), true);
}

/**
 * Add samples from an iterator to a histogram.
 * If the caller has already calculated statistics from the iterator,
 * those can be passed it so this function doesn't need to do it again.
 */
template <typename It>
void HistogramBuilder::add(It begin, It end)
{
  // verify that we need to do anything
  if (begin == end) {
    return;
  }

  // Short cuts for simple cases that can be handled more efficiently.
  if (std::numeric_limits<typename std::iterator_traits<It>::value_type>::is_integer &&
      std::numeric_limits<typename std::iterator_traits<It>::value_type>::min() == this->gethisto().getmin() &&
      std::numeric_limits<typename std::iterator_traits<It>::value_type>::max() == this->gethisto().getmax())
  {
    if (this->gethisto().getsize() == 256 &&
        std::numeric_limits<typename std::iterator_traits<It>::value_type>::min()+128 == 0 &&
        std::numeric_limits<typename std::iterator_traits<It>::value_type>::max() == +127)
    {
      // This can be dramatically (10x) faster, since we don't need to increment statistics for each sample.
      fastAddSignedChar256(begin, end);
      return;
    }
  }

  // Calculate statistics and histogram information in a single pass.
  StatisticData tempstats;
  tryAdd<It>(begin, end, &tempstats);
  stats_ += tempstats;
}

/**
 * Internal version of Add(), handles both the no-convert case and
 * the case where the caller passes a range to convert the samples to.
 */
template <typename It>
void HistogramBuilder::tryAdd(It begin, It end, StatisticData *localstats)
{
  double combo_offset = 0, combo_scale = 1; // from storage value to bin number.

  // conversion from coding range to bin number.
  histo_.calculateConversionFactors(&combo_offset, &combo_scale);

  // Performance note: If I add +1.5 to offset here, then instead of the
  // RoundD2I later (round to nearest 0-based bin index) I would instead be doing
  // a simple assign from float to int, rounding towards zero, to get as a result
  // a 1-based bin index. Then subtract 1 to get the 0-based bin index.
  // This is more efficient on x64 and when /arch:sse2 is given, but less
  // efficient in plain x86 mode. The reason I wouldn't just add 0.5 is that
  // it would not handle the range [-1, -0.5] correctly as they would map to
  // bin 0 and not an invalid bin number.

  // More on the topic of rounding: If the value is so huge that it won't fit
  // in a 32-bit integer then the RoundD2I result is technically undefined.
  // I am assuming that the result won't fall in the range of valid bin numbers,
  // and that *appears* to be true on the IA32 architecture. But the assumption
  // is fragile. "n" must not be declared as unsigned; we need to cast it.
  // See HardwareOptimized.h for a longer discussion.

  // To avoid having to iterate twice over the input, both the statistics
  // and the histogram information will be updated in the same loop.
  // Furthermore, the statistics will be collected in temporary variables
  // so the compiler has a chance of placing them in registers.
  // Yes, this is a bit messy because it muddles the separation of
  // Histogram and Statistics, and it duplicates logic from the latter.
  std::int64_t cnt = 0;
  std::int64_t inf = 0;
  double sum = 0;
  double ssq = 0;
  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();

  for (It i = begin; i != end; ++i) {
    register double value = *i;

    // The test for std::isfinite is not needed if *i is an integral type.
    if (XXX_DISABLE_NAN_CHECK || std::numeric_limits<typename std::iterator_traits<It>::value_type>::is_integer || std::isfinite(value)) {

      int n = RoundD2I(combo_offset + combo_scale * value);
      if (n >= 0 && n < gethisto().getsize())
        histo_.addToBinNumberUNSAFE(n);

      // Always accumulate statistics. localstats is almost always != 0.
      sum += value;
      ssq += value*value;
      if (min > value)
        min = value;
      if (max < value)
        max = value;
      ++cnt;
    } else {
      ++inf;
    }
  }

  if (localstats != NULL) {
    *localstats = StatisticData(cnt, inf, sum, ssq, min, max);
  }
}

} // end namespace
