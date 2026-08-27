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
#pragma once

#include "../declspec.h"
#include "roundandclip.h"

#include <memory>
#include <math.h>
#include <cmath>
#include <string>
#include <cstdint>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file histogramdata.h
 * \brief provides class \ref InternalZGY::HistogramData.
 */

/**
 * This could have been a base class of HistogramData, but currently isn't.
 * Only ExpandableBuilder needs to deal with histograms without actual data.
 */
class OPENZGY_TEST_API HistogramBase
{
public:
  typedef std::int32_t size_type;

  HistogramBase(size_type nbins, double min, double max)
    : nbins_(nbins)
    , min_(min)
    , max_(max)
  {
  }

  size_type getsize() const { return nbins_; } /**< Number of bins. */
  double     getmin() const { return min_; }   /**< Center value of first bin */
  double     getmax() const { return max_; }   /**< Center value of last bin */
  double     getbinwidth() const { return (max_-min_) / (nbins_-1); }

protected:
  size_type  nbins_;
  double     min_;
  double     max_;
};

/**
 * \brief A histogram for a data set.
 *
 * \internal
 * KEEP THIS TEXT IN SYNC WITH OpenZGY::SampleHistogram
 * \endinternal
 *
 * The histogram is described by the fixed total number of bins, the
 * center value of the samples in the first bin, and the center value
 * of the samples in the last bin.
 *
 * The width of each bin is given by (max - min) / (nbins - 1).
 * Bin 0 holds samples with values min_ +/-binwidth/2.
 *
 * This means that the total range of samples that can be represented in the
 * histogram is actually min-binwidth/2 to max+binwidth/2, which is slightly
 * larger than just min..max _unless_ each of the bins can only hold a single
 * value (e.g. data converted from 8-bit storage, in 256 bins).
 *
 * This definition has some subtle effects, best illustrated by a few examples.
 *
 * If the source data is ui8_t, the logical range is 0..255. Assume nbins=256.
 * This gives binwidth=1, and the true range of the histogram -0.5..+255.5.
 * But since the input is integral, the actual range is just 0..255 inclusive.
 * Try to fill the histogram with evenly distrubuted random data and you end
 * up with each bin having roughly the same number of elements.
 *
 * Now consider ui16_t, range 0..65535 and nbins is still 256. This gives
 * binwidth=257, not 256. The true range of the histogram is -128.5..+65663.5.
 * Try to fill the histogram with evenly distrubuted random data and you end
 * up with the first and the last bin having approximately half as many
 * elements as all the others. This is not really a problem, but may seem
 * a bit surprising.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_TEST_API HistogramData
{
public:
  typedef std::int32_t size_type;
  typedef std::int64_t count_type;

  HistogramData(size_type nbins, double min, double max);
  HistogramData(size_type nbins, double min, double max, bool fixed);
  HistogramData(const count_type* bins, int nbins, double min, double max);
  HistogramData(const HistogramData& other);
  HistogramData& operator=(const HistogramData& other);
  HistogramData& operator+=(const HistogramData& other);
  HistogramData& operator-=(const HistogramData& other);
  HistogramData& operator*=(count_type factor);
  bool operator==(const HistogramData& other) const = delete; // YAGNI
  bool operator!=(const HistogramData& other) const = delete; // YAGNI

  // Queries
  size_type getsize()              const { return nbins_; }      /**< Get the number of bins in the histogram. */
  const count_type* getbins()      const { return bins_.get(); } /**< Get the actual bins. */
  double     getmin()              const { return min_; }        /**< Center value of data in first bin */
  double     getmax()              const { return max_; }        /**< Center value of data in last bin */
  bool       isfixed()             const { return fixed_; }
  count_type get(double value)     const;                        /**< Get sample count for this value */
  count_type getcount()            const;                        /**< Get sample count for all values */
  std::string toString(bool verbose = false) const;
  void       toGraph(std::ostream& os) const;
  double     getbinwidth()         const { return (max_-min_) / (nbins_-1); }

  // Operations
  void scale(double oldmin, double oldmax, double newmin, double newmax);

  // Internal, to be called from HistogramBuilder only.
public:
  void clear();
  static void getLinearTransform(double *offset, double *scale, double oldmin, double oldmax, double newmin, double newmax);
  static void combineLinearTransform(double *A, double *B, double A1, double B1, double A2, double B2);
  void calculateConversionFactors(double* A, double* B) const; // Called from get(), AddBins(), TryAdd().
  count_type* dangerousGetBins() { return bins_.get(); }  /**< Will corrupt data if binno outside range */
  void setNewRange(double min, double max) {min_ = min; max_ = max;}

private:
  void addBins(const count_type* bins, int nbins, double min, double max, bool add); // Called from +=, -=.
  inline void addOne(double value, count_type factor, double A, double B); // Called from AddBins() only.

protected:
  double            min_;       /**< Lower end of histogram value range */
  double            max_;       /**< Upper end of histogram value range */
  size_type         nbins_;     /**< Number of bins, guaranteed to be >=2 */
  std::unique_ptr<count_type[]> bins_;
  bool              fixed_;
};

} // end namespace
