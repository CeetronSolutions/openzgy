// Copyright 2017-2021, Schlumberger
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

//Adapted from: Zgy/ArrayBasic/ArrayTile

#include "lodalgo.h"
#include "lodsampling.h"
#include "roundandclip.h"
#include "databuffer.h"
#include "fancy_timers.h"
#include "environment.h"
#include "mtguard.h"
#include "exception.h"

#include <memory.h>
#include <math.h>
#include <assert.h>
#include <limits>
#include <typeinfo>
#include <memory>
#include <atomic>
#include <omp.h>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \brief Weighted arithmetic average of 8 neighboring samples.
 *
 * Each sample is weighted agaist how common this value is in the survey
 * as a whole. Rare values (likely very high of very low) get higher
 * priority. This prevents tiles from looking more and more washed out
 * the more they get decimated.
 *
 * TODO-Performance: Profile this function and try to speed it up.
 * It is by far the most expensize of the decimators, ~10x the cost
 * of LowPass. By default it is used for lod 2+ only. So it gets called
 * with just 1/7 of the data passed to LowPass. In sum the code still
 * spends more CPU cycles on this than on LowPass.
 *
 * \details Thread safety:
 * Modification may lead to a data race. This should not be an issue.
 * Instances are meant to be short lived and used only from one place.
 */
template<typename T>
class LodWeightedAverage
{
  // Information passed to constructor
  const std::int64_t* const hist_;
  const int bins_;
  const double minHist_;
  const double maxHist_;
  // pre-calculated in constructor
  double hsum_;
  double factor_;
  double offset_;

public:
  LodWeightedAverage(const std::int64_t *hist, int bins, double minHist, double maxHist) :
    hist_(hist), bins_(bins), minHist_(minHist), maxHist_(maxHist)
  {
    // find total number of samples to base weighting on
    hsum_ = 0.0;
    for (int dsti = 0; dsti < bins_; dsti++){
      hsum_ += static_cast<double>(hist_[dsti]);
    }

    // Avoid divide by zero in weights below
    if (hsum_ < 1.0){
      hsum_ = 1.0;
    }

    // Linear transform from sample value to bin number.
    // Note that result is supposed to be rounded to nearest integer.
    factor_ = std::isfinite(maxHist) && std::isfinite(minHist) && maxHist>minHist && bins>1 ? ((bins-1)/(maxHist-minHist)) : 0.0;
    offset_ = -minHist * factor_;
  }

  static const LodAlgorithm algorithm = LodAlgorithm::WeightedAverage;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    T data[8] = { s0, s1, s2, s3, s4, s5, s6, s7};

    // calculate weighted sum for 2x2x2 part of src brick
    double weight(0.0), sum(0.0);
    for (int ii=0; ii<8; ++ii) {

      // get value from src brick
      const T val = data[ii];

      if (!IsFiniteT(val))
          continue;

      // Get histogram bin number for this sample, then look up this
      // sample's frequency. See HistogramData::get().
      double histv;
      int hVal= RoundD2I(val * factor_ + offset_);
      if (hVal >= 0 && hVal < bins_) {
        histv = static_cast<double>(hist_[hVal]);
      } else {
        histv = 0.0;  // outside histogram, frequency definitely 0
      }
      // If frequency is reported as 0, this is an inconsistency.
      // Since we know it occurs at least one time.
      if (histv<1.0){ histv = 1.0; }

      // update weighted sum
      sum += val * hsum_ / histv;

      // update weight
      weight += hsum_ / histv;
    }

    double wval;
    if (weight != 0) // use weighted value, at least one finite value found.
      wval = sum / weight;
    else // all 8 values were infinite or NaN, just pick one of them.
      wval = s0;

    if (wval > maxHist_){ wval = maxHist_; }
    if (wval < minHist_){ wval = minHist_; }

    // assign weighted value to destination brick
    return static_cast<T>(wval);
  }
};

/**
 * \brief Arithmetic average of 8 neighboring integer samples.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodAverage
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Average;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    double sum =
      static_cast<double>(s0) +
      static_cast<double>(s1) +
      static_cast<double>(s2) +
      static_cast<double>(s3) +
      static_cast<double>(s4) +
      static_cast<double>(s5) +
      static_cast<double>(s6) +
      static_cast<double>(s7);
    return static_cast<T>(sum/8);
  }
};

/**
 * \brief Arithmetic average of 8 neighboring float samples, ignoring NaN.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<>
class LodAverage<float>
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Average;
  float operator()(float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7)
  {
    double sum = 0;
    int count = 0;
    if (std::isfinite(s0)) { sum += s0; ++count; }
    if (std::isfinite(s1)) { sum += s1; ++count; }
    if (std::isfinite(s2)) { sum += s2; ++count; }
    if (std::isfinite(s3)) { sum += s3; ++count; }
    if (std::isfinite(s4)) { sum += s4; ++count; }
    if (std::isfinite(s5)) { sum += s5; ++count; }
    if (std::isfinite(s6)) { sum += s6; ++count; }
    if (std::isfinite(s7)) { sum += s7; ++count; }
    if (count == 0)
      return s0; // all infinite, so this we return.
    return static_cast<float>(sum/count);
  }
};

/**
 * \brief Arithmetic average, excluding zero, of integer data.
 *
 * \details This is not useful for seismic data. It might be handy
 * if samples are by nature integer values, not just floating point
 * numbers (such as seismic amplitude) crammed into a too-small integer.
 *
 * Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodAverageNon0
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::AverageNon0;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    double sum =
      static_cast<double>(s0)+
      static_cast<double>(s1)+
      static_cast<double>(s2)+
      static_cast<double>(s3)+
      static_cast<double>(s4)+
      static_cast<double>(s5)+
      static_cast<double>(s6)+
      static_cast<double>(s7);
    int count = 0;
    if (s0 != 0) ++count;
    if (s1 != 0) ++count;
    if (s2 != 0) ++count;
    if (s3 != 0) ++count;
    if (s4 != 0) ++count;
    if (s5 != 0) ++count;
    if (s6 != 0) ++count;
    if (s7 != 0) ++count;
    if (count == 0)
      return 0;
    else
      return static_cast<T>(sum / count);
  }
};

/**
 * \brief Arithmetic average, excluding NaN and zero, of float data.
 *
 * \details CAVEAT: This is not useful for seismic data. It might be handy
 * if samples are by nature integer values, not just floating point
 * numbers (such as seismic amplitude) crammed into a too-small integer.
 * This means that running the algorithm on float data is probably a mistake.
 *
 * Thread safety: Safe. Instance holds no data.
 */
template<>
class LodAverageNon0<float>
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::AverageNon0;
  float operator()(float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7)
  {
    double sum = 0;
    int count = 0;
    if (std::isfinite(s0) && s0 != 0) { sum += s0; ++count; }
    if (std::isfinite(s1) && s1 != 0) { sum += s1; ++count; }
    if (std::isfinite(s2) && s2 != 0) { sum += s2; ++count; }
    if (std::isfinite(s3) && s3 != 0) { sum += s3; ++count; }
    if (std::isfinite(s4) && s4 != 0) { sum += s4; ++count; }
    if (std::isfinite(s5) && s5 != 0) { sum += s5; ++count; }
    if (std::isfinite(s6) && s6 != 0) { sum += s6; ++count; }
    if (std::isfinite(s7) && s7 != 0) { sum += s7; ++count; }
    if (count == 0) {
      // All infinite or zero. If at least one zero return that, else return any of the infinite values.
      if (s0 == 0 || s1 == 0 || s2 == 0 || s3 == 0 || s4 == 0 || s5 == 0 || s6 == 0 || s7 == 0)
        return 0;
      else
        return s0;
    }
    return static_cast<float>(sum / count);
  }
};

/**
 * brief Median of 8 neighboring integer samples.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodMedian
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Median;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    T data[8] = { s0, s1, s2, s3, s4, s5, s6, s7};

    // Bubble sort the array of 8 source numbers.
    // Yes, Bubble sort. The array is short enough
    // that this should not be a problem.
    // I might not even bother bailing out as soon
    // as the array is proven to be sorted.
    for (int ii=7; ii>0; --ii)
      for (int jj=0; jj<ii; ++jj)
        if (data[jj]>data[jj+1])
          {
            T tmp = data[jj+1];
            data[jj+1] = data[jj];
            data[jj] = tmp;
          }
    return static_cast<T>(((float)data[3] + (float)data[4])/2.0);
  }
};

/**
 * \brief Median of 8 neighboring float samples, ignoring NaN.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<>
class LodMedian<float>
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Median;
  float operator()(float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7)
  {
    float data[8] = { s0, s1, s2, s3, s4, s5, s6, s7};
    int count = 0;
    if (std::isfinite(s0)) { data[count++] = s0; };
    if (std::isfinite(s1)) { data[count++] = s1; };
    if (std::isfinite(s2)) { data[count++] = s2; };
    if (std::isfinite(s3)) { data[count++] = s3; };
    if (std::isfinite(s4)) { data[count++] = s4; };
    if (std::isfinite(s5)) { data[count++] = s5; };
    if (std::isfinite(s6)) { data[count++] = s6; };
    if (std::isfinite(s7)) { data[count++] = s7; };

    for (int ii=count-1; ii>0; --ii)
      for (int jj=0; jj<ii; ++jj)
        if (data[jj]>data[jj+1])
          {
            float tmp = data[jj+1];
            data[jj+1] = data[jj];
            data[jj] = tmp;
          }
    if (count == 0)
      return s0; // All NaN
    else if (count%2 == 1)
      return data[count/2]; // odd number of bins - return the center one
    else
      return (data[count/2-1] + data[count/2])/2;
  }
};

/**
 * \brief Minimum of 8 neighboring samples, excluding NaN
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodMinimum
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Minimum;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    T result = s0;
    // Invert the test, so if result is NaN the test will succeed
    // and we will favor returning a result that is not NaN.
    if (!(result <= s1)) result = s1;
    if (!(result <= s2)) result = s2;
    if (!(result <= s3)) result = s3;
    if (!(result <= s4)) result = s4;
    if (!(result <= s5)) result = s5;
    if (!(result <= s6)) result = s6;
    if (!(result <= s7)) result = s7;
    return result;
  }
};

/**
 * \brief Maximum of 8 neighboring samples, excluding NaN
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodMaximum
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Maximum;
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    T result = s0;
    if (!(result >= s1)) result = s1;
    if (!(result >= s2)) result = s2;
    if (!(result >= s3)) result = s3;
    if (!(result >= s4)) result = s4;
    if (!(result >= s5)) result = s5;
    if (!(result >= s6)) result = s6;
    if (!(result >= s7)) result = s7;
    return result;
  }
};

/**
 * \brief Set the low resolution data to all zeros.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodAllZero
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::AllZero;
  T operator()(T, T, T, T, T, T, T, T)
  {
    return static_cast<T>(0);
  }
};

/**
 * \brief Use just one of the 8 input values.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodDecimate
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::Decimate;
  T operator()(T s0, T, T, T, T, T, T, T)
  {
    return s0;
  }
};

/**
 * \brief Use just one of the 8 integral input values.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<typename T>
class LodDecimateSkipNaN
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::DecimateSkipNaN;
  T operator()(T s0, T, T, T, T, T, T, T)
  {
    return s0;
  }
};

/**
 * \brief Use just one of the 8 float input values, looking for one not NaN.
 * \details Thread safety: Safe. Instance holds no data.
 */
template<>
class LodDecimateSkipNaN<float>
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::DecimateSkipNaN;
  float operator()(float s0, float s1, float s2, float s3, float s4, float s5, float s6, float s7)
  {
    if (std::isfinite(s0)) return s0;
    if (std::isfinite(s1)) return s1;
    if (std::isfinite(s2)) return s2;
    if (std::isfinite(s3)) return s3;
    if (std::isfinite(s4)) return s4;
    if (std::isfinite(s5)) return s5;
    if (std::isfinite(s6)) return s6;
    if (std::isfinite(s7)) return s7;
    return s0; // All are infinite, so this is what we return.
  }
};

#ifdef _MSC_VER
#pragma warning (push)
#pragma warning (disable: 4127) // "conditional expression is constant" deliberate in template code
#endif

/**
 * \brief Most frequent value.
 *
 * Return the value that occurs most frequently in the input.
 * If there is a tie, return the first one.
 * NaN or 0 are returned if and only if all inputs have this value.
 * If the input is a mix of 0 and NaN only, the result is 0.
 *
 * Since there are tests for equality, this algorithm probably only
 * makes sense for integral data.
 *
 * Thread safety: Safe. Instance holds no data.
 */
template<typename T, bool SkipNaN, bool SkipZero>
class LodMostFrequent
{
public:
  static const LodAlgorithm algorithm = LodAlgorithm::DecimateSkipNaN; // TODO-Low, this is wrong, and depends on template args.
  T operator()(T s0, T s1, T s2, T s3, T s4, T s5, T s6, T s7)
  {
    T in[8] = {s0, s1, s2, s3, s4, s5, s6, s7};
    signed char counts[8] = {0};
    for (int jj=1; jj<8; ++jj)    // how many inputs after the first one
      if (in[0]==in[jj])          // have the same value as the first?
        ++counts[0];              // end result will be in range 0..7
    if (counts[0] == 7)           // all 8 values are equal and not NaN
      return in[0];               // so this is a short cut
    for (int ii=1; ii<8; ++ii)    // start at 1, counts[0] done above
      for (int jj=ii+1;jj<8;++jj) // how many inputs above this one
        if (in[ii]==in[jj])       // have the same numerical value?
          ++counts[ii];           // NaN will always end up with count 0
    if (SkipNaN)                  // typically true if T is floating point
      for (int ii=0; ii<8; ++ii)  // caller use std::numeric_limits to check
        if (!(IsFiniteT(in[ii]))) // this checks both NaN and +/- inf
          counts[ii] = -2;        // strong bias against NaN in result
    if (SkipZero)                 // set depending on application's choice
      for (int ii=0; ii<8; ++ii)  // note: comparing against storage 0 not
        if (in[ii] == 0)          // converted 0, maybe we should change that?
          counts[ii] = -1;        // weaker bias against 0 than against NaN
    int ok = 0;                   // locate entry with the highest count
    for (int ii=1; ii<8; ++ii)    // note that all loops might be unrolled
      if (counts[ii]>counts[ok])  // could almost avoid building the in[] array
        ok=ii;                    // found a better result
    return in[ok];                // the one place where in[] should be array.
  }
};

#ifdef _MSC_VER
#pragma warning (pop) // restore "conditional expression is constant"
#endif

/**
 * Generic LOD calculation that can be used for all algorithms where
 * the resulting sample depends only on the 8 surrounding samples
 * in the source brick. I.e. currently everything except the lowpass
 * algorithm that needs a window of size 10. The actual algorithm to
 * be used is passed as a functor. If any calculation needs to be done
 * up front, that can be handled by the functor's constructor.
 *
 * TODO-Medium: Support for files containing just a single section or
 * skice. Possibly also with bricksize 1 in one of the dimensions
 *
 * Called from: GenLodImpl::_calculate() --> GenLodImpl::_decimate()
 *
 * There are several ways the caller can choose to handle any padding
 * area between the survey edge and the edge of the last brick.
 * The _caller_ chooses what it sends. Currently (1d), planned (1a).
 * There are several ways of handling the last sections and/or slice
 * if the survey size is odd. _This_ function chooses what to do.
 *
 * If you at this point are thinking "who cares?", note that there are
 * plans to allow a ZGY file to store a single slice or a single
 * section. In that case all the samples fall in the "odd size" group.
 * So to get decimation to work at all we will need one of the fancier
 * options. Or sidestep the issue by claiming that those files have no
 * low resolution data stored.
 *
 * 1) The input buffer has been shrunk to fit the actual data available.
 *    Both the source and target strides describe a c-contiguous buffer.
 *
 *    Handling of odd sizes:
 *
 *    1a) Use the available 4, 2, or 1 samples.
 *    1b) Use just one sample and the trivial Decimate algorithm.
 *    1c) Fill the edges with a default value to be specified by caller.
 *   *1d) Fill the edges with zero and hope nobody notices.
 *    1e) Don't touch the output samples at the edge. Caller should
 *        initialize the entire buffer to desired defaultvalue.
 *
 * 2) The input buffer includes the padding area and the source and
 *    target strides have been set up to describe a non-contiguous
 *    buffer layout. Size refers to the used part of the cube. The
 *    caller must initialize the entire buffer to the desired
 *    defaultvalue because won't know how large the buffer really is.
 *    At least not in the slowest dimension.
 *
 *    Handling of odd sizes: All options above would work, but c and d
 *    are rather pointless.
 *
 * 3) The buffers, sizes, and strides describe the entire brick which
 *    means this decimation function won't know what is real data and
 *    what is padding. The caller must ensure that the *input* buffer
 *    has sensible values in the padding area. Garbage values will
 *    cause the edges of the result to contain garbage if real size
 *    is odd. Also, processing the garbage sample values might cause
 *    arithmetic overflow, signalling NaN, and other mayhem.
 *    Also, a brick size of 1 will still need special handling.
 *    (Other than that case, sizes will always be even).
 *
 *    Handling of odd sizes: In this case it is the caller that decides,
 *    it depends on how the caller pads the input. It doesn't matter
 *    which of abcde is implemented here because we will never see any
 *    odd sizes. Except for bricksize 1; in that case it still matters
 *    what happens in this function. Bottom line: I didn't realize it
 *    before starting to write this (way too verbose) explanation.
 *    Alternative (3) is simply a very bad choice.
 */
template <typename T, typename F>
static void
createGenericLevelOfDetail(
    T* dst,
    const std::array<std::int64_t,3>&dsize,
    const std::array<std::int64_t,3>&dstride,
    const T* src,
    const std::array<std::int64_t,3>&ssize,
    const std::array<std::int64_t,3>&sstride,
    F function)
{
  //printf("@createGenericLevelOfDetail<%s,%s> Algo=%d\n",
  //        typeid(T).name(), typeid(F).name(), (int)function.algorithm);

  // The samples of the output where all 8 inputs are available.
  const std::array<std::int64_t,3> core
    {std::min(ssize[0]/2,dsize[0]),
     std::min(ssize[1]/2,dsize[1]),
     std::min(ssize[2]/2,dsize[2])};

  // The samples of the output where between 1 and 4 inputs are available.
  const std::array<std::int64_t,3> more
    {std::min((ssize[0]+1)/2,dsize[0]),
     std::min((ssize[1]+1)/2,dsize[1]),
     std::min((ssize[2]+1)/2,dsize[2])};

  // Caller is not expected to send a too large output buffer. If it
  // does so anyway then initialize the entire output buffer to zero.
  // Performance is irrelevant since this isn't supposed to happen.
  const bool needpad = more[0]<dsize[0] || more[1]<dsize[1] || more[2]<dsize[2];
  if (needpad)
    for (std::int64_t ii = 0; ii < dsize[0]; ++ii)
      for (std::int64_t jj = 0; jj < dsize[1]; ++jj)
        for (std::int64_t kk = 0; kk < dsize[2]; ++kk)
          dst[ii*dstride[0] + jj*dstride[1] + kk*dstride[2]] = 0;

  // Given offset of first input, where are the other 8?
  const std::array<std::int64_t,8> offsets
    {0*sstride[0] + 0*sstride[1] + 0*sstride[2],
     0*sstride[0] + 0*sstride[1] + 1*sstride[2],
     0*sstride[0] + 1*sstride[1] + 0*sstride[2],
     0*sstride[0] + 1*sstride[1] + 1*sstride[2],
     1*sstride[0] + 0*sstride[1] + 0*sstride[2],
     1*sstride[0] + 0*sstride[1] + 1*sstride[2],
     1*sstride[0] + 1*sstride[1] + 0*sstride[2],
     1*sstride[0] + 1*sstride[1] + 1*sstride[2]};

  for (std::int64_t ii = 0; ii < core[0]; ++ii) {
    for (std::int64_t jj = 0; jj < core[1]; ++jj) {
      const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1];
      T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1];
      for (std::int64_t kk = 0; kk < core[2]; ++kk,
             cursrc += 2*sstride[2], curdst += dstride[2]) {
        // Invoke the functor to calculate a single sample from 8 input samples.
        // The reason the cast is safe: While the product of size[] is int64,
        // size in each dimension is limited by int32. The size is declared
        // as int64_t[] simply to reduce the risk of forgetting to widen if
        // computing e.g. size[0]*size[1]*size[2].
        *curdst = function(cursrc[offsets[0]], cursrc[offsets[1]],
                           cursrc[offsets[2]], cursrc[offsets[3]],
                           cursrc[offsets[4]], cursrc[offsets[5]],
                           cursrc[offsets[6]], cursrc[offsets[7]]);
      }
    }
  }

  // 2d: Handle last il, if count was odd, where both xl and slices ok.
  if (more[0] > core[0]) {
    const std::int64_t ii = core[0]; {
      for (std::int64_t jj = 0; jj < core[1]; ++jj) {
        for (std::int64_t kk = 0; kk < core[2]; ++kk) {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last inline, offsets 4,5,6,7 are invalid.
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[1]],
                             cursrc[offsets[2]], cursrc[offsets[3]],
                             cursrc[offsets[4-4]], cursrc[offsets[5-4]],
                             cursrc[offsets[6-4]], cursrc[offsets[7-4]]);
        }
      }
    }
  }

  // 2d: Handle last xl, if count was odd, where both il, slice ok
  if (more[1] > core[1]) {
    for (std::int64_t ii = 0; ii < core[0]; ++ii) {
      const std::int64_t jj = core[1]; {
        for (std::int64_t kk = 0; kk < core[2]; ++kk) {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last crossline, offsets 2,3,6,7 are invalid.
          *curdst = function(cursrc[offsets[0]],   cursrc[offsets[1]],
                             cursrc[offsets[2-2]], cursrc[offsets[3-2]],
                             cursrc[offsets[4]],   cursrc[offsets[5]],
                             cursrc[offsets[6-2]], cursrc[offsets[7-2]]);
        }
      }
    }
  }

  // 2d: Handle last slice, if count was odd, where both il, xl ok
  if (more[2] > core[2]) {
    for (std::int64_t ii = 0; ii < core[0]; ++ii) {
      for (std::int64_t jj = 0; jj < core[1]; ++jj) {
        const std::int64_t kk = core[2]; {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last slice, all odd-numbered offsets are invalid.
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[1-1]],
                             cursrc[offsets[2]], cursrc[offsets[3-1]],
                             cursrc[offsets[4]], cursrc[offsets[5-1]],
                             cursrc[offsets[6]], cursrc[offsets[7-1]]);
        }
      }
    }
  }

  // 1d: Last il and xl, only Z varies.
  if (more[0] > core[0] && more[1] > core[1]) {
    const std::int64_t ii = core[0]; {
      const std::int64_t jj = core[1]; {
        for (std::int64_t kk = 0; kk < core[2]; ++kk) {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last inline, offsets 4,5,6,7 are invalid.
          // At the last crossline, offsets 2,3,6,7 are invalid.
          // So all that is left here is 0,1.
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[1]],
                             cursrc[offsets[0]], cursrc[offsets[1]],
                             cursrc[offsets[0]], cursrc[offsets[1]],
                             cursrc[offsets[0]], cursrc[offsets[1]]);
        }
      }
    }
  }

  // 1d: Last il and Z, only xl varies.
  if (more[0] > core[0] && more[2] > core[2]) {
    const std::int64_t ii = core[0]; {
      for (std::int64_t jj = 0; jj < core[1]; ++jj) {
        const std::int64_t kk = core[2]; {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last inline, offsets 4,5,6,7 are invalid.
          // At the last slice, all odd-numbered offsets are invalid.
          // What is left here is 0,2.
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[2]], cursrc[offsets[2]],
                             cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[2]], cursrc[offsets[2]]);
        }
      }
    }
  }

  // 1d: Last xl and Z, only il varies.
  if (more[1] > core[1] && more[2] > core[2]) {
    for (std::int64_t ii = 0; ii < core[0]; ++ii) {
      const std::int64_t jj = core[1]; {
        const std::int64_t kk = core[2]; {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          // At the last crossline, offsets 2,3,6,7 are invalid.
          // At the last slice, all odd-numbered offsets are invalid.
          // What is left here is 0,4.
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[4]], cursrc[offsets[4]],
                             cursrc[offsets[4]], cursrc[offsets[4]]);
        }
      }
    }
  }

  // 0d: Last il, xl and Z, single sample only.
  if (more[0] > core[0] && more[1] > core[1] && more[2] > core[2]) {
    const std::int64_t ii = core[0]; {
      const std::int64_t jj = core[1]; {
        const std::int64_t kk = core[2]; {
          const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1] + 2*kk*sstride[2];
          T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1] +   kk*dstride[2];
          *curdst = function(cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[0]], cursrc[offsets[0]],
                             cursrc[offsets[0]], cursrc[offsets[0]]);
        }
      }
    }
  }
}

/**
 * \brief Lowpass decimation vertically, simple decimation horizontally.
 *
 * Unlike the other algorithms this one needs full traces as input.
 * Or at least prefers full traces, to avoid brick artifacts vertically.
 */
template <typename T>
static void
createGenericLevelOfDetailLoPass(
    T* dst,
    const std::array<std::int64_t,3>&dsize,
    const std::array<std::int64_t,3>&dstride,
    const T* src,
    const std::array<std::int64_t,3>&ssize,
    const std::array<std::int64_t,3>&sstride)
{
  //printf("@createLODLoPass<%s> size (%d, %d, %d) -> (%d, %d, %d)\n",
  //       typeid(T).name(),
  //       (int)ssize[0], (int)ssize[1], (int)ssize[2],
  //       (int)dsize[0], (int)dsize[1], (int)dsize[2]);

  for (int ii=0; ii<3; ++ii)
    if ((ssize[ii]+1) / 2 != dsize[ii])
      throw OpenZGY::Errors::ZgyInternalError("createLODLoPass buffer size mismatch.");

  std::unique_ptr<double[]> vecSrc (new double[std::max(ssize[2], 2*dsize[2])]);
  std::unique_ptr<double[]> vecDst (new double[dsize[2]]);

  for (std::int64_t ii = 0; ii < dsize[0]; ++ii) {
    for (std::int64_t jj = 0; jj < dsize[1]; ++jj) {
      const T* cursrc = src + 2*ii*sstride[0] + 2*jj*sstride[1];
      T*       curdst = dst +   ii*dstride[0] +   jj*dstride[1];
      // Copy a single trace of source data. If not enough source
      // (need twice the output size) then replicate the last value.
      // This supports the case where ssize[2] is odd.
      for (std::int64_t kk = 0; kk < ssize[2]; ++kk)
        vecSrc[kk] = cursrc[kk*sstride[2]];
      for (std::int64_t kk = ssize[2]; kk < 2*dsize[2]; ++kk)
        vecSrc[kk] = vecSrc[ssize[2]-1];
      // Filter and downsample vector.
      // The reason the cast is safe: While the product of size[] is int64,
      // size in each dimension is limited by int32. The size is declared
      // as int64_t[] simply to reduce the risk of forgetting to widen if
      // computing e.g. size[0]*size[1]*size[2].
      LodSampling::downSample1D(
        vecDst.get(),
        static_cast<int>(dsize[2]),
        vecSrc.get(),
        static_cast<int>(2*dsize[2]));
      // Convert back from double (which the lowpass filter uses)
      // to the storage valuetype. RoundAndClip() takes care of
      // clipping to the limits of T, just in case the result ended
      // up with slightly larger amplitude than the original and
      // that result also overflowed the target type. No clipping
      // is done for floats.
      for (std::int64_t kk = 0; kk < dsize[2]; ++kk)
        curdst[kk*dstride[2]] = RoundAndClip<T>(vecDst[kk]);
    }
  }
}

/**
 * Create a single LOD brick based on the data from one LOD level below.
 * The valuetype of the source and target must be the same, but can be almost
 * any scalar type.
 *
 * Histogram information is passed in for use by some of the algorithms.
 * The histogram should be in "storage" units. That only affects its
 * limits though. TODO-Worry make sure these are set correctly. The old accessor
 * assumed that for integral data the limits were simply the entire range.
 *
 * When doing lowpass filtering of integral data, the actual calculation is
 * done using doubles.  We need to clip the result to the range of the data
 * value type to prevent overflows, since the lowpass filter can give a result
 * slightly outside the input range.
 *
 * When doing lowpass filtering of float data there will be no clipping.
 * Some samples might end up slightly outside the actual value range of the
 * input data. Note behavior change: The old accessor did some clipping to
 * avoid that case. Which just added to the confusion because the range it
 * clipped to might be inaccurate.
 */
template <typename T>
static void
createLevelOfDetail(
    T* dst,
    const std::array<std::int64_t,3>&dsize,
    const std::array<std::int64_t,3>&dstride,
    const T* src,
    const std::array<std::int64_t,3>&ssize,
    const std::array<std::int64_t,3>&sstride,
    LodAlgorithm algorithm,
    const std::int64_t *hist, int bins, double minHist, double maxHist)
{
    switch (algorithm) {
    case LodAlgorithm::LowPass:
      createGenericLevelOfDetailLoPass(dst, dsize, dstride, src, ssize, sstride);
      break;

    case LodAlgorithm::WeightedAverage:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodWeightedAverage<T>(hist, bins, minHist, maxHist));
      break;

    case   LodAlgorithm::Average:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodAverage<T>());
      break;

    case   LodAlgorithm::Median:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMedian<T>());
      break;

    case   LodAlgorithm::Minimum:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMinimum<T>());
      break;

    case   LodAlgorithm::Maximum:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMaximum<T>());
      break;

    //case   LodAlgorithm::MinMax:
      //createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMinMax<T>());
      //break;

    case   LodAlgorithm::Decimate:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodDecimate<T>());
      break;

    case   LodAlgorithm::DecimateSkipNaN:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodDecimateSkipNaN<T>());
      break;

    //case   LodAlgorithm::DecimateRandom:
      //createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodDecimateRandom<T>());
      //break;

    case   LodAlgorithm::AllZero:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodAllZero<T>());
      break;

    //case   LodAlgorithm::WhiteNoise:
      //createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodWhiteNoise<T>());
      //break;

    case   LodAlgorithm::MostFrequent:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMostFrequent<T, !std::numeric_limits<T>::is_integer, false>());
      break;

    case   LodAlgorithm::MostFrequentNon0:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodMostFrequent<T, !std::numeric_limits<T>::is_integer, true>());
      break;

    case   LodAlgorithm::AverageNon0:
      createGenericLevelOfDetail(dst, dsize, dstride, src, ssize, sstride, LodAverageNon0<T>());
      break;

          default:
      throw OpenZGY::Errors::ZgyInternalError("Unsupported decimation type");
    }
}

/**
 * New version with "modernized" types.
 * Not sure if it belongs in this file or higher up.
 */
template <typename T>
static void
createLodT(const std::shared_ptr<DataBuffer>& result,
           const std::shared_ptr<const DataBuffer>& input,
           LodAlgorithm algorithm,
           const std::int64_t* hist,
           std::int32_t bincount,
           double histogram_min,
           double histogram_max)
{
  typedef DataBufferNd<T,3> buffer_t;
  typedef const DataBufferNd<T,3> constbuffer_t;
  std::shared_ptr<buffer_t> resultT =
    std::dynamic_pointer_cast<buffer_t>(result);
  std::shared_ptr<constbuffer_t> inputT =
    std::dynamic_pointer_cast<constbuffer_t>(input);
  if (result == nullptr || input == nullptr)
    throw OpenZGY::Errors::ZgyInternalError("createLodT null buffer(s) passed.");
  if (resultT == nullptr)
    throw OpenZGY::Errors::ZgyInternalError("createLodT wrong result data type.");
  if (inputT == nullptr)
    throw OpenZGY::Errors::ZgyInternalError("createLodT wrong input data type.");
  createLevelOfDetail(resultT->data(), resultT->size3d(), resultT->stride3d(),
                      inputT->data(), inputT->size3d(), inputT->stride3d(),
                      algorithm,
                      hist, bincount, histogram_min, histogram_max);
}

/**
 * \brief Main entry point for low resolution compute.
 *
 * Create a single low resolution brick 1/8th the size of the input.
 * This is the single threaded version.
 */
static void
createLodST(const std::shared_ptr<DataBuffer>& result,
            const std::shared_ptr<const DataBuffer>& input,
            LodAlgorithm algorithm,
            const std::int64_t* hist,
            std::int32_t bincount,
            double histogram_min,
            double histogram_max)
{
  static SummaryPrintingTimerEx timer("createLod");
  SimpleTimerEx tt(timer);
  switch (result->datatype()) {
  case RawDataType::SignedInt8:    createLodT<std::int8_t>  (result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::UnsignedInt8:  createLodT<std::uint8_t> (result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::SignedInt16:   createLodT<std::int16_t> (result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::UnsignedInt16: createLodT<std::uint16_t>(result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::SignedInt32:   createLodT<std::int32_t> (result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::UnsignedInt32: createLodT<std::uint32_t>(result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::Float32:       createLodT<float>        (result, input, algorithm, hist, bincount, histogram_min, histogram_max); break;
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Run the low resolution compute on just a part of the provided
 * data buffer. The part will be one of more slices in the slowest
 * changing direction. If passing the last 3 parameters as "entire survey"
 * then this is equivalent to calling createLodST() directly.
 */
static void
createLodPart(const std::shared_ptr<DataBuffer>& result,
              const std::shared_ptr<const DataBuffer>& input,
              LodAlgorithm algorithm,
              const std::int64_t* hist,
              std::int32_t bincount,
              double histogram_min,
              double histogram_max,
              int slice_dim,
              std::int64_t slice_beg,
              std::int64_t slice_count)
{
  // If at the end of the survey, we might process fewer slices. Or zero
  // slices e.g. if slices_per_thread was rounded up to make it even.
  const std::int64_t survey_end = input->size3d()[slice_dim];
  slice_count = std::min(slice_beg + slice_count, survey_end) - slice_beg;
  if (slice_count <= 0)
    return;
  // Cannot happen. Unless a programmer messes up...
  // Note that slice_count may be odd if we are at the end of the survey.
  // In that case, remember that the output size rounds up.
  if (slice_beg % 2 != 0)
    throw OpenZGY::Errors::ZgyInternalError("Slice start must be even.");
  // Create a view on the buffers so the actual lod generator still
  // thinks we start at zero.
  auto ibuf = input->slice1(slice_dim, slice_beg, slice_count);
  auto obuf = result->slice1(slice_dim, slice_beg/2, (slice_count+1)/2);
  createLodST(obuf, ibuf, algorithm,
              hist, bincount, histogram_min, histogram_max);
}

/**
 * \brief Main entry point for low resolution compute.
 *
 * Create a single low resolution brick 1/8th the size of the input.
 * This is the multi-threaded version.
 *
 * Caveat: If parallelizing is also done at a higher level then
 * consider omp_set_max_active_levels(), especially if the higher
 * level might only be able to use 2 or 4 threads.
 *
 * One argument for parallelizing on a higher level instead
 * is that we might then be able to read and compute at the
 * same time. N/A for local file access because buffer cache
 * prefectch can give us that automatically. I think.
 *
 * TODO-Low fall back to single threaded if the input is small
 * (e.g. less than 16 slices or less than 64^3 samples) or if
 * the caller asks for cheaper algorithms. Low priority because
 * those cases are (a) rare and (b) should run fast anyway, even
 * if single thread *might* give a slight improvement.
 *
 * Using a simple "copy" app between two local ssd files I measured
 * the time for lod compute as 4% of finalize time or 2% of total time
 * with 8 threads. Compared to 23% / 13% for a single thread. So the
 * MT case had a 10% overall speedup on this test.
 *
 * Switching to the much cheaper "Decimate" algorithm (just pick one
 * of every 8 samples) then even in the single threaded case the
 * algorithm only accounts for < 1% so there is not much to gain.
 *
 * | Algorithm(s) used for Onnia | Time     |
 * | :-------------------------- | -------: |
 * | Decimate                    |    8 sec |
 * | Average                     |   37 sec |
 * | LowPass                     |   72 sec |
 * | WeightedAverage             |  732 sec |
 * | Decimate,WeightedAverage    |  102 sec |
 * | LowPass,Decimate            |   64 sec |
 * | LowPass,WeightedAverage     |  159 sec |
 */
static void
createLodMT(const std::shared_ptr<DataBuffer>& result,
            const std::shared_ptr<const DataBuffer>& input,
            LodAlgorithm algorithm,
            const std::int64_t* hist,
            std::int32_t bincount,
            double histogram_min,
            double histogram_max)
{
  static SummaryPrintingTimerEx timerST("createLod[ST]");
  static SummaryPrintingTimerEx timerMT("createLod[MT]");
  const auto isize   = input->size3d();
#if 1
  // If this is a 2d slice, regardless of which dimension had size 1,
  // there is probably not enough data to warrant slicing it up and
  // processing it in multiple threads.
  if (isize[0] == 1 || isize[1] == 1 || isize[2] == 1) {
    SimpleTimerEx tt(timerST);
    createLodST(result, input, algorithm,
                hist, bincount, histogram_min, histogram_max);
    return;
  }
#endif
  // Now that DataBuffer is always c-contiguous, this is simpler.
  // TODO-Test: The 2d case case where size[2]==1.
  const int slowest_dim = isize[0] == 1 ? 1 : 0;
  // Technically this should work with any buffer layouts. But that
  // means a lot of extra testing, less efficiency, and caveats such
  // as updating a byte might be implemented as read/modify/write of
  // 4 or 8 bytes. So, always slice the slowest dim.
  SimpleTimerEx tt(timerMT);
  MTGuard guard("lod", -1);
#pragma omp parallel
  {
    std::int64_t slices_per_thread = (isize[slowest_dim]-1) / omp_get_num_threads() + 1;
    if (slices_per_thread % 2 == 1)
      ++slices_per_thread;
#if 0
    if (omp_get_thread_num() == 0)
      std::cout << "LOD compute " << isize[slowest_dim]
                << " slices using " << omp_get_num_threads()
                << " threads with " << slices_per_thread
                << " slices per thread."
                << std::endl;
#endif
    guard.run([&](){
      createLodPart(result, input, algorithm,
                    hist, bincount, histogram_min, histogram_max,
                    slowest_dim,
                    omp_get_thread_num() * slices_per_thread,
                    slices_per_thread);
    });
  }
  guard.finished();
}

/**
 * \brief Main entry point for low resolution compute.
 *
 * Create a single low resolution brick 1/8th the size of the input.
 * Decides whether to enable multi-threading or not.
 */
void
createLod(const std::shared_ptr<DataBuffer>& result,
          const std::shared_ptr<const DataBuffer>& input,
          LodAlgorithm algorithm,
          const std::int64_t* hist,
          std::int32_t bincount,
          double histogram_min,
          double histogram_max)
{
  // Lowpass on a trace just a few samples long makes no sense, and
  // would hit a corner case and an assert.
  if (algorithm == LodAlgorithm::LowPass && input->size3d()[2] < 5)
    algorithm = LodAlgorithm::Decimate;

  // TODO-Worry: A lot more can go wrong in the MT case.
  // Have I tested this enough?
  createLodMT(result, input, algorithm,
              hist, bincount, histogram_min, histogram_max);
}

} // namespace
