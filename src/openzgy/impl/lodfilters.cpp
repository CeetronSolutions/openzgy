// Copyright 2017-2022, Schlumberger
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

//Adapted from: Zgy/Common/Sampling.cpp

#include "lodfilters.h"
#include "roundandclip.h"

#include <float.h> // needed for _finite() on Windows
#include <algorithm>

namespace InternalZGY {
#if 0
}
#endif

namespace
{
  /**
   * Return the average of two values.
   * If one of them is NaN or Infinite, return the other one.
   * If both are NaN/Infinite, return the first.
   */
  template<typename T>
  inline T average(T a, T b)
  {
    T val = a/2 + b/2;
    if (!IsFiniteT<T>(val))
      val = b;
    if (!IsFiniteT<T>(val))
      val = a;
    return val;
  }
}

/*
  Parks-McClellan FIR Filter Design

  Filter type: Low pass
  Passband: 0 - 0.22
  Order: 9
  Passband ripple: 1.0 dB
  Transition band: 0.07
  Stopband attenuation: 20.0 dB

  Coefficients:

  a[0] =   0.07996591908598821
  a[1] =  -0.06968050331585399
  a[2] =  -0.10596589191287473
  a[3] =   0.12716426995479813
  a[4] =   0.4496382067876111
  a[5] =   0.4496382067876111
  a[6] =   0.12716426995479813
  a[7] =  -0.10596589191287473
  a[8] =  -0.06968050331585399
  a[9] =   0.07996591908598821

  UPDATED:

  The filter was defined without requiring a 0 dB DC response,
  which is really bad for the purpose we are using it for.
  We need the sum of the coefficients to be exactly 1.0
  Instead of designing a new filter, multiply each coefficient
  with the same fudge factor to make the sum 1.0.

  0.46728086247062456
  0.13215387136350137
 -0.11012372306899207
 -0.07241458842975841
  0.08310357766462453

  UPDATED:

  For the purpose of computing low resolution bricks, doing the
  calculation in single precision float ought to be more than enough.
  Also, the old code had several helper variables to make it easier
  to replace the filter with something of a different rank. That
  won't happen -- and it wouldn't work anyway to just change the
  variables. So, do more hard coding of "magic" sizes.

  0.46728086f
  0.13215387f
 -0.11012372f
 -0.07241458f
  0.08310357f
*/


namespace {
#if 0
}
#endif

/**
 * Run a lowpass filter on the provided trace. The output will have
 * half the number of samples as the input, rounded upwards.
 *
 * Normally the trace is mirrored to get 4 samples before the start
 * of the trace and 4 samples past the last trace. This is needed
 * by the lowpass filter. The caller can avoid the mirroring by
 * providing the 4 before-start and/or after-end samples.
 * All 4 samples must be provided. If caller has fewer, pass null.
 *
 * The reason for this feature is that the caller might have read the
 * data to be decimated as raw bricks. For performance reasons it is
 * better to skip the step where the buffers are reshaped to make
 * the traces contiguous. Decimation will be run on (typically) 64
 * traces input / 32 traces output at a time. The feature avoids
 * getting artefacts every 32 traces.
 *
 * OpenZGY currently reshapes the buffers, so it won't use this feature.
 * This might change shortly.
 */
template<typename T, typename U>
void downSampleT(T* dst, int sizeDst, const T* src, int sizeSrc, const T* src4before, const T* src4after)
{
  // Low pass filter (See description of filter parameters above.)
  U a[] = {
    (U) 0.46728086247062456,
    (U) 0.13215387136350137,
    (U)-0.11012372306899207,
    (U)-0.07241458842975841,
    (U) 0.08310357766462453};

  // Technically unneeded test, but in almost all cases it will fail
  // which means that on average it saves a few cpu cycles.
  if (sizeSrc < 4 || sizeSrc != sizeDst*2) {
    // Never need more than two input values for one output.
    // Sort of. A very minor issue: The algorithm could have used
    // a few extra source samples so that it could avoid the
    // reflection done at the end of the trace. But, really,
    // a too small output buffer just shouldn't happen in the
    // first place.
    sizeSrc = std::max(0, std::min(sizeSrc, sizeDst*2));
    // Cannot produce more than one output value from two inputs,
    // except for allowing an odd number of source bytes and
    // fudging the last output value.
    sizeDst = std::max(0, std::min(sizeDst, (sizeSrc+1)/2));

    // assert (sizeSrc == sizeDst*2 || sizeSrc == sizeDst*2 - 1);

    // Handle corner cases.
    // These are rare enough that testing for src4before and
    // src4after isn't useful.
    switch (sizeSrc) {
    case 0:
      return;
    case 1:
      dst[0] = src[0];
      return;
    case 2:
      dst[0] = average<T>(src[0], src[1]);
      return;
    case 3:
      dst[0] = average<T>(src[0], src[1]);
      dst[1] = average<T>(src[1], src[2]);
      return;
    default:
      break;
    }

    // Handle odd input size when the destination size was rounded up.
    // When using the filter for lowres computation it doesn't really
    // matter how we choose the output as long as it isn't completely
    // bogus. It only affects the bottom few samples in each trace.
    if (sizeSrc & 1) {
      // Already know that sizeSrc == 2*sizeDst - 1, i.e. one missing.
      dst[sizeDst-1] = average<T>(src[sizeSrc-2], src[sizeSrc-1]);
      --sizeSrc;
      --sizeDst;
    }
  } // end of corner cases.

  // Edges are mirrored to reduce edge artefacts.
  // Start edge mirror filtering
  if (!src4before) {
    dst[0] = RoundAndClip<T>
               (a[0]*((U)src[0]+(U)src[1]) +
                a[1]*((U)src[1]+(U)src[2]) +
                a[2]*((U)src[2]+(U)src[3]) +
                a[3]*((U)src[3]+(U)src[4]) +
                a[4]*((U)src[4]+(U)src[5]));
    dst[1] = RoundAndClip<T>
               (a[0]*((U)src[2]+(U)src[3]) +
                a[1]*((U)src[1]+(U)src[4]) +
                a[2]*((U)src[0]+(U)src[5]) +
                a[3]*((U)src[1]+(U)src[6]) +
                a[4]*((U)src[2]+(U)src[7]));
  }
  else {
    dst[0] = RoundAndClip<T>
               (a[0]*((U)src[0]+(U)src[1]) +
                a[1]*((U)src4before[3]+(U)src[2]) +
                a[2]*((U)src4before[2]+(U)src[3]) +
                a[3]*((U)src4before[1]+(U)src[4]) +
                a[4]*((U)src4before[0]+(U)src[5]));
    dst[1] = RoundAndClip<T>
               (a[0]*((U)src[2]+(U)src[3]) +
                a[1]*((U)src[1]+(U)src[4]) +
                a[2]*((U)src[0]+(U)src[5]) +
                a[3]*((U)src4before[3]+(U)src[6]) +
                a[4]*((U)src4before[2]+(U)src[7]));
  }

  // If any of the first 5 values were NaN or Infinite,
  // fall back to using a simpler algorithm.
  if (!IsFiniteT(dst[0]))
    dst[0] = average(src[0], src[1]);

  if (!IsFiniteT(dst[1]))
    dst[1] = average(src[2], src[3]);

  const int e = sizeSrc-1;

  if (!src4after) {
    dst[sizeDst-2] = RoundAndClip<T>
               (a[0]*((U)src[e-3]+(U)src[e-2]) +
                a[1]*((U)src[e-4]+(U)src[e-1]) +
                a[2]*((U)src[e-5]+(U)src[e-0]) +
                a[3]*((U)src[e-6]+(U)src[e-1]) +
                a[4]*((U)src[e-7]+(U)src[e-2]));
    dst[sizeDst-1] = RoundAndClip<T>
               (a[0]*((U)src[e-1]+(U)src[e-0]) +
                a[1]*((U)src[e-2]+(U)src[e-1]) +
                a[2]*((U)src[e-3]+(U)src[e-2]) +
                a[3]*((U)src[e-4]+(U)src[e-3]) +
                a[4]*((U)src[e-5]+(U)src[e-4]));
  }
  else {
    dst[sizeDst-2] = RoundAndClip<T>
               (a[0]*((U)src[e-3]+(U)src[e-2]) +
                a[1]*((U)src[e-4]+(U)src[e-1]) +
                a[2]*((U)src[e-5]+(U)src[e-0]) +
                a[3]*((U)src[e-6]+(U)src4after[0]) +
                a[4]*((U)src[e-7]+(U)src4after[1]));
    dst[sizeDst-1] = RoundAndClip<T>
               (a[0]*((U)src[e-1]+(U)src[e-0]) +
                a[1]*((U)src[e-2]+(U)src4after[0]) +
                a[2]*((U)src[e-3]+(U)src4after[1]) +
                a[3]*((U)src[e-4]+(U)src4after[2]) +
                a[4]*((U)src[e-5]+(U)src4after[3]));
  }

  if (!IsFiniteT(dst[sizeDst-2]))
    dst[sizeDst-2] = average(src[e-3], src[e-2]);
  if (!IsFiniteT(dst[sizeDst-1]))
    dst[sizeDst-1] = average(src[e-1], src[e-0]);

  // Filter and downsample vector (excluding edges).
  for (int i=2, ii=4; i<sizeDst-2; ++i, ii+=2){

    dst[i] = RoundAndClip<T>(
      a[0]*((U)src[ii-0]+(U)src[ii+1]) +
      a[1]*((U)src[ii-1]+(U)src[ii+2]) +
      a[2]*((U)src[ii-2]+(U)src[ii+3]) +
      a[3]*((U)src[ii-3]+(U)src[ii+4]) +
      a[4]*((U)src[ii-4]+(U)src[ii+5]));

    // Verify the last iteration: i = sizeDst-3, ii = 2*sizeDst-6
    // which means that the last output produced is dst[sizeDst-3]
    // (leaving two missing) and the last input consumed is
    // src[(2*sizeDst-6)+5] = src[2*sizeDst-1] = src[sizeSrc-1].
    // The last two output values are handled by the edge compute.

    // If any of the values ended up infinite, use a simpler approach.
    // Performance note: We could skip this test if caller can
    // guarantee that src does not contain and NaN values, e.g.
    // because it was converted from an integral type. But we
    // probably wouldn't gain much by it.
    if (!IsFiniteT(dst[i]))
      dst[i] = average(src[ii], src[ii+1]);
  }
}
} // namespace

/**
 * This overload is the one currently used. Caller copies the data to a
 * temporary buffer and copies it back afterwards. And it handles most
 * corner cases itself. It should be more efficient to use one of the
 * other overloads and operate directly on the input and output cubes.
 * That requires step[2]==1, but that ought to be a safe assumption.
 */
void
LodFilters::downSample1D(double* dst, int sizeDst, const double* src, int sizeSrc, const double* src4before, const double* src4after)
{
  downSampleT<double, double>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(float* dst, int sizeDst, const float* src, int sizeSrc, const float* src4before, const float* src4after)
{
  // Using float intermediate values (the second template parameter) means
  // we don't need to cast from T to U which happens 10 times for each sample.
  // TODO-Performance: The cost in this case is that RoundAndClip(double)
  // now costs 2 casts (to double and back). An overloaded RoundAndClip(float)
  // would avoid this. But might not be worth the trouble.
  downSampleT<float, float>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::int32_t* dst, int sizeDst, const std::int32_t* src, int sizeSrc, const std::int32_t* src4before, const std::int32_t* src4after)
{
  // With int32 input, intermediate value must be double.
  // TODO YAGNI: int32 and all unsigned are not supported
  // and are unlikely to ever be. Remove those overloads,
  // and also remove their use in lodalgo.cpp
  downSampleT<std::int32_t, double>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::int16_t* dst, int sizeDst, const std::int16_t* src, int sizeSrc, const std::int16_t* src4before, const std::int16_t* src4after)
{
  // TODO Performance: Measure whether float or double is faster
  // for the intermediate calculation. Precision wise it doesn't
  // matter for 8- and 16-bit final result.
  downSampleT<std::int16_t, float>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::int8_t* dst, int sizeDst, const std::int8_t* src, int sizeSrc, const std::int8_t* src4before, const std::int8_t* src4after)
{
  downSampleT<std::int8_t, float>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::uint32_t* dst, int sizeDst, const std::uint32_t* src, int sizeSrc, const std::uint32_t* src4before, const std::uint32_t* src4after)
{
  // With int32 input, intermediate value must be double.
  downSampleT<std::uint32_t, double>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::uint16_t* dst, int sizeDst, const std::uint16_t* src, int sizeSrc, const std::uint16_t* src4before, const std::uint16_t* src4after)
{
  downSampleT<std::uint16_t, float>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}

void
LodFilters::downSample1D(std::uint8_t* dst, int sizeDst, const std::uint8_t* src, int sizeSrc, const std::uint8_t* src4before, const std::uint8_t* src4after)
{
  downSampleT<std::uint8_t, float>(dst, sizeDst, src, sizeSrc, src4before, src4after);
}


} // end namespace
