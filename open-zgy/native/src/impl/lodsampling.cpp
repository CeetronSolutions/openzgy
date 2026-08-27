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

//Adapted from: Zgy/Common/Sampling.cpp

#include "lodsampling.h"
#include "roundandclip.h"

#include <assert.h>
#include <float.h> // needed for _finite() on Windows

namespace InternalZGY {
#if 0
}
#endif

namespace
{
  /**
   * Return the average of two doubles.
   * If one of them is NaN or Infinite, return the other one.
   * If both are NaN/Infinite, return the first.
   */
  inline double average(double a, double b)
  {
    double val = 0.5*(a+b);
    if (!std::isfinite(val))
      val = static_cast<double>(b);
    if (!std::isfinite(val))
      val = static_cast<double>(a);
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

  0.46728086f
  0.13215387f
 -0.11012372f
 -0.07241458f
  0.08310357f
*/


void LodSampling::downSample1D(double* dst, int sizeDst, const double* src, int sizeSrc)
{
  // Low pass filter (See description of filter parameters above.)
  // Note: The filter can only be replaced by a filter with even number of coefficients.
  double a[] = {
    0.46728086247062456,
    0.13215387136350137,
   -0.11012372306899207,
   -0.07241458842975841,
    0.08310357766462453};

  const int FILTER_LENGTH_HALF = sizeof(a)/sizeof(double);
  const int FILTER_LENGTH_HALF_MINUS_ONE = FILTER_LENGTH_HALF - 1;

  assert ( sizeSrc >= FILTER_LENGTH_HALF_MINUS_ONE );
  assert ( sizeSrc == sizeDst*2);

  // Edges are mirrored to reduce edge artefacts.
  double dstEdge[FILTER_LENGTH_HALF_MINUS_ONE];

  // Start edge mirror filtering
  dstEdge[0] =  a[0]*(src[0]+src[1]) +
                a[1]*(src[1]+src[2]) +
                a[2]*(src[2]+src[3]) +
                a[3]*(src[3]+src[4]) +
                a[4]*(src[4]+src[5]);
  // If any of the first 5 values were NaN or Infinite,
  // fall back to using a simpler algorithm.
  if (!std::isfinite(dstEdge[0]))
    dstEdge[0] = average(src[0], src[1]);

  dstEdge[1] =  a[0]*(src[2]+src[3]) +
                a[1]*(src[1]+src[4]) +
                a[2]*(src[0]+src[5]) +
                a[3]*(src[1]+src[6]) +
                a[4]*(src[2]+src[7]);
  if (!std::isfinite(dstEdge[1]))
    dstEdge[1] = average(src[2], src[3]);

  // End edge mirror filtering
  const int e = sizeSrc-1;
  dstEdge[2] =  a[0]*(src[e-3]+src[e-2]) +
                a[1]*(src[e-4]+src[e-1]) +
                a[2]*(src[e-5]+src[e-0]) +
                a[3]*(src[e-6]+src[e-1]) +
                a[4]*(src[e-7]+src[e-2]);
  if (!std::isfinite(dstEdge[2]))
    dstEdge[2] = average(src[e-3], src[e-2]);

  dstEdge[3] =  a[0]*(src[e-1]+src[e-0]) +
                a[1]*(src[e-2]+src[e-1]) +
                a[2]*(src[e-3]+src[e-2]) +
                a[3]*(src[e-4]+src[e-3]) +
                a[4]*(src[e-5]+src[e-4]);
  if (!std::isfinite(dstEdge[3]))
    dstEdge[3] = average(src[e-1], src[e-0]);

  // Set edge values
  for (int i=0; i<(FILTER_LENGTH_HALF_MINUS_ONE/2); ++i){
    dst[i] = dstEdge[i];
    dst[sizeDst-1-i] = dstEdge[FILTER_LENGTH_HALF-2-i];
  }

  // Filter and downsample vector (excluding edges).
  for (int i=2, ii=4; i<sizeDst-2; ++i, ii+=2){

    double val = a[0]*(src[ii]+src[ii+1]);
    for (int n=1;n<FILTER_LENGTH_HALF;++n){
      val += a[n]*(src[ii-n]+src[ii+1+n]);
    }

    // If any of the values ended up infinite, use a simpler approach.
    // Performance note: We could skip this test if caller can
    // guarantee that src does not contain and NaN values, e.g.
    // because it was converted from an integral type. But we
    // probably wouldn't gain much by it.
    if (!std::isfinite(val))
      val = average(src[ii], src[ii+1]);

    dst[i] = val;
  }
}

} // end namespace
