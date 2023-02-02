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

//Adapted from: Zgy/Common/Sampling.cpp

#pragma once

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file lodsampling.h
 * \brief Static methods for downsampling. Used by lodalgo.cpp only.
 */

/**
 * \brief Static methods for downsampling. Used by lodalgo.cpp only.
 * \details Thread safety: Safe because this is a static class with no data.
 */
class LodSampling
{
  LodSampling() = delete;
  LodSampling(const LodSampling&) = delete;
  LodSampling& operator=(const LodSampling&) = delete;
  public:

 /**
  * Downsample a vector with a factor of 2.
  * The vector is low pass filtered and then down sampled.
  * \param[out] dst Out vector with downsampled result
  * \param[in] sizeDst Size of output vector. It must be sizeDst*2==sizeSrc
  * \param[in] src In vector
  * \param[in] sizeSrc Size of src vector. It must be sizeSrc==sizeDst*2 and >=5
  */
  static void downSample1D(double* dst, int sizeDst, const double* src, int sizeSrc);
};

} // namespace
