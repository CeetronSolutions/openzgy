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

#pragma once

#include "../declspec.h"

#include <cstdint>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file lodfilters.h
 * \brief Static methods for downsampling. Used by lodalgo.cpp only.
 * \details
 * This is an improved version of class LodSampling.
 */

/**
 * \brief Static methods for downsampling. Used by lodalgo.cpp only.
 * \details Thread safety: Safe because this is a static class with no data.
 */
class OPENZGY_TEST_API LodFilters
{
  LodFilters() = delete;
  LodFilters(const LodFilters&) = delete;
  LodFilters& operator=(const LodFilters&) = delete;
  public:

 /**
  * Downsample a vector with a factor of 2.
  * The vector is low pass filtered and then down sampled.
  * \param[out] dst Out vector with downsampled result
  * \param[in] sizeDst Size of output vector. It must be sizeDst*2==sizeSrc
  * \param[in] src In vector
  * \param[in] sizeSrc Size of src vector. It must be sizeSrc==sizeDst*2 and >=5
  */
  static void downSample1D(double* dst, int sizeDst, const double* src, int sizeSrc, const double* src4before = nullptr, const double* src4after = nullptr);
  static void downSample1D(float*,         int, const float*,         int, const float* src4before = nullptr, const float* src4after = nullptr);
  static void downSample1D(std::int32_t*,  int, const std::int32_t*,  int, const std::int32_t*  src4before = nullptr, const std::int32_t*  src4after = nullptr);
  static void downSample1D(std::int16_t*,  int, const std::int16_t*,  int, const std::int16_t*  src4before = nullptr, const std::int16_t*  src4after = nullptr);
  static void downSample1D(std::int8_t*,   int, const std::int8_t*,   int, const std::int8_t*   src4before = nullptr, const std::int8_t*   src4after = nullptr);
  static void downSample1D(std::uint32_t*, int, const std::uint32_t*, int, const std::uint32_t* src4before = nullptr, const std::uint32_t* src4after = nullptr);
  static void downSample1D(std::uint16_t*, int, const std::uint16_t*, int, const std::uint16_t* src4before = nullptr, const std::uint16_t* src4after = nullptr);
  static void downSample1D(std::uint8_t*,  int, const std::uint8_t*,  int, const std::uint8_t*  src4before = nullptr, const std::uint8_t*  src4after = nullptr);
};

} // namespace
