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
#pragma once

#include "../declspec.h"
#include "enum.h"

#include <cstdint>
#include <memory>
#include <array>

/**
 * \file lodalgo.h
 * \brief Decimation algorithms to output low resolution bricks.
 */

namespace InternalZGY {
#if 0
}
#endif

class DataBuffer;
template<typename T, int NDim> class DataBufferNd;

/**
 * Possible algorithms to generate LOD bricks.
 * This is work in progress, and not all these
 * algorithms will get implemented. The ones
 * that are not needed will eventually get trimmed
 * from this enum.
 */
enum class LodAlgorithm
{
  LowPass          =  0, ///< Lowpass Z / decimate XY
  WeightedAverage  =  1, ///< Weighted averaging (depends on global stats)
  Average          =  2, ///< Simple averaging
  Median           =  3, ///< Somewhat more expensive averaging
  Minimum          =  4, ///< Minimum value
  Maximum          =  5, ///< Maximum value
  //MinMax         =  6, ///< Checkerboard of minimum and maximum values
  Decimate         =  7, ///< Simple decomation, use first sample
  DecimateSkipNaN  =  8, ///< Use first sample that is not NaN
  //DecimateRandom =  9, ///< Random decimation using a fixed seed
  AllZero          = 10, ///< Just fill the LOD brick with zeroes
  //WhiteNoise     = 11, ///< Fill with white noise
  MostFrequent     = 12, ///< The value that occurs most frequently
  MostFrequentNon0 = 13, ///< The non-zero value that occurs most frequently
  AverageNon0      = 14, ///< Average value, but treat 0 as NaN.
};

void OPENZGY_TEST_API createLod(const std::shared_ptr<DataBuffer>& result,
               const std::shared_ptr<const DataBuffer>& input,
               LodAlgorithm algorithm,
               const std::int64_t* hist,
               std::int32_t bincount,
               double histogram_min,
               double histogram_max);

void OPENZGY_TEST_API clearLodTimers(bool show);

} // namespace
