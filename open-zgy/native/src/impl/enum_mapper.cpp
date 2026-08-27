// Copyright 2017-2024, Schlumberger
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

#include "enum_mapper.h"
#include "enum.h"    // InternalZGY enums
#include "lodalgo.h" // InternalZGY::LodAlgorithm
#include "../api.h"  // OpenZGY enums, should really be in a separate file.
#include "../exception.h"

namespace OpenZGY { namespace Impl {
#if 0
}}
#endif

/**
 * Map between enums used in the public API and the internal ones that
 * might change without notice and might be used to define the actual
 * numbers written to file.
 *
 * As a general rule, if an invalid enum tag is encountered while
 * converting from public to internal then this will throw an exceptiom
 * because it would be a user error. When converting from internal to
 * public the error is probably bad data encountered on the file.
 * The mapping function might just return "unknown", leaving to the
 * caller to decide whether this should be silently ignored.
 */
OpenZGY::SampleDataType
EnumMapper::mapRawDataTypeToSampleDataType(InternalZGY::RawDataType value)
{
  using InternalZGY::RawDataType;
  using OpenZGY::SampleDataType;
  switch (value) {
  case RawDataType::SignedInt8:    return SampleDataType::int8;
  case RawDataType::UnsignedInt8:  return SampleDataType::unknown;
  case RawDataType::SignedInt16:   return SampleDataType::int16;
  case RawDataType::UnsignedInt16: return SampleDataType::unknown;
  case RawDataType::SignedInt32:   return SampleDataType::unknown;
  case RawDataType::UnsignedInt32: return SampleDataType::unknown;
  case RawDataType::Float32:       return SampleDataType::float32;
  case RawDataType::IbmFloat32:    return SampleDataType::unknown;
  default:                         return SampleDataType::unknown;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
InternalZGY::RawDataType
EnumMapper::mapSampleDataTypeToRawDataType(OpenZGY::SampleDataType value)
{
  using InternalZGY::RawDataType;
  using OpenZGY::SampleDataType;
  switch (value) {
  case SampleDataType::int8:    return RawDataType::SignedInt8;
  case SampleDataType::int16:   return RawDataType::SignedInt16;
  case SampleDataType::float32: return RawDataType::Float32;
  case SampleDataType::unknown:
    throw Errors::ZgyUserError("SampleDataType::unknown is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for SampleDataType");
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
OpenZGY::UnitDimension
EnumMapper::mapRawHorizontalDimensionToUnitDimension(InternalZGY::RawHorizontalDimension value)
{
  using InternalZGY::RawHorizontalDimension;
  using OpenZGY::UnitDimension;
  switch (value) {
  default:
  case RawHorizontalDimension::Unknown:  return UnitDimension::unknown;
  case RawHorizontalDimension::Length:   return UnitDimension::length;
  case RawHorizontalDimension::ArcAngle: return UnitDimension::arcangle;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * The distinction between TWT and OWT, which is stored on the file, is lost.
 * To my knowledge there is no code that recognizes OWD in ZGY files anyway,
 * so it is better to not expose this to the API.
 */
OpenZGY::UnitDimension
EnumMapper::mapRawVerticalDimensionToUnitDimension(InternalZGY::RawVerticalDimension value)
{
  using InternalZGY::RawVerticalDimension;
  using OpenZGY::UnitDimension;
  switch (value) {
  default:
  case RawVerticalDimension::Unknown:    return UnitDimension::unknown;
  case RawVerticalDimension::Depth:      return UnitDimension::length;
  case RawVerticalDimension::SeismicTWT: return UnitDimension::time;
  case RawVerticalDimension::SeismicOWT: return UnitDimension::time;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * An explicit "unknown" is allowed, but passing e.g. "time" in this
 * context (which is supposedly horozontal) will raise an exception.
 */
InternalZGY::RawHorizontalDimension
EnumMapper::mapUnitDimensionToRawHorizontalDimension(OpenZGY::UnitDimension value)
{
  using InternalZGY::RawHorizontalDimension;
  using OpenZGY::UnitDimension;
  switch (value) {
  case UnitDimension::unknown:  return RawHorizontalDimension::Unknown;
  case UnitDimension::length:   return RawHorizontalDimension::Length;
  case UnitDimension::arcangle: return RawHorizontalDimension::ArcAngle;
  case UnitDimension::time:
    throw Errors::ZgyUserError("UnitDimension::time is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for UnitDimension");
  }
};

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * An explicit "unknown" is allowed, but passing e.g. "arcangle" in this
 * context (which is supposedly vertical) will raise an exception.
 */
InternalZGY::RawVerticalDimension
EnumMapper::mapUnitDimensionToRawVerticalDimension(OpenZGY::UnitDimension value)
{
  using InternalZGY::RawVerticalDimension;
  using OpenZGY::UnitDimension;
  switch (value) {
  case UnitDimension::unknown: return RawVerticalDimension::Unknown;
  case UnitDimension::time:    return RawVerticalDimension::SeismicTWT;
  case UnitDimension::length:  return RawVerticalDimension::Depth;
  case UnitDimension::arcangle:
    throw Errors::ZgyUserError("UnitDimension::arcangle is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for UnitDimension");
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
InternalZGY::LodAlgorithm
EnumMapper::mapDecimationTypeToLodAlgorithm(OpenZGY::DecimationType value)
{
  using InternalZGY::LodAlgorithm;
  using OpenZGY::DecimationType;
  switch (value) {
  case DecimationType::LowPass:          return LodAlgorithm::LowPass;
  case DecimationType::LowPassNew:       return LodAlgorithm::LowPassNew;
  case DecimationType::WeightedAverage:  return LodAlgorithm::WeightedAverage;
  case DecimationType::Average:          return LodAlgorithm::Average;
  case DecimationType::Median:           return LodAlgorithm::Median;
  case DecimationType::Minimum:          return LodAlgorithm::Minimum;
  case DecimationType::Maximum:          return LodAlgorithm::Maximum;
  //case DecimationType::MinMax:         return LodAlgorithm::MinMax;
  case DecimationType::Decimate:         return LodAlgorithm::Decimate;
  case DecimationType::DecimateSkipNaN:  return LodAlgorithm::DecimateSkipNaN;
  //case DecimationType::DecimateRandom: return LodAlgorithm::DecimateRandom;
  case DecimationType::AllZero:          return LodAlgorithm::AllZero;
  //case DecimationType::WhiteNoise:     return LodAlgorithm::WhiteNoise;
  case DecimationType::MostFrequent:     return LodAlgorithm::MostFrequent;
  case DecimationType::MostFrequentNon0: return LodAlgorithm::MostFrequentNon0;
  case DecimationType::AverageNon0:      return LodAlgorithm::AverageNon0;
  default:
    throw Errors::ZgyUserError("Invalid enum tag for DecimationType");
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
std::vector<InternalZGY::LodAlgorithm>
EnumMapper::mapDecimationTypeToLodAlgorithm(const std::vector<OpenZGY::DecimationType>& values)
{
  using InternalZGY::LodAlgorithm;
  using OpenZGY::DecimationType;
  if (!values.size())
    return std::vector<InternalZGY::LodAlgorithm>
      {LodAlgorithm::LowPass, LodAlgorithm::WeightedAverage};
  std::vector<InternalZGY::LodAlgorithm> result;
  for (DecimationType value : values)
    result.push_back(mapDecimationTypeToLodAlgorithm(value));
  return result;
}

}}
