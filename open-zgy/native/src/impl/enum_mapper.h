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

/**
 * \file impl/enum_mapper.cpp
 * \brief convert between internal and external data types.
 *
 * \details

 * Used by the API layer (api.cpp and writerargs.cpp) which puts it in
 * the OpenZGY namespace. But the code should not be visible to the
 * application, so it is moved down to OpenZGY::Impl. Note, not
 * InternalZGY because that is supposed to imply that all types and
 * methods in the API layer are off limits.

 * The plan was to have all code in below the OpenZGY namespace in
 * src/ instead of src/impl/. But on the other hand, header files in
 * src/ *.h are assumed to be accessible from the application. So I put
 * this in src/impl/.
 */

#pragma once

#include "../declspec.h"

#include "enum.h" // The internal types
#include <vector>

namespace OpenZGY {
  enum class SampleDataType; // api.h
  enum class UnitDimension;  // api.h
  enum class DecimationType; // api.h
}

namespace InternalZGY {
  enum class RawDataType;            // impl/enum.h
  enum class RawHorizontalDimension; // impl/enum.h
  enum class RawVerticalDimension;   // impl/enum.h
  enum class LodAlgorithm;           // impl/lodalgo.h
}

namespace OpenZGY { namespace Impl {
#if 0
}}
#endif

/**
 * \brief convert between internal and external data types.
 *
 * \details
 * Thread safe because the class only contains static methods with no state.
 *
 * \internal

 * Used by api.cpp (for the deprecated ZgyWriterArgs) and
 * writerargs.cpp (for ZgyWriterArgsV3). Also from the unit tests.
 * TODO write unit tests.
 */
class OPENZGY_TEST_API EnumMapper
{
  EnumMapper() = delete;
  EnumMapper(const EnumMapper&) = delete;
  EnumMapper& operator=(const EnumMapper&) = delete;
public:
  static SampleDataType mapRawDataTypeToSampleDataType(InternalZGY::RawDataType);
  static UnitDimension  mapRawHorizontalDimensionToUnitDimension(InternalZGY::RawHorizontalDimension);
  static UnitDimension  mapRawVerticalDimensionToUnitDimension(InternalZGY::RawVerticalDimension);
  static InternalZGY::RawDataType            mapSampleDataTypeToRawDataType(SampleDataType);
  static InternalZGY::RawHorizontalDimension mapUnitDimensionToRawHorizontalDimension(UnitDimension);
  static InternalZGY::RawVerticalDimension   mapUnitDimensionToRawVerticalDimension(UnitDimension);
  static InternalZGY::LodAlgorithm mapDecimationTypeToLodAlgorithm(DecimationType value);
  static std::vector<InternalZGY::LodAlgorithm> mapDecimationTypeToLodAlgorithm(const std::vector<DecimationType>& values);
};

}} // namespace OpenZGY::Impl
