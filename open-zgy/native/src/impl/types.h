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

#pragma once

#include "enum.h"
#include <cstdint>

/**
 * \file types.h
 * \brief Types and enums only used internally.
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * Helper for enum RawDataType in enum.h:
 * convert a templated type into a RawDataType enum.
 */
template<typename T> struct RawDataTypeTraits {};
template<> struct RawDataTypeTraits<std::int8_t>   { static const RawDataType datatype = RawDataType::SignedInt8; };
template<> struct RawDataTypeTraits<std::uint8_t>  { static const RawDataType datatype = RawDataType::UnsignedInt8; };
template<> struct RawDataTypeTraits<std::int16_t>  { static const RawDataType datatype = RawDataType::SignedInt16; };
template<> struct RawDataTypeTraits<std::uint16_t> { static const RawDataType datatype = RawDataType::UnsignedInt16; };
template<> struct RawDataTypeTraits<std::int32_t>  { static const RawDataType datatype = RawDataType::SignedInt32; };
template<> struct RawDataTypeTraits<std::uint32_t> { static const RawDataType datatype = RawDataType::UnsignedInt32; };
template<> struct RawDataTypeTraits<float>         { static const RawDataType datatype = RawDataType::Float32; };

/**
 * Helper for enum RawDataType in enum.h:
 * Given an enum tag, return traits for the C++ type it represents.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class RawDataTypeDetails
{
public:
  std::size_t size;
  double lowest;
  double highest;
  bool is_integer;
  bool is_signed;
  explicit RawDataTypeDetails(RawDataType type);
private:
  RawDataTypeDetails(std::size_t size_in, double lowest_in, double highest_in, bool is_integer_in, bool is_signed_in);
  template<typename T> explicit RawDataTypeDetails(T* dummy);
  static const RawDataTypeDetails& get(RawDataType type);
};

} // namespace
