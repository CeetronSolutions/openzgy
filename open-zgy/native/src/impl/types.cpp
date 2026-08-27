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

#include "types.h"
#include "enum.h"
#include "../exception.h"
#include <cstdint>
#include <limits>

namespace InternalZGY {
#if 0
}
#endif

RawDataTypeDetails::RawDataTypeDetails(RawDataType type)
{
  *this = get(type);
}

RawDataTypeDetails::RawDataTypeDetails(
        std::size_t size_in,
        double lowest_in,
        double highest_in,
        bool is_integer_in,
        bool is_signed_in)
  : size(size_in)
  , lowest(lowest_in)
  , highest(highest_in)
  , is_integer(is_integer_in)
  , is_signed(is_signed_in)
{
}

template<typename T>
RawDataTypeDetails::RawDataTypeDetails(T* dummy)
  : size(sizeof(T))
  , lowest    (std::numeric_limits<T>::lowest())
  , highest   (std::numeric_limits<T>::max())
  , is_integer(std::numeric_limits<T>::is_integer)
  , is_signed (std::numeric_limits<T>::is_signed)
{
}

const RawDataTypeDetails&
RawDataTypeDetails::get(RawDataType type)
{
  static RawDataTypeDetails SignedInt8_((std::int8_t*)nullptr);
  static RawDataTypeDetails UnsignedInt8_((std::uint8_t*)nullptr);
  static RawDataTypeDetails SignedInt16_((std::int16_t*)nullptr);
  static RawDataTypeDetails UnsignedInt16_((std::uint16_t*)nullptr);
  static RawDataTypeDetails SignedInt32_((std::int32_t*)nullptr);
  static RawDataTypeDetails UnsignedInt32_((std::uint32_t*)nullptr);
  static RawDataTypeDetails Float32_((float*)nullptr);
  static RawDataTypeDetails IbmFloat32_(4, -7.237e+75, +7.237e+75, false, true);
  switch (type) {
  case RawDataType::SignedInt8:    return SignedInt8_;
  case RawDataType::UnsignedInt8:  return UnsignedInt8_;
  case RawDataType::SignedInt16:   return SignedInt16_;
  case RawDataType::UnsignedInt16: return UnsignedInt16_;
  case RawDataType::SignedInt32:   return SignedInt32_;
  case RawDataType::UnsignedInt32: return UnsignedInt32_;
  case RawDataType::Float32:       return Float32_;
  case RawDataType::IbmFloat32:    return IbmFloat32_;
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

} // namespace
