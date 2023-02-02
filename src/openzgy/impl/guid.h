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

#pragma once

#include <cstdint>
#include <cstddef>
#include <array>
#include <string>
#include <ostream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file guid.h
 * \brief Simplified %GUID handling. Only big endian, random number guids.
 * \details
 * See doc/uuid.md
 */

/**
 * \brief Simplified %GUID handling. Only big endian, random number guids.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class GUID
{
public:
  typedef std::array<std::uint8_t,16> guid_bytes_t;

public:
  GUID();
  explicit GUID(const guid_bytes_t& in);
  explicit GUID(std::nullptr_t);
  std::string toString() const;
  void copyTo(std::uint8_t *ptr, std::int64_t len);

private:
  static guid_bytes_t generate();
  static std::string format(const guid_bytes_t& guid);

private:
  guid_bytes_t _data;
};

namespace Formatters
{
  extern std::ostream& operator<<(std::ostream& os, const ::InternalZGY::GUID& guid);
}

} // namespace

