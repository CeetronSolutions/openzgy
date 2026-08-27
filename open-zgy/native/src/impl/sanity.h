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

#pragma once

#include "../declspec.h"
#include "enum.h"

#include <cstdint>

namespace InternalZGY {
#if 0
}
#endif

class OPENZGY_TEST_API SanityCheck
{
  SanityCheck() = delete; // Static class

private:
  static bool isPowerOfTwo(std::int64_t n);

  static bool willUnsignedMultiplyOverflow(
       std::uint64_t a,
       std::uint64_t b,
       std::uint64_t lim);

  static bool isNotPositiveOrWillMultiplyOverflow(
       std::int64_t a,
       std::int64_t b,
       std::uint64_t lim);

  static bool isNotPositiveOrWillMultiplyOverflow(
       const std::array<std::int64_t,3>& in,
       std::int64_t bytes_per_sample,
       std::uint64_t lim);

  static void checkValidBrickSize(
       const std::array<std::int64_t,3>& in,
       const std::int64_t bytes_per_sample);

  static void checkValidSurveySize(
       const std::array<std::int64_t,3>& in,
       const std::int64_t bytes_per_sample);

public:
  static void checkValidBrickAndSurveySize(
     const std::array<std::int64_t,3>& bricksize,
     const std::array<std::int64_t,3>& size,
     RawDataType dtype);
};

} // namespace
