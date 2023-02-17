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

//#define BIG_ENDIAN_ARCH 1

#include "test_all.h"
#include "../impl/structaccess.h"

#include <iostream>
//#include <memory>
#include <cstdint>

//using namespace OpenZGY;
using namespace InternalZGY;

namespace {
#if 0
}
#endif

/**
 * byteswapping single-byte variables should be a no-op.
 */
static void
test_byteswap_8()
{
  std::int8_t  a[2] { 42, -15 };
  std::uint8_t b    { 25 };
  byteswapAlways(a, 2);
  byteswapAlways(&b);
  TEST_CHECK(a[0] == 42);
  TEST_CHECK(a[1] == -15);
  TEST_CHECK(b == 25);
}

static void
test_byteswap_16()
{
  std::int16_t  a[2] { 0x1234, 0x5678 };
  std::uint16_t b    { 0x9ABC };
  byteswapAlways(a, 2);
  byteswapAlways(&b);
  TEST_CHECK(a[0] == 0x3412);
  TEST_CHECK(a[1] == 0x7856);
  TEST_CHECK(b == 0xBC9A);
}

static void
test_byteswap_32()
{
  std::int32_t  a { 0x12345678 };
  std::uint32_t b { 0x9ABC };
  byteswapAlways(&a);
  byteswapAlways(&b);
  TEST_CHECK(a == 0x78563412);
  TEST_CHECK(b == 0xBC9A0000);
}

static bool
compare64(std::int64_t a, std::uint64_t b)
{
  const std::int64_t signbit = static_cast<std::int64_t>(1) << 63;
  bool aneg = (a & signbit) != 0;
  bool bneg = (b & signbit) != 0;
  return (a & ~signbit) == static_cast<std::int64_t>(b & ~signbit) && aneg==bneg;
}

static void
test_byteswap_64()
{
  std::int64_t  a { 0x0BADC0DE };
  std::uint64_t b { 0xDEADBEEF0000F00D };
  byteswapAlways(&a);
  byteswapAlways(&b);
  const std::int64_t signbit = static_cast<std::int64_t>(1) << 63;
  TEST_CHECK((a & ~signbit) == (0xDEC0AD0B00000000 & ~signbit));
  TEST_CHECK((a & ~signbit) != 0);
  TEST_CHECK(b == 0x0DF00000EFBEADDE);
}

static void
test_big_little()
{
  std::int64_t  a[2] { 0x00DEAD0000BEEF00, 0x1122334455667788 };
  std::uint64_t b    { 0xDEADBEEF0000F00D };
  byteswapV1Long(a,2);
  byteswapV1Long(&b);
#if BIG_ENDIAN_ARCH
  TEST_CHECK(compare64(a[0], 0x00EFBE0000ADDE00));
  TEST_CHECK(compare64(a[1], 0x8877665544332211));
  TEST_CHECK(b    == 0x0DF00000EFBEADDE);
#else
  TEST_CHECK(compare64(a[0], 0x00BEEF0000DEAD00));
  TEST_CHECK(compare64(a[1], 0x5566778811223344));
  TEST_CHECK(b    == 0x0000F00DDEADBEEF);
#endif
  byteswapV1Long(a,2);
  byteswapV1Long(&b);
  TEST_CHECK(a[0] == 0x00DEAD0000BEEF00);
  TEST_CHECK(a[1] == 0x1122334455667788);
  TEST_CHECK(b    == 0xDEADBEEF0000F00D);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("struct.byteswap_8",       test_byteswap_8);
      register_test("struct.byteswap_16",      test_byteswap_16);
      register_test("struct.byteswap_32",      test_byteswap_32);
      register_test("struct.byteswap_64",      test_byteswap_64);
      register_test("struct.big_little",       test_big_little);
    }
  } dummy;
} // namespace for registration
