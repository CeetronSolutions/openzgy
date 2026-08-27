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

#include "test_all.h"
#include "../impl/roundandclip.h"

#include <iostream>
#include <sstream>
#include <memory>

using namespace InternalZGY;

#define CPPUNIT_ASSERT(a) TEST_CHECK((a))
#define CPPUNIT_ASSERT_EQUAL(a,b) TEST_CHECK((a)==(b))
#define CPPUNIT_ASSERT_DOUBLES_EQUAL(a,b,eps) TEST_CHECK(std::abs((a)-(b)) <= (eps))

namespace {
#if 0
}
#endif

static void test_Floating()
{
  // Assign to floats should not change anything.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14f, RoundAndClip<float>(3.14), 0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14, RoundAndClip<double>(3.14), 0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14, (double)RoundAndClip<long double>(3.14), 0.0001);

  // any nan argument should be ignored.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14f, RoundAndClip<float>(3.14, 99), 0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14, RoundAndClip<double>(3.14, 99), 0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14, (double)RoundAndClip<long double>(3.14, 99), 0.0001);

  const double nan      = std::numeric_limits<double>::quiet_NaN();
  CPPUNIT_ASSERT(RoundAndClip<float>(nan, 99) != 99);
  CPPUNIT_ASSERT(RoundAndClip<double>(nan, 99) != 99);
  CPPUNIT_ASSERT(RoundAndClip<long double>(nan, 99) != 99);
}

template <typename T>
void testIntegral()
{
  // Input arguments
  const double nan      = std::numeric_limits<double>::quiet_NaN();
  const double abithigh = 42.1; // should be rounded down
  const double abitlow  = 41.9; // should be rounded up
  const double toohigh  = static_cast<double>(std::numeric_limits<T>::max()) + 42.0;
  const double toolow   = static_cast<double>(std::numeric_limits<T>::min()) - 42.0;
  const double veryhigh = std::numeric_limits<double>::infinity();
  const double verylow  = -std::numeric_limits<double>::infinity();
  const double nearhigh = static_cast<double>(std::numeric_limits<T>::max() - 1);
  const double nearlow  = static_cast<double>(std::numeric_limits<T>::min() + 1);

  // Not a const, so we hopefully avoid the "expression is const" later.
  bool is_signed = std::numeric_limits<T>::is_signed;

  // Expected results, expressed in the correct valuetype.
  const T fortytwo = static_cast<T>(42);
  // For unsigned types, -42 is of course outside the valid range
  // so the expected result is zero.
  const T minusfortytwo = is_signed ? static_cast<T>(-42) : 0;
  const T min = std::numeric_limits<T>::min();
  const T max = std::numeric_limits<T>::max();
  const T newnan = 99;

#if 0
  printf("\n\nRange is %llx to %llx\n", min, max);
  printf("%14g -> %llx\n", nan,      (long long)RoundAndClip<T>(nan));
  printf("%14g -> %llx\n", abithigh, (long long)RoundAndClip<T>(abithigh));
  printf("%14g -> %llx\n", abitlow,  (long long)RoundAndClip<T>(abitlow));
  printf("%14g -> %llx\n", -abithigh, (long long)RoundAndClip<T>(-abithigh));
  printf("%14g -> %llx\n", -abitlow,  (long long)RoundAndClip<T>(-abitlow));
  printf("%14g -> %llx\n", toohigh,  (long long)RoundAndClip<T>(toohigh));
  printf("%14g -> %llx\n", toolow,   (long long)RoundAndClip<T>(toolow));
  printf("%14g -> %llx\n", veryhigh, (long long)RoundAndClip<T>(veryhigh));
  printf("%14g -> %llx\n", verylow,  (long long)RoundAndClip<T>(verylow));
  printf("%14g -> %llx\n", nearhigh, (long long)RoundAndClip<T>(nearhigh));
  printf("%14g -> %llx\n", nearlow,  (long long)RoundAndClip<T>(nearlow));
#endif

  //CPPUNIT_ASSERT_EQUAL(0, RoundAndClip<T>(nan)); // undefined!
  CPPUNIT_ASSERT_EQUAL(fortytwo, RoundAndClip<T>(abithigh));
  CPPUNIT_ASSERT_EQUAL(fortytwo, RoundAndClip<T>(abitlow));
  CPPUNIT_ASSERT_EQUAL(max, RoundAndClip<T>(toohigh));
  CPPUNIT_ASSERT_EQUAL(min, RoundAndClip<T>(toolow));
  CPPUNIT_ASSERT_EQUAL(max, RoundAndClip<T>(veryhigh));
  CPPUNIT_ASSERT_EQUAL(min, RoundAndClip<T>(verylow));
  CPPUNIT_ASSERT_EQUAL(minusfortytwo, RoundAndClip<T>(-abithigh));
  CPPUNIT_ASSERT_EQUAL(minusfortytwo, RoundAndClip<T>(-abitlow));
  // The next one is tricky for unsigned int. RoundD2I doesn't handle it.
  // The test will fail for 64 bit integers, because RoundAndClip
  // takes a double argument which is not precise enough.
  // To avoid this, the expected and actual results are both
  // cast to double before cheching, to reduce the precision.
  CPPUNIT_ASSERT_EQUAL((double)nearhigh, (double)RoundAndClip<T>(nearhigh));
  CPPUNIT_ASSERT_EQUAL((double)nearlow, (double)RoundAndClip<T>(nearlow));

  // Try again, with the version that has explicit tests for NaN
  CPPUNIT_ASSERT_EQUAL(newnan, RoundAndClip<T>(nan, newnan));
  CPPUNIT_ASSERT_EQUAL(fortytwo, RoundAndClip<T>(abithigh, newnan));
  CPPUNIT_ASSERT_EQUAL(fortytwo, RoundAndClip<T>(abitlow, newnan));
  CPPUNIT_ASSERT_EQUAL(max, RoundAndClip<T>(toohigh, newnan));
  CPPUNIT_ASSERT_EQUAL(min, RoundAndClip<T>(toolow, newnan));
  CPPUNIT_ASSERT_EQUAL(max, RoundAndClip<T>(veryhigh, newnan));
  CPPUNIT_ASSERT_EQUAL(min, RoundAndClip<T>(verylow, newnan));
  CPPUNIT_ASSERT_EQUAL(minusfortytwo, RoundAndClip<T>(-abithigh, newnan));
  CPPUNIT_ASSERT_EQUAL(minusfortytwo, RoundAndClip<T>(-abitlow, newnan));
  CPPUNIT_ASSERT_EQUAL((double)nearhigh, (double)RoundAndClip<T>(nearhigh, newnan));
  CPPUNIT_ASSERT_EQUAL((double)nearlow, (double)RoundAndClip<T>(nearlow, newnan));
}


#if 0
#if defined(_MSC_VER) && !defined(_MANAGED)
#define NOINLINE __declspec(noinline) // before function declaration
#elif defined(__GNUC__)
#define NOINLINE __attribute__ ((noinline)) // after type, before decl?
#else
#define NOINLINE
#endif
#endif

template<typename T>
bool try_isfinite(T value)
{
  return IsFiniteT(value);
}

template<typename T>
bool try_isnan(T value)
{
  return IsNanT(value);
}

template <typename T>
void test_fin()
{
  T minf = -std::numeric_limits<T>::infinity();
  T pinf = +std::numeric_limits<T>::infinity();
  T qnan = std::numeric_limits<T>::quiet_NaN();
  T snan = std::numeric_limits<T>::signaling_NaN();
  T low1 = std::numeric_limits<T>::lowest();
  T high = std::numeric_limits<T>::max();
  T zero = 0;
  T pipi = static_cast<T>(3.14);

  CPPUNIT_ASSERT(!try_isfinite(minf));
  CPPUNIT_ASSERT(!try_isfinite(pinf));
  CPPUNIT_ASSERT(!try_isfinite(qnan));
  CPPUNIT_ASSERT(!try_isfinite(snan));
  CPPUNIT_ASSERT(try_isfinite(low1));
  CPPUNIT_ASSERT(try_isfinite(high));
  CPPUNIT_ASSERT(try_isfinite(zero));
  CPPUNIT_ASSERT(try_isfinite(pipi));

  CPPUNIT_ASSERT(!try_isnan(minf));
  CPPUNIT_ASSERT(!try_isnan(pinf));
  CPPUNIT_ASSERT(try_isnan(qnan));
  CPPUNIT_ASSERT(try_isnan(snan));
  CPPUNIT_ASSERT(!try_isnan(low1));
  CPPUNIT_ASSERT(!try_isnan(high));
  CPPUNIT_ASSERT(!try_isnan(zero));
  CPPUNIT_ASSERT(!try_isnan(pipi));

  CPPUNIT_ASSERT(!std::isfinite(minf));
  CPPUNIT_ASSERT(!std::isfinite(pinf));
  CPPUNIT_ASSERT(!std::isfinite(qnan));
  CPPUNIT_ASSERT(!std::isfinite(snan));
  CPPUNIT_ASSERT(std::isfinite(low1));
  CPPUNIT_ASSERT(std::isfinite(high));
  CPPUNIT_ASSERT(std::isfinite(zero));
  CPPUNIT_ASSERT(std::isfinite(pipi));

  CPPUNIT_ASSERT(!std::isnan(minf));
  CPPUNIT_ASSERT(!std::isnan(pinf));
  CPPUNIT_ASSERT(std::isnan(qnan));
  CPPUNIT_ASSERT(std::isnan(snan));
  CPPUNIT_ASSERT(!std::isnan(low1));
  CPPUNIT_ASSERT(!std::isnan(high));
  CPPUNIT_ASSERT(!std::isnan(zero));
  CPPUNIT_ASSERT(!std::isnan(pipi));
}

static void test_fin_int()
{
  CPPUNIT_ASSERT(try_isfinite(std::int64_t(42)));
  CPPUNIT_ASSERT(!try_isnan(std::int64_t(42)));
  CPPUNIT_ASSERT(try_isfinite(std::int32_t(42)));
  CPPUNIT_ASSERT(!try_isnan(std::int32_t(42)));
  CPPUNIT_ASSERT(try_isfinite(std::int16_t(42)));
  CPPUNIT_ASSERT(!try_isnan(std::int16_t(42)));
  CPPUNIT_ASSERT(try_isfinite(std::int8_t(42)));
  CPPUNIT_ASSERT(!try_isnan(std::int8_t(42)));
  CPPUNIT_ASSERT(try_isfinite(std::uint64_t(42)));
  CPPUNIT_ASSERT(!try_isnan(std::uint64_t(42)));
}

static void test_fin_float()
{
  test_fin<float>();
}

static void test_fin_double()
{
  test_fin<double>();
}

static void test_fin_ldouble()
{
  test_fin<long double>();
}

static void test_fin_const_float()
{
  test_fin<const float>();
}

static void test_fin_const_double()
{
  test_fin<const double>();
}

static void test_fin_const_ldouble()
{
  test_fin<const long double>();
}

static void test_char()
{
  testIntegral<char>();
}

static void test_schar()
{
  testIntegral<signed char>();
}

static void test_uchar()
{
  testIntegral<unsigned char>();
}

static void test_short()
{
  testIntegral<short>();
}

static void test_ushort()
{
  testIntegral<unsigned short>();
}

static void test_int()
{
  testIntegral<int>();
}

static void test_uint()
{
  testIntegral<unsigned int>();
}

static void test_longlong()
{
  testIntegral<long long int>();
}

static void test_ulonglong()
{
  testIntegral<unsigned long long int>();
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("roundandclip.float",      test_Floating);
      register_test("roundandclip.char",       test_char);
      register_test("roundandclip.schar",      test_schar);
      register_test("roundandclip.uchar",      test_uchar);
      register_test("roundandclip.short",      test_short);
      register_test("roundandclip.ushort",     test_ushort);
      register_test("roundandclip.int",        test_int);
      register_test("roundandclip.uint",       test_uint);
      register_test("roundandclip.longlong",   test_longlong);
      register_test("roundandclip.ulonglong",  test_ulonglong);
      register_test("roundandclip.f_int",      test_fin_int);
      register_test("roundandclip.f_float",    test_fin_float);
      register_test("roundandclip.f_double",   test_fin_double);
      register_test("roundandclip.f_ldouble",  test_fin_ldouble);
      register_test("roundandclip.f_cfloat",   test_fin_const_float);
      register_test("roundandclip.f_cdouble",  test_fin_const_double);
      register_test("roundandclip.f_cldouble", test_fin_const_ldouble);
    }
  } dummy;
} // namespace for registration

namespace ReadThisAssemblyCode
{
  // My bit-fiddling code.
  // Currently not handling long double, so that
  // one should fall back to library version.
  bool NewIsFinite(float value)
  {
    return InternalZGY::IsFiniteT(value);
  }

  bool NewIsFinite(double value)
  {
    return InternalZGY::IsFiniteT(value);
  }

  bool NewIsFinite(long double value)
  {
    return InternalZGY::IsFiniteT(value);
  }

  bool NewIsNan(float value)
  {
    return InternalZGY::IsNanT(value);
  }

  bool NewIsNan(double value)
  {
    return InternalZGY::IsNanT(value);
  }

  bool NewIsNan(long double value)
  {
    return InternalZGY::IsNanT(value);
  }

  // STL functions, not always optimal.
  bool StdIsFinite(float value)
  {
    return std::isfinite(value);
  }

  bool StdIsFinite(double value)
  {
    return std::isfinite(value);
  }

  bool StdIsFinite(long double value)
  {
    return std::isfinite(value);
  }

  bool StdIsNan(float value)
  {
    return std::isnan(value);
  }

  bool StdIsNan(double value)
  {
    return std::isnan(value);
  }

  bool StdIsNan(long double value)
  {
    return std::isnan(value);
  }

  // IsNan is in principle the same as x!=x
  bool PlainIsNan(float value)
  {
    return value != value;
  }

  bool PlainIsNan(double value)
  {
    return value != value;
  }

  bool PlainIsNan(long double value)
  {
    return value != value;
  }
}
