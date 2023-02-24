// Copyright 2017-2023, Schlumberger
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
#include "test_utils.h"
#include "../impl/edgebrick.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <cmath>

using namespace InternalZGY;

namespace {
#if 0
}
#endif

// Brick size 4 * 4 * 5 = 80
// ROI size   2 * 3 * 2 = 12
//
// First inline   Second inline
// 1  3  5  x      7  9 11  x
// 2  4  6  x      8 10 12  x
// x  x  x  x      x  x  x  x
// x  x  x  x      x  x  x  x
// x  x  x  x      x  x  x  x

static const std::int16_t my_values[]
{
   1,  2, 99, 99, 99,
   3,  4, 99, 99, 99,
   5,  6, 99, 99, 99,
  99, 99, 99, 99, 99,
   7,  8, 99, 99, 99,
   9, 10, 99, 99, 99,
  11, 12, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99,
  99, 99, 99, 99, 99
};

static const std::array<std::int64_t, 3> my_size{ 4,4,5 };
static const std::array<std::int64_t, 3> my_roi{ 2,3,2 };


static void check_contents(const EdgeBrick<const std::int16_t>& brick)
{
  std::int64_t count{ 0 };
#if 0
  for (const auto& it : brick) {
    if (!TEST_EQUAL(it, count + 1))
      break;
    ++count;
  }
#endif
  for (auto it = brick.begin(); it != brick.end(); ++it) {
    if (!TEST_EQUAL(*it, count + 1))
      break;
    ++count;
  }
  TEST_EQUAL(count, 12);
}

/**
 * forward_iterator must be:
 *  - default constructible,
 *  - copy constructible,
 *  - copy assignable,
 *  - destructible.
 */
static void test_construct()
{
  BrickIterator<std::int16_t> b1;
  BrickIterator<std::int16_t> b2(b1);
  b1 = b2;
}

/**
 * Can be compared for equivalence using the equality / inequality operators
 * when both iterator values iterate over the same underlying sequence.
 *
 * Can be incremented if in a dereferenceable state.
 * The result is either also dereferenceable or a past - the - end iterator.
 *
 * Can be dereferenced as an rvalue if in a dereferenceable state).
 *
 * Iterators that compare equal, keep comparing equal when both are increased.
 */
static void test_compare()
{
  BrickIterator<std::int16_t> b1;
  BrickIterator<std::int16_t> b2;
  TEST_CHECK(b1 == b2);
  TEST_CHECK(!(b1 != b2));

  EdgeBrick<const std::int16_t> brick(my_values, my_size, my_roi);
  check_contents(brick);
  auto b3 = brick.begin();
  auto b4 = brick.begin();
  TEST_CHECK(b3 == b4 && !(b3 != b4));
  ++b3;
  TEST_CHECK(b3 != b4 && !(b3 == b4));
  TEST_EQUAL(*b3, 2);
  TEST_EQUAL(*b4, 1);
  b4++;
  TEST_CHECK(b3 == b4 && !(b3 != b4));
  TEST_CHECK(b4 != brick.end());
  for (int ii = 0; ii < 11; ++ii) {
    if (!TEST_CHECK(b4 != brick.end()))
      break;
    ++b4;
  }
#if 0
  if (verbose()) {
    std::cout << "Should be at end: (" << b4.pos_[0] << "," << b4.pos_[1] << "," << b4.pos_[2] << ") offset = " << (b4.data_ptr_ - b4.data_) << std::endl;
    auto b5 = brick.end();
    std::cout << "Actual end:       (" << b5.pos_[0] << "," << b5.pos_[1] << "," << b5.pos_[2] << ") offset = " << (b5.data_ptr_ - b5.data_) << std::endl;
  }
#endif
  TEST_CHECK(b4 == brick.end());
}

/**
 * For mutable iterators:
 * Can be dereferenced as an lvalue if in a dereferenceable state.
 */
static void test_update()
{
  std::vector<std::int16_t> data(my_values, my_values + 80);
  EdgeBrick<std::int16_t> brick(data.data(), my_size, my_roi);
  for (auto& it : brick) {
    ++it;
  }
  TEST_EQUAL(*brick.begin(), 2);
  TEST_EQUAL(data[0], 2);
  TEST_EQUAL(data[1], 3);
  TEST_EQUAL(data[2], 99);

  for (auto it = brick.begin(); it != brick.end(); ++it) {
    *it = *it * -1;
  }
  TEST_EQUAL(*brick.begin(), -2);
  TEST_EQUAL(data[0], -2);
  TEST_EQUAL(data[1], -3);
  TEST_EQUAL(data[2], 99);
}

/**
 * In C++11, Lvalues must be swappable.
 */
static void test_swap()
{
  EdgeBrick<const std::int16_t> brick(my_values, my_size, my_roi);
  auto a = brick.begin();
  auto b = brick.begin();
  ++b;
  b++;
  *b++;
  TEST_EQUAL(*a, 1);
  TEST_EQUAL(*b, 4);
  std::swap(a, b);
  TEST_EQUAL(*a, 4);
  TEST_EQUAL(*b, 1);
}

/**
 * Degenerate case: roi == size.
 */
static void test_trivial()
{
  std::vector<std::int16_t> data;
  for (int ii = 0; ii < my_roi[0] * my_roi[1] * my_roi[2]; ++ii)
    data.push_back(ii + 1);
  EdgeBrick<const std::int16_t> brick(data.data(), my_roi, my_roi);
  check_contents(brick);
}

static void test_misc()
{
  // Size zero in some dimensions, positive in others.
  // This is silly, but should behave as an empty set.
  EdgeBrick<const std::int16_t> brick(nullptr, std::array<std::int64_t, 3>{1, 0, 3}, std::array<std::int64_t, 3>{1, 0, 2});
  TEST_CHECK(brick.begin() == brick.end());
}

/**
* Test with int8 instead of int16
*/
static void test_int8()
{
  std::vector<std::int8_t> data;
  for (int ii = 0; ii < 80; ++ii)
    data.push_back((std::int8_t)my_values[ii]);
  EdgeBrick<const std::int8_t> brick(data.data(), my_size, my_roi);
  EdgeBrick<const std::int16_t> check(my_values, my_size, my_roi);
  auto b = brick.begin();
  auto c = check.begin();
  for (; b != brick.end() && c != check.end(); ++b, ++c)
    if (!TEST_EQUAL(*b, *c))
      break;
  TEST_CHECK(b == brick.end());
  TEST_CHECK(c == check.end());
}

/**
 * Test with float instead of int16
 */
static void test_float()
{
  std::vector<float> data;
  for (int ii = 0; ii < 80; ++ii)
    data.push_back(my_values[ii]);
  EdgeBrick<float> brick(data.data(), my_size, my_roi);
  EdgeBrick<const std::int16_t> check(my_values, my_size, my_roi);
  auto b = brick.begin();
  auto c = check.begin();
  for (; b != brick.end() && c != check.end(); ++b, ++c)
    if (!TEST_EQUAL_FLOAT(*b, *c, 0.0001f))
      break;
  TEST_CHECK(b == brick.end());
  TEST_CHECK(c == check.end());
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("edgebrick_construct",     test_construct);
      register_test("edgebrick_compare",       test_compare);
      register_test("edgebrick_update",        test_update);
      register_test("edgebrick_swap",          test_swap);
      register_test("edgebrick_trivial",       test_trivial);
      register_test("edgebrick_misc",          test_misc);
      register_test("edgebrick_int8",          test_int8);
      register_test("edgebrick_float",         test_float);
    }
  } dummy;
} // namespace for registration
