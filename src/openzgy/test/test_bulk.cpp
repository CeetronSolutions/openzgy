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

#include "test_all.h"
#include "../impl/bulk.h"
#include "../impl/databuffer.h"

#include <iostream>
#include <memory>

using namespace InternalZGY;

namespace Test {
  class TestBulk {
  public:
    static void test_padding_const();
    static void test_padding_edge();
    template<typename T>
    static void test_padding_samples();
  };
}

namespace {
  template<typename T>
  static std::shared_ptr<DataBufferNd<T,3>>
  make_testdata(const std::array<std::int64_t,3> size,
                const std::array<std::int64_t,3> used)
  {
    auto data = std::make_shared<DataBufferNd<T,3>>(size);
    data->fill(0);
    T *raw = data->data();
    for (std::int64_t ii=0; ii<used[0]; ++ii)
      for (std::int64_t jj=0; jj<used[1]; ++jj)
        for (std::int64_t kk=0; kk<used[2]; ++kk)
          raw[ii*size[1]*size[2] + jj*size[2] + kk] =
            static_cast<T>(
              (kk == used[2]-1) ? 2 : // bottom edge
              (jj == used[1]-1) ? 3 : // last crossline
              (ii == used[0]-1) ? 4 : // last inline
              1);
    return data;
  }
  template<typename T>
  static std::array<int,7>
  count_testdata(const std::shared_ptr<DataBufferNd<T,3>>& data)
  {
    std::array<int,7> counts{0};
    T *raw = data->data();
    T *end = raw + data->allocsize();
    for (const T *ptr = raw; ptr < end; ++ptr)
      ++counts[*ptr < 0 || *ptr >= 6 ? 6 : int(*ptr)];
    if (false)
      std::cout << "\n"
                << "Counts: unwritten     " << counts[0] << "\n"
                << "        inner region  " << counts[1] << "\n"
                << "        bottom slice  " << counts[2] << "\n"
                << "        last xl slice " << counts[3] << "\n"
                << "        last il slice " << counts[4] << "\n"
                << "        const padded  " << counts[5] << "\n"
                << "        garbage       " << counts[6] << "\n"
        ;
    return counts;
  }
} // namespace

void
Test::TestBulk::test_padding_const()
{
  const std::array<std::int64_t,3> size{32, 64, 128};
  const std::array<std::int64_t,3> used{7, 13, 17};
  auto data = make_testdata<float>(size, used);

  ZgyInternalBulk::_setPaddingToConst(data, used, 5, 0);
  ZgyInternalBulk::_setPaddingToConst(data, used, 5, 1);
  ZgyInternalBulk::_setPaddingToConst(data, used, 5, 2);

  // Instead of computing the expected value of every sample
  // (which would replicate most of the logic I am testing)
  // I will check only the result histogram.
  std::array<int,7> counts = count_testdata<float>(data);
  TEST_CHECK(counts[0] == 0); // From initial fill, not overwritten.
  TEST_CHECK(counts[1] == (used[0]-1)*(used[1]-1)*(used[2]-1)); // inner
  TEST_CHECK(counts[2] == used[0]*used[1]);     // One "bottom" slice.
  TEST_CHECK(counts[3] == used[0]*(used[2]-1)); // One crossline, except bottom.
  TEST_CHECK(counts[4] == (used[1]-1)*(used[2]-1)); // One inline, except...
  TEST_CHECK(counts[5] == size[0]*size[1]*size[2] - used[0]*used[1]*used[2]);
  TEST_CHECK(counts[6] == 0); // Garbage value we never wrote.
}

void
Test::TestBulk::test_padding_edge()
{
  const std::array<std::int64_t,3> size{32, 64, 128};
  const std::array<std::int64_t,3> used{7, 13, 17};
  const std::array<std::int64_t,3> padd{8, 16, 20};
  auto data = make_testdata<float>(size, used);

  ZgyInternalBulk::_setPaddingToEdge(data, used, 4, 0);
  ZgyInternalBulk::_setPaddingToEdge(data, used, 4, 1);
  ZgyInternalBulk::_setPaddingToEdge(data, used, 4, 2);

  std::array<int,7> counts = count_testdata<float>(data);
  TEST_CHECK(counts[0] == size[0]*size[1]*size[2] - padd[0]*padd[1]*padd[2]);
  TEST_CHECK(counts[1] == (used[0]-1)*(used[1]-1)*(used[2]-1)); // inner
  TEST_CHECK(counts[2] == 4*padd[0]*padd[1]);     // 4 "bottom" slice.
  TEST_CHECK(counts[3] == 4*padd[0]*(used[2]-1)); // 4 crossline, except bottom.
  TEST_CHECK(counts[4] == 2*(used[1]-1)*(used[2]-1)); // 2 inline, except...
  TEST_CHECK(counts[5] == 0);
  TEST_CHECK(counts[6] == 0); // Garbage value we never wrote.
}

template<typename T>
void
Test::TestBulk::test_padding_samples()
{
  const std::array<std::int64_t,3> size{32, 64, 128};
  const std::array<std::int64_t,3> used{7, 13, 17};
  const std::array<std::int64_t,3> padd{8, 16, 20};
  auto data = make_testdata<T>(size, used);

  ZgyInternalBulk::_setPaddingSamples(data, used, 5, compressor_t());

  // Expected result identical to test_padding_edge except all the
  // "unwritten" samples are now "const padded".
  std::array<int,7> counts = count_testdata<T>(data);
  TEST_CHECK(counts[5] == size[0]*size[1]*size[2] - padd[0]*padd[1]*padd[2]);
  TEST_CHECK(counts[1] == (used[0]-1)*(used[1]-1)*(used[2]-1)); // inner
  TEST_CHECK(counts[2] == 4*padd[0]*padd[1]);     // 4 "bottom" slice.
  TEST_CHECK(counts[3] == 4*padd[0]*(used[2]-1)); // 4 crossline, except bottom.
  TEST_CHECK(counts[4] == 2*(used[1]-1)*(used[2]-1)); // 2 inline, except...
  TEST_CHECK(counts[0] == 0);
  TEST_CHECK(counts[6] == 0); // Garbage value we never wrote.
}

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("bulk.padding_const",      Test::TestBulk::test_padding_const);
      register_test("bulk.padding_edge",       Test::TestBulk::test_padding_edge);
      register_test("bulk.padding_samples<f>", Test::TestBulk::test_padding_samples<float>);
      register_test("bulk.padding_samples<h>", Test::TestBulk::test_padding_samples<std::int16_t>);
      register_test("bulk.padding_samples<b>", Test::TestBulk::test_padding_samples<std::int8_t>);
    }
  } dummy;
} // namespace for registration
