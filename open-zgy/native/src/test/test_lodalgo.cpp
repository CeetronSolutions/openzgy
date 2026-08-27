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
#include "test_utils.h"
#include "../impl/lodalgo.h"
#include "../impl/databuffer.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <cstdint>
#include <cmath>
#include <omp.h>

using namespace InternalZGY;

namespace {
#if 0
}
#endif

class PushNumThreads
{
  int _old_count;
public:
  explicit PushNumThreads(int n)
  {
    _old_count = omp_get_max_threads();
    omp_set_num_threads(n);
    //std::cout << "Threadcount changed from " << _old_count
    //        << " to " << omp_get_max_threads() << std::endl;
  }
  ~PushNumThreads()
  {
    omp_set_num_threads(_old_count);
  }
};

// Use a brick size of 16x16x16 and a cube size of 2x2x2 bricks.
// The test data will be very monotonous, every 2x2x2 area
// effectively holding the same pattern. For most algorithms
// this means that the resulting LOD brick will consist of
// all identical values.

template <typename T>
static void
fillCube(DataBufferNd<T,3> *databuffer, double s0, double s1, double s2, double s3, double s4, double s5, double s6, double s7)
{
  T *cube = databuffer->data();
  const double ss[8] {s0, s1, s2, s3, s4, s5, s6, s7};
  T data[8] {};
  for (int ii=0; ii<8; ++ii) {
    if ((ss[ii] < 0 && !std::numeric_limits<T>::is_signed) || ss[ii] == -1) {
      if (std::numeric_limits<T>::has_quiet_NaN)
        data[ii] = std::numeric_limits<T>::quiet_NaN();
      else
        data[ii] = 0;
    }
    else {
      data[ii] = static_cast<T>(ss[ii]);
    }
  }

  const std::array<std::int64_t,3> size = databuffer->size3d();
  int pos = 0;
  for (int ii=0; ii<size[0]; ++ii)
    for (int jj=0; jj<size[1]; ++jj)
      for (int kk=0; kk<size[2]; ++kk)
        cube[pos++] = data[4*(ii&1)+2*(jj&1)+1*(kk&1)];
}

//static const int bricksize[3] = {16, 16, 16};
//static const int brickstride[3] = {16*16, 16, 1}; // i varies slowest
static const std::int64_t histogram[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

// The test cube will be filled with 1 2 3 4 5 6 and two NaNs (for float)
// or two zeroes (for int). The average is 21/6 = 3.5 for float data
// and 21/8 = 2.62 for integral data. Currently the integral case will
// round the result towards zero.
// For the LowPass and MinMax filters, the result is not expected to be
// the same value for all cells. So I'll just use a looser test to see
// thet the result is in the same range as the input.

static struct testdata
{
  const char *name;
  LodAlgorithm algorithm;
  double expectmin;
  double expectmax;
  int imin;
  int imax;
} tests[] = {
  { "LowPass",         LodAlgorithm::LowPass,         0.0, 9.0, 0, 9 },
  { "LowPassNew",      LodAlgorithm::LowPassNew,      0.0, 9.0, 0, 9 },
  { "WeightedAverage", LodAlgorithm::WeightedAverage, 3.4, 3.6, 2, 3 },
  { "Average",         LodAlgorithm::Average,         3.4, 3.6, 2, 3 },
  { "Median",          LodAlgorithm::Median,          3.4, 3.6, 2, 3 },
  { "Minimum",         LodAlgorithm::Minimum,         0.99, 1.01, 0, 0 },
  { "Maximum",         LodAlgorithm::Maximum,         5.99, 6.01, 6, 6 },
  //{ "MinMax",          LodAlgorithm::MinMax,          0.99, 6.01, 0, 6 },
  { "Decimate",        LodAlgorithm::Decimate,        -1,  -1, 0, 0 },
  { "DecimateSkipNaN", LodAlgorithm::DecimateSkipNaN, 1.99, 2.01, 0, 0 },
  //{ "DecimateRandom",  LodAlgorithm::DecimateRandom,  0.99, 6.01, 1, 6 },
  { "AllZero",         LodAlgorithm::AllZero,         0.00, 0.00, 0, 0 },
  //{ "WhiteNoise",      LodAlgorithm::WhiteNoise,      -1.0, 17.0, 0, 0 },
  { "MostFrequent",    LodAlgorithm::MostFrequent,    2.00, 2.00, 0, 0 },
  { "MostFrequentNon0",LodAlgorithm::MostFrequentNon0,2.00, 2.00, 2, 2 },
  { "AverageNon0",     LodAlgorithm::AverageNon0,     3.4, 3.6, 2, 3 },
};

template <typename T>
static void
testArrayTile()
{
  auto src = std::make_shared<DataBufferNd<T,3>>(index3_t{32,32,32});
  auto dst = std::make_shared<DataBufferNd<T,3>>(index3_t{16,16,16});

  for (unsigned int testno=0; testno<sizeof(tests)/sizeof(tests[0]); ++testno) {

    // A -1 passed to fillCube() means put a NaN in that place for float data
    // or a 0 for integral data. Similarly, a -1 in expectmin means we
    // expect to se a NaN value in the result.

    fillCube(src.get(), -1.0, 2.0, -1.0, 1.0, 3.0, 4.0, 5.0, 6.0);
    //fillCube(dst.get(), 0.0, 0.0, 0.0, 0.0, 99.0, 99.0, 99.0, 99.0);
    fillCube(dst.get(), 0.0, 0.0, 0.0, 0.0, 99.0, 98.0, 97.0, 96.0);

    if (verbose())
      printf("Testing %s %s\n", typeid(T).name(), tests[testno].name);

    createLod(dst, src, tests[testno].algorithm, 0, histogram, 16, 0.0, 16.0);

    int print_warnings{10};
    const T* dptr = dst->data();
    for (int ii=0; ii<16; ++ii)
      for (int jj=0; jj<16; ++jj)
        for (int kk=0; kk<16; ++kk) {
          double val = dptr[ii*16*16+jj*16+kk];
          double min = tests[testno].expectmin;
          double max = tests[testno].expectmax;
          if (std::numeric_limits<T>::is_integer) {
            min = tests[testno].imin;
            max = tests[testno].imax;
          }

          if (min >= 0 && val >= min && val <= max)
            ;
          else if (min<0 && !std::isfinite(val))
            ;
          else {
            if (--print_warnings >= 0) {
              printf("Mismatch at [%d,%d,%d] %s %s expect %g - %g got %g\n",
                     ii, jj, kk, typeid(T).name(),
                     tests[testno].name, min, max, val);
              TEST_CHECK(min == val);
            }
          }
        }
  }
}
#if 0
template void TestArrayBasic::TestArrayTile<float>();
template void TestArrayBasic::TestArrayTile<char>();
template void TestArrayBasic::TestArrayTile<unsigned char>();
template void TestArrayBasic::TestArrayTile<short>();
template void TestArrayBasic::TestArrayTile<unsigned short>();
template void TestArrayBasic::TestArrayTile<int>();
template void TestArrayBasic::TestArrayTile<unsigned int>();
#endif

static struct testdata2
{
  const char *name;
  LodAlgorithm algorithm;
} tests2[] = {
  { "LowPass",         LodAlgorithm::LowPass         },
  { "LowPassNew",      LodAlgorithm::LowPassNew      },
  { "WeightedAverage", LodAlgorithm::WeightedAverage },
  { "Average",         LodAlgorithm::Average         },
  { "Median",          LodAlgorithm::Median          },
  { "Minimum",         LodAlgorithm::Minimum         },
  { "Maximum",         LodAlgorithm::Maximum         },
  //{ "MinMax",          LodAlgorithm::MinMax          },
  { "Decimate",        LodAlgorithm::Decimate        },
  { "DecimateSkipNaN", LodAlgorithm::DecimateSkipNaN },
  { "MostFrequent",    LodAlgorithm::MostFrequent    },
  { "MostFrequentNon0",LodAlgorithm::MostFrequentNon0},
  { "AverageNon0",     LodAlgorithm::AverageNon0     },
  //{ "DecimateRandom",  LodAlgorithm::DecimateRandom },
  //{ "AllZero",         LodAlgorithm::AllZero        },
  //{ "WhiteNoise",      LodAlgorithm::WhiteNoise     },
};


/**
 * When the entire brick is filled with NaN values,
 * the LOD brick should be the same. Test this for all the
 * algorithms except those ignoring the input completely.
 * This particular test was added because it uncovered
 * a compiled bug in g++ 3.4.6 in Release mode.
 * So the test might not be terribly relevant today.
 */
static void
test_NaN()
{
  auto src = std::make_shared<DataBufferNd<float,3>>(index3_t{32,32,32});
  auto dst = std::make_shared<DataBufferNd<float,3>>(index3_t{16,16,16});
  src->fill(std::numeric_limits<float>::quiet_NaN());

  for (unsigned int testno=0; testno<sizeof(tests2)/sizeof(tests2[0]); ++testno) {
    createLod(dst, src, tests2[testno].algorithm, 0, histogram, 16, 0.0, 16.0);
    //printf("RESULT: %15s -> %g\n", tests2[testno].name, dst[0]);
    TEST_CHECK(!std::isfinite(dst->data()[0]));
  }
}

static void
test_int8()
{
  // Using a prime number will hit some corner cases because the
  // slices are very unlikely to divide evenly into the thread count.
  PushNumThreads pnt(7);
  testArrayTile<std::int8_t>();
}

static void
test_int16()
{
  PushNumThreads pnt(13);
  testArrayTile<std::int16_t>();
}

static void
test_float()
{
  PushNumThreads pnt(8);
  testArrayTile<float>();
}

/**
 * Code copied from GenLodImpl::decimateOneInputBrick().
 *
 * The caller has a multiple of 8 input data bricks to be decimated,
 * with each group of 8 ordered dim2 fastest, dim0 slowest. I.e. same
 * ordering of bricks as ordering of samples inside each brick.
 * And one output brick per group.
 *
 * Return the offset into the decimated output brick, in samples,
 * corresponding to the first sample in inputs[index].
 */
static std::int64_t
offsetInOutput(int index, const index3_t& bs)
{
  //const std::int64_t vbrick = index / 8;
  const std::int64_t subi = (index / 4) % 2;
  const std::int64_t subj = (index / 2) % 2;
  const std::int64_t subk = index % 2;
  return (subk * bs[2] / 2 +
          subj * bs[2] * bs[1] / 2 +
          subi * bs[2] * bs[1] * bs[0] / 2);
}

static void
test_scalar()
{
  const index3_t bs{8,16,32};
  const RawDataType dtype{RawDataType::SignedInt16};
  std::shared_ptr<DataBuffer> dst = DataBuffer::makeNewBuffer3d(bs, dtype);
  dst->fill(99);
  const auto get = [dst,bs](std::int64_t ii, std::int64_t jj, std::int64_t kk)
                   {
                     TEST_CHECK(ii >= 0 && ii < bs[0]);
                     TEST_CHECK(jj >= 0 && jj < bs[1]);
                     TEST_CHECK(kk >= 0 && kk < bs[2]);
                     return static_cast<std::int16_t*>(dst->voidData().get())
                       [ii*bs[1]*bs[2] + jj*bs[2] + kk];
                   };
  for (int ii=0; ii<8; ++ii)
  {
    std::shared_ptr<DataBuffer> src =
      DataBuffer::makeScalarBuffer3d(ii+1, bs, dtype);
    std::int64_t offset = offsetInOutput(ii, bs);
    createLod(dst, src, LodAlgorithm::Decimate, offset, nullptr, 0, 0, 0);
  }

  TEST_EQUAL(get(0,0,0),  1);
  TEST_EQUAL(get(0,0,1),  1);
  TEST_EQUAL(get(0,0,2),  1);
  TEST_EQUAL(get(0,0,16), 2);
  TEST_EQUAL(get(0,8,0),  3);
  TEST_EQUAL(get(0,8,16), 4);
  TEST_EQUAL(get(4,0,0),  5);
  TEST_EQUAL(get(4,0,16), 6);
  TEST_EQUAL(get(4,8,0),  7);
  TEST_EQUAL(get(4,8,16), 8);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("lodalgo.int8",            test_int8);
      register_test("lodalgo.int16",           test_int16);
      register_test("lodalgo.float",           test_float);
      register_test("lodalgo.int32",           testArrayTile<std::int32_t>);
      register_test("lodalgo.uint32",          testArrayTile<std::uint32_t>);
      register_test("lodalgo.uint16",          testArrayTile<std::uint16_t>);
      register_test("lodalgo.uint8",           testArrayTile<std::uint8_t>);
      register_test("lodalgo.nan",             test_NaN);
      register_test("lodalgo.scalar",          test_scalar);
    }
  } dummy;
} // namespace for registration
