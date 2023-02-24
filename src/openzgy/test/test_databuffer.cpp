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
#include "../api.h"
#include "../exception.h"
#include "../impl/databuffer.h"

#include <iostream>
#include <iomanip>
#include <memory>
#include <omp.h>

using namespace OpenZGY;
using namespace InternalZGY;
using Test_Utils::must_throw;

template<typename T>
static void
test_databuffer_construct()
{
  typedef std::array<int64_t,3> size3i_t;
  {
    //** Buffer represents a constant value. No allocated data.
    DataBufferNd<T,3> buffer(static_cast<T>(3.5), size3i_t{2, 5, 10});
    TEST_CHECK(buffer.isScalar());
    TEST_CHECK(buffer.scalarValue() == static_cast<T>(3.5));
    TEST_CHECK(buffer.scalarAsDouble() == buffer.scalarValue());
    TEST_CHECK(buffer.data() != nullptr);
    TEST_CHECK(buffer.allocsize() == 1);
    TEST_CHECK(buffer.totalsize() == 2*5*10);
    TEST_CHECK(buffer.itemsize() == sizeof(T));
    TEST_CHECK(buffer.size3d()[0] == 2);
    TEST_CHECK(buffer.size3d()[1] == 5);
    TEST_CHECK(buffer.size3d()[2] == 10);
    TEST_CHECK(buffer.stride3d()[0] == 0);
    TEST_CHECK(buffer.stride3d()[1] == 0);
    TEST_CHECK(buffer.stride3d()[2] == 0);
  }

  {
    //** Allocate a new buffer and take ownership of it.
    DataBufferNd<T,3> buffer(size3i_t{64, 50, 100});
    buffer.fill(0);
    TEST_CHECK(!buffer.isScalar());
    TEST_CHECK(buffer.data() != nullptr);
    TEST_CHECK(buffer.data() == buffer.voidData().get());
    const DataBufferNd<T,3>& cbuffer(buffer);
    TEST_CHECK(buffer.data() == cbuffer.data());
    TEST_CHECK(buffer.data() == cbuffer.voidData().get());
    TEST_CHECK(buffer.allocsize() == 64*50*100);
    TEST_CHECK(buffer.totalsize() == 64*50*100);
    TEST_CHECK(buffer.itemsize() == sizeof(T));
    TEST_CHECK(buffer.size3d()[0] == 64);
    TEST_CHECK(buffer.size3d()[1] == 50);
    TEST_CHECK(buffer.size3d()[2] == 100);
    TEST_CHECK(buffer.stride3d()[0] == 100*50);
    TEST_CHECK(buffer.stride3d()[1] == 100);
    TEST_CHECK(buffer.stride3d()[2] == 1);
  }

  {
    std::shared_ptr<T> mydata(new T[1000], std::default_delete<T[]>());
    DataBufferNd<T,3> buffer(mydata, size3i_t{10, 10, 10});
    TEST_CHECK(!buffer.isScalar());
    TEST_CHECK(buffer.data() == mydata.get());
    TEST_CHECK(buffer.allocsize() == 10*10*10);
    TEST_CHECK(buffer.totalsize() == 10*10*10);
    TEST_CHECK(buffer.itemsize() == sizeof(T));
  }

  {
    // Use the convenience function to make a scalar of the correct type.
    // This is a bit backwards since the code first converts T to a
    // RawDataType, so the function has to convert it back.
    RawDataType dtype =  RawDataTypeTraits<T>::datatype;
    std::shared_ptr<DataBuffer> buffer =
      DataBuffer::makeScalarBuffer3d(3.14, size3i_t{64, 50, 100}, dtype);
    TEST_CHECK(buffer->isScalar());
    TEST_CHECK(buffer->scalarAsDouble() == static_cast<double>(static_cast<T>(3.14)));
  }

  {
    // Use the convenience function to make a buffer of the correct type.
    RawDataType dtype =  RawDataTypeTraits<T>::datatype;
    const size3i_t size{3, 5, 7};
    const std::int64_t count = size[0]*size[1]*size[2];
    std::shared_ptr<T> data(new T[count], std::default_delete<T[]>());
    std::fill(data.get(), data.get() + count, (T)1);
    data.get()[42] = 42;
    std::shared_ptr<DataBuffer> buffer =
      DataBuffer::makeDataBuffer3d(data, count*sizeof(T), size, dtype);
    TEST_CHECK(!buffer->isScalar());
    // These should actually be separate tests, for range(), clear(), clone()
    auto r = buffer->range();
    TEST_CHECK(r.first == 1);
    TEST_CHECK(r.second == 42);
    auto cloned = buffer->clone();
    TEST_CHECK(buffer->voidData() != cloned->voidData()); // Deep clone
    r = buffer->range();
    TEST_CHECK(r.first == 1);   // Check data was copied.
    TEST_CHECK(r.second == 42);
    TEST_CHECK(data.use_count() != 1); // Because our "data" still in scope.
    buffer->clear();
    TEST_CHECK(data.use_count() == 1);  // The databuffer should have released it.
  }
}

static void
test_databuffer_range()
{
  typedef std::array<int64_t,3> size3i_t;
  const size3i_t size{10, 12, 19};
  const size3i_t usedsize{3, 5, 7};
  DataBufferNd<float,3> buffer(size);
  std::fill(buffer.data(), buffer.data() + (10*12*19), -999.25f);
  int value{0};
  for (int ii=0; ii<usedsize[0]; ++ii)
    for (int jj=0; jj<usedsize[1]; ++jj)
      for (int kk=0; kk<usedsize[2]; ++kk)
        buffer.data()[ii*size[1]*size[2] + jj*size[2] + kk] = (float)++value;
  auto r1 = buffer.range(nullptr); // entire buffer
  auto r2 = buffer.range(usedsize.data());  // valid part
  TEST_EQUAL(r1.first, -999.25);
  TEST_EQUAL(r1.second, 3*5*7);
  TEST_EQUAL(r2.first, 1);
  TEST_EQUAL(r2.second, 3*5*7);
}

template<typename T>
static void
test_databuffer_scale()
{
  const std::array<double,2> factors{2.0, 5.0};
  const std::array<std::int64_t,3> size{12, 50, 100};
  typedef DataBufferNd<T,3> intbuffer_t;
  typedef DataBufferNd<float,3> floatbuffer_t;
  typedef std::numeric_limits<T> limits;
  RawDataType datatype =  RawDataTypeTraits<T>::datatype;
  {
    // Constant value
    auto buffer = std::make_shared<intbuffer_t>((T)7, size);
    auto scaled1 = buffer ? buffer->scaleToFloat(factors) : nullptr;
    auto scaled2 = scaled1 ? scaled1->scaleToFloat(factors) : nullptr;
    auto scaled = std::dynamic_pointer_cast<floatbuffer_t>(scaled1);
    TEST_CHECK(buffer != nullptr && buffer->datatype() == datatype);
    if (!limits::is_integer) {
      TEST_CHECK(scaled1 == nullptr);
    }
    else {
      TEST_CHECK(scaled != nullptr);
      TEST_CHECK(buffer != scaled1);
      TEST_CHECK(!scaled2);
      if (scaled != nullptr) {
        TEST_CHECK(scaled->datatype() == RawDataType::Float32);
        TEST_CHECK(scaled->isScalar());
        TEST_CHECK(scaled->scalarValue() == 7 * 2 + 5);
        TEST_CHECK(scaled->allocsize() == 1);
        TEST_CHECK(scaled->totalsize() == 12*50*100);
        TEST_CHECK(scaled->safesize() == size);
      }

      // Convert back to integer
      auto check1 = scaled1->scaleToStorage(factors, datatype);
      auto check = std::dynamic_pointer_cast<intbuffer_t>(check1);
      TEST_CHECK(check1 != nullptr);
      TEST_CHECK(check != nullptr);
      if (check != nullptr) {
        TEST_CHECK(check->datatype() == datatype);
        TEST_CHECK(check->isScalar());
        TEST_CHECK(check->scalarValue() == 7);
        TEST_CHECK(check->allocsize() == 1);
        TEST_CHECK(check->totalsize() == 12*50*100);
        TEST_CHECK(check->safesize() == size);
      }
    }
  }

  {
    //** Regular buffer
    auto buffer = std::make_shared<intbuffer_t>(size);
    buffer->fill(0);
    buffer->data()[0] = 0;
    buffer->data()[1] = 1;
    buffer->data()[2] = 7;
    auto scaled1 = buffer ? buffer->scaleToFloat(factors) : nullptr;
    auto scaled2 = scaled1 ? scaled1->scaleToFloat(factors) : nullptr;
    auto scaled = std::dynamic_pointer_cast<floatbuffer_t>(scaled1);
    TEST_CHECK(buffer != nullptr && buffer->datatype() == datatype);
    if (!limits::is_integer) {
      TEST_CHECK(scaled1 == nullptr);
    }
    else {
      TEST_CHECK(scaled != nullptr);
      TEST_CHECK(buffer != scaled1);
      TEST_CHECK(!scaled2);
      if (scaled != nullptr) {
        TEST_CHECK(scaled->datatype() == RawDataType::Float32);
        TEST_CHECK(!scaled->isScalar());
        TEST_CHECK(scaled->allocsize() == scaled->totalsize());
        TEST_CHECK(scaled->totalsize() == 12*50*100);
        TEST_CHECK(scaled->itemsize() == sizeof(float));
        TEST_CHECK(scaled->safesize() == size);
        TEST_CHECK(scaled->data()[0] == 5);
        TEST_CHECK(scaled->data()[1] == 7);
        TEST_CHECK(scaled->data()[2] == 7 * 2 + 5);

        // Add some more floats to check rounding
        scaled->data()[3] = 18.9f;
        scaled->data()[4] = 19.1f;
        if (limits::is_integer && limits::is_signed) {
          scaled->data()[5] = -16.9f;
          scaled->data()[6] = -17.1f;
        }
      }

      // Convert back to integer
      auto check1 = scaled1->scaleToStorage(factors, datatype);
      auto check = std::dynamic_pointer_cast<intbuffer_t>(check1);
      TEST_CHECK(check1 != nullptr);
      TEST_CHECK(check != nullptr);
      if (check != nullptr) {
        TEST_CHECK(check->datatype() == datatype);
        TEST_CHECK(!check->isScalar());
        TEST_CHECK(check->allocsize() == check->allocsize());
        TEST_CHECK(check->totalsize() == 12*50*100);
        TEST_CHECK(check->itemsize() == sizeof(T));
        TEST_CHECK(check->safesize() == size);
        TEST_CHECK(check->data()[0] == 0);
        TEST_CHECK(check->data()[1] == 1);
        TEST_CHECK(check->data()[2] == 7);
        TEST_CHECK(check->data()[3] == 7);
        TEST_CHECK(check->data()[4] == 7);
        if (limits::is_integer && limits::is_signed) {
          TEST_CHECK(check->data()[5] == -11);
          TEST_CHECK(check->data()[6] == -11);
        }
      }
    }
  }

  {
    // Negative test: Convert an integral buffer to storage should throw.
    // Convert float buffer to float is a no-op and should return nullptr.
    auto buffer = std::make_shared<intbuffer_t>((T)42, size);
    TEST_CHECK(buffer != nullptr);
    if (buffer != nullptr) {
      if (limits::is_integer) {
        try {
          auto scaled = buffer->scaleToStorage(factors, datatype);
          TEST_CHECK(false && "Did not get exception in bad convert");
        }
        catch (const OpenZGY::Errors::ZgyInternalError&)
        {
        }
      }
      else {
        auto scaled = buffer->scaleToStorage(factors, datatype);
        TEST_CHECK(scaled == nullptr);
      }
    }
  }
}

template<typename T>
static std::int64_t
countequal(float value, const T* data, std::int64_t n)
{
  std::int64_t result = 0;
  while (n--)
    if (*data++ == value)
      ++result;
  return result;
}

template<typename T>
static void
test_databuffer_copyfrom()
{
  //Testing:
  //copyFrom(const DataBuffer* src,
  //       const std::int64_t *srcorig,
  //       const std::int64_t *dstorig,
  //       const std::int64_t *cpyorig,
  //       const std::int64_t *cpysize)

  // A user has requested reading data from survey's origin and with
  // size (100, 90, 110). The buffer is allocated by us.
  DataBufferNd<T,3> target(std::array<std::int64_t,3>{100,90,110});
  target.fill(0);
  std::array<std::int64_t,3> targetorigin{0, 0, 0};

  // 64^3 samples worth of 99 has been read from (0,0,64) in the survey
  // which means it doesn't go to the start of the user's buffer and
  // the user isn't going to make use of more than 46 samples vertically.
  DataBufferNd<T,3> brick(std::array<std::int64_t,3>{64,64,64});
  std::fill(brick.data(), brick.data()+brick.allocsize(), (T)99);
  std::int64_t orig1[3]{0,0,64};
  target.copyFrom(&brick, orig1, targetorigin.data(), nullptr, nullptr);

  TEST_CHECK(target.data()[0] == 0); // data[0,0,0] was not touched.
  TEST_CHECK(target.data()[64] == 99); // data[0,0,64] should have been.
  TEST_CHECK(countequal(99,target.data(),target.totalsize()) == 64*64*46);

  // 100 samples implicitly constant 3 at survey address (7,9,13).
  // User will get all of these.
  DataBufferNd<T,3> constbrick(3, std::array<std::int64_t,3>{2,5,10});
  std::int64_t orig2[3]{7,9,13};
  target.copyFrom(&constbrick, orig2, targetorigin.data(), nullptr, nullptr);

  // As above, but the constant is now the answer to life, the universe, and
  // everything. And the additional constraint is now set to only copy the
  // first sample. Note that the constraint is specified survey-relative.
  DataBufferNd<T,3> constbrick2(42, std::array<std::int64_t,3>{2,5,10});
  std::int64_t justone[3]{1,1,1};
  target.copyFrom(&constbrick2, orig2, targetorigin.data(),
                  orig2, justone);

  TEST_CHECK(target.data()[7*90*110 + 9*110 + 13] == 42);
  TEST_CHECK(target.data()[7*90*110 + 9*110 + 14] == 3);
  TEST_CHECK(countequal(42,target.data(),target.totalsize()) == 1);
  TEST_CHECK(countequal(3,target.data(),target.totalsize()) == 2*5*10-1);
}

/**
 * Used as a performance test only. There is no assert on a specific
 * maximum time. Instead, look at the test log which gives the elapsed
 * time of every test.
 */
template<typename T, int LOOPS>
static void
test_databuffer_copyfromST()
{
  DataBufferNd<T,3> target(std::array<std::int64_t,3>{64,64,64});
  DataBufferNd<T,3> brick1(std::array<std::int64_t,3>{64,64,64});
  DataBufferNd<T,3> brick2(std::array<std::int64_t,3>{64,64,64});
  std::fill(target.data(), target.data()+target.allocsize(), (T)0);
  std::fill(brick1.data(), brick1.data()+brick1.allocsize(), (T)99);
  std::fill(brick2.data(), brick2.data()+brick2.allocsize(), (T)42);
  const std::int64_t origin[3]{0, 0, 0};

  for (int ii=0; ii<LOOPS; ii += 2) {
    target.copyFrom(&brick1, origin, origin, nullptr, nullptr);
    target.copyFrom(&brick2, origin, origin, nullptr, nullptr);
  }
}

/**
 * Used as a performance test only. There is no assert on a specific
 * maximum time. Instead, look at the test log which gives the elapsed
 * time of every test.
 */
template<typename T, int LOOPS,int THREADS>
static void
test_databuffer_copyfromMT()
{
  DataBufferNd<T,3> target(std::array<std::int64_t,3>{64,64,64});
  DataBufferNd<T,3> brick1(std::array<std::int64_t,3>{64,64,64});
  DataBufferNd<T,3> brick2(std::array<std::int64_t,3>{64,64,64});
  std::fill(brick1.data(), brick1.data()+brick1.allocsize(), (T)99);
  std::fill(brick2.data(), brick2.data()+brick2.allocsize(), (T)42);
  const std::int64_t origin[3]{0, 0, 0};
  int actual_threads = 0;
#pragma omp parallel num_threads(THREADS)
  {
    DataBufferNd<T,3> target(std::array<std::int64_t,3>{64,64,64});
    if (omp_get_thread_num() == 0) {
      actual_threads = omp_get_num_threads();
    }
#pragma omp for
    for (int ii=0; ii<LOOPS; ii += 2) {
      target.copyFrom(&brick1, origin, origin, nullptr, nullptr);
      target.copyFrom(&brick2, origin, origin, nullptr, nullptr);
    }
  }
  TEST_EQUAL(actual_threads, THREADS);
}

static void
test_databuffer_copytoscalar()
{
  const std::array<std::int64_t,3>zero{0, 0, 0};
  const std::array<std::int64_t,3>srcsize{100, 90, 110};
  const std::array<std::int64_t,3>dstsize{100, 90, 120};
  DataBufferNd<float,3> source(srcsize);
  DataBufferNd<float,3> target(42, dstsize);
  source.fill(0);
  source.data()[0] = 1;
  source.data()[100*90*110 - 1] = 2;
  must_throw("copy to scalar buffer", [&](){
    target.copyFrom(&source, zero.data(), zero.data(), nullptr, nullptr);
  });
}

template<typename T>
static void
test_databuffer_tostring()
{
  DataBufferNd<T,3> data(std::array<std::int64_t,3>{11,22,33});
  TEST_CHECK(data.toString().find("(11, 22, 33)") != std::string::npos);
  TEST_CHECK(data.toString().find("(726, 33, 1)") != std::string::npos);
}

template<typename T>
static void
test_databuffer_issame()
{
  std::array<std::int64_t,3> size{32,64,128};
  std::array<std::int64_t,3> used{7,13,17};
  DataBufferNd<T,3> data(size);
  data.fill(42);
  data.data()[111] = 43;
  TEST_CHECK(!data.isAllSame(size.data()));
  TEST_CHECK(data.isAllSame(used.data()));
  data.data()[0] = 1;
  TEST_CHECK(!data.isAllSame(used.data()));
}

static void
test_databuffer_isnan()
{
  std::array<std::int64_t,3> size{32,64,128};
  std::array<std::int64_t,3> used{7,13,17};
  DataBufferNd<float,3> data(size);
  data.fill(std::numeric_limits<float>::quiet_NaN());
  data.data()[111] = 43;
  TEST_CHECK(!data.isAllSame(size.data()));
  TEST_CHECK(data.isAllSame(used.data()));
  data.data()[0] = 1;
  TEST_CHECK(!data.isAllSame(used.data()));
}

template<typename T>
static void
test_databuffer_slice()
{
  typedef std::array<std::int64_t,3> ndsize_t;
  typedef DataBufferNd<T,3> buffer_t;
  buffer_t buff(ndsize_t{32,64,128});
  buff.fill(0);
  TEST_CHECK(buff.voidData().use_count() == 2);
  buff.fill(1);

  std::shared_ptr<buffer_t> part;
  part = std::dynamic_pointer_cast<buffer_t>(buff.slice1(0, 7, 3));
  if (verbose()) {
    std::cout << "\n"
              << "Slice input:  " << buff.toString() << "\n"
              << "Slice output: " << part->toString() << "\n";
  }
  TEST_CHECK(part->stride3d() == buff.stride3d());
  TEST_CHECK(part->size3d() == (ndsize_t{3,64,128}));
  TEST_CHECK(part->allocsize() == part->totalsize());
  TEST_CHECK(buff.voidData().use_count() == 3);
  part->data()[0] = 42;
  TEST_CHECK(buff.data()[7*64*128] == 42);

  buffer_t scalar(42, ndsize_t{32,64,128});
  TEST_CHECK(scalar.isScalar());
}

/**
 * Instead of registering 7 tests * 7 value types, just register
 * test_databuffer_all with 7 valuetypes. It reduces clutter in the log.
 */
template<typename T>
static void
test_databuffer_all()
{
  test_databuffer_construct<T> ();
  test_databuffer_scale<T>();
  test_databuffer_copyfrom<T>();
  test_databuffer_tostring<T>();
  test_databuffer_issame<T>();
  test_databuffer_isnan(); // only makes sense for float
  test_databuffer_slice<T>();
}

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("databuffer.range",  test_databuffer_range);
      register_test("databuffer.copytoscalar", test_databuffer_copytoscalar);
#if 1
      register_test("databuffer.all<f>", test_databuffer_all<float>);
      register_test("databuffer.all<i>", test_databuffer_all<std::int32_t>);
      register_test("databuffer.all<I>", test_databuffer_all<std::uint32_t>);
      register_test("databuffer.all<h>", test_databuffer_all<std::int16_t>);
      register_test("databuffer.all<H>", test_databuffer_all<std::uint16_t>);
      register_test("databuffer.all<b>", test_databuffer_all<std::int8_t>);
      register_test("databuffer.all<B>", test_databuffer_all<std::uint8_t>);
#endif
#if 0
      register_test("databuffer.construct<f>", test_databuffer_construct<float>);
      register_test("databuffer.scale<h>",     test_databuffer_scale<std::int16_t>);
      register_test("databuffer.copyfrom<f>",  test_databuffer_copyfrom<float>);
      register_test("databuffer.tostring<f>",  test_databuffer_tostring<float>);
      register_test("databuffer.issame<f>",    test_databuffer_issame<float>);
      register_test("databuffer.isnan",        test_databuffer_isnan);
      register_test("databuffer.slice<h>",     test_databuffer_slice<std::int16_t>);
#endif
      // Raise the iteration count at least by a factor of 100 to get
      // meaningful timing results. Do *not* do this in a valgrind
      // run; that might take 20 minutes per test.
      register_test("databuffer.copyfromST",
                    test_databuffer_copyfromST<float,200>);
      register_test("databuffer.copyfrom02",
                    test_databuffer_copyfromMT<float,1000,2>);
      register_test("databuffer.copyfrom04",
                    test_databuffer_copyfromMT<float,1000,4>);
      register_test("databuffer.copyfrom08",
                    test_databuffer_copyfromMT<float,1000,8>);
      register_test("databuffer.copyfrom16",
                    test_databuffer_copyfromMT<float,1000,16>);
      register_test("databuffer.copyfrom24",
                    test_databuffer_copyfromMT<float,1000,24>);
      register_test("databuffer.copyfrom32",
                    test_databuffer_copyfromMT<float,1000,32>);
      register_test("databuffer.copyfrom48",
                    test_databuffer_copyfromMT<float,1000,48>);
      register_test("databuffer.copyfrom64",
                    test_databuffer_copyfromMT<float,1000,64>);
      register_test("databuffer.copyfrom96",
                    test_databuffer_copyfromMT<float,1000,96>);
    }
  } dummy;
}
