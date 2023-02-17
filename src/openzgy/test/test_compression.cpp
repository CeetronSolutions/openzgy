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
#include "../impl/compression.h"
#include "../exception.h"
#include "../api.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace InternalZGY;
using Test_Utils::must_throw;

namespace Test {
  class MockCompressPlugin;
  class MockCompressInstaller;
  class TestCompression {
  public:
    static void test_registration();
    static void test_called();
    static void test_zfp();
  };
}

namespace {
  template<typename T>
  static rawdata_t
  makeRawData(int n)
  {
    std::shared_ptr<T> ptr(new T[n], std::default_delete<T[]>());
    memset(ptr.get(), 0, sizeof(T)*n);
    return rawdata_t{ptr,sizeof(T)*n};
  }
}

class Test::MockCompressPlugin
{
public:
  static int compress_count;
  static int decompress_count;

public:
  static compressor_t getCompressor(const std::vector<std::string>&)
  {
    return [](const rawdata_t& data, const index3_t& shape) {
             return compress(data, shape);
           };
  }

  /**
   * An aray of exactly 256 floats can be compressed
   * at an unmentionably low quality. Anything else
   * will be rejected.
   */
  static rawdata_t compress(const rawdata_t& data, const index3_t& shape)
  {
    TEST_CHECK(data.second == shape[0]*shape[1]*shape[2]*(int)sizeof(float));
    ++compress_count;
    if (data.second == 256*sizeof(float)) {
      std::shared_ptr<char> ptr(new char[64], std::default_delete<char[]>());
      memset(ptr.get(), 0, 64);
#ifdef _WIN32
      strcpy_s(ptr.get(), 64, "Mock");
#else
      strcpy(ptr.get(), "Mock");
#endif
      return rawdata_t{ptr, 64};
    }
    return rawdata_t{nullptr, 0};
#if 0
    // If I need a real compressor, here is a starting point.
    // Probably overkill for a test; just exercise ZFP instead.
    const float *src = static_cast<const float*>(data.first.get());
    const std::int64_t size = data.second / sizeof(float);
    if (size < 64)
      return rawdata_t{nullptr, 0};
    const std::int64_t outsize = (size/2+1) + 2;
    float *dst = new float[outsize];
    rawdata_t result{std::shared_ptr<float>(dst), outsize};
    dst[0] = 42; // my magic number
    dst[1] = size; // decompress must be able to strip trailing garbage
    //Oops! not handling odd input size correctly.
    for (std::int64_t ii = 0, jj = 2; ii < size; ii += 2, jj += 1)
      dst[jj] = (src[ii] + src[ii+1]) / 2;
    return result;
#endif
  }

  /**
   * A compressed data block of 64 bytes or more,
   * starting with the magic "Mock", can be decompressed
   * albeit with zero snr. The result will be 256 zeros.
   * Anything else will be rejected.
   */
  static rawdata_t decompress(const rawdata_t& cdata,
                              const BrickStatus& status,
                              const index3_t& shape)
  {
    ++decompress_count;
    const char *ptr = static_cast<const char*>(cdata.first.get());
    if (ptr && cdata.second >= 64 && 0==strncmp(ptr, "Mock", 4)) {
      return makeRawData<float>(256);
    }
    return rawdata_t{nullptr, 0};
  }
};

int Test::MockCompressPlugin::compress_count = 0;
int Test::MockCompressPlugin::decompress_count = 0;

class Test::MockCompressInstaller
{
public:
  MockCompressInstaller()
  {
    CompressFactoryImpl::registerCompressor("Mock", &MockCompressPlugin::getCompressor);
    CompressFactoryImpl::registerDecompressor("Mock", &MockCompressPlugin::decompress);
    MockCompressPlugin::compress_count = 0;
    MockCompressPlugin::decompress_count = 0;
  }
  ~MockCompressInstaller()
  {
    CompressFactoryImpl::registerCompressor("Mock", nullptr);
    CompressFactoryImpl::registerDecompressor("Mock", nullptr);
  }
};

void
Test::TestCompression::test_registration()
{
  {
    MockCompressInstaller mi;
    std::vector<std::string> list1 = CompressFactoryImpl::knownCompressors();
    std::vector<std::string> list2 = CompressFactoryImpl::knownDecompressors();
    TEST_CHECK(std::find(list1.begin(), list1.end(), std::string("Mock")) != list1.end());
    TEST_CHECK(std::find(list2.begin(), list2.end(), std::string("Mock")) != list2.end());
    auto cfn = CompressFactoryImpl::getCompressor("Mock", std::vector<std::string>());
    TEST_CHECK(bool(cfn));
  }
  // "mi" out of scope, so the registration ought to be removed.
  std::vector<std::string> list1 = CompressFactoryImpl::knownCompressors();
  std::vector<std::string> list2 = CompressFactoryImpl::knownDecompressors();
  TEST_CHECK(std::find(list1.begin(), list1.end(), std::string("Mock")) == list1.end());
  TEST_CHECK(std::find(list2.begin(), list2.end(), std::string("Mock")) == list2.end());
  must_throw("not recognized",
             [&](){
               CompressFactoryImpl::getCompressor("Mock", std::vector<std::string>());
             });
}

void
Test::TestCompression::test_called()
{
  MockCompressInstaller mi;
  auto cfn = CompressFactoryImpl::getCompressor("Mock", std::vector<std::string>());
  TEST_CHECK(bool(cfn));
  TEST_CHECK(Test::MockCompressPlugin::compress_count == 0);
  TEST_CHECK(Test::MockCompressPlugin::decompress_count == 0);

  // Small block, cannot compress.
  rawdata_t r1 = cfn(makeRawData<float>(2), index3_t{1,1,2});
  TEST_CHECK(Test::MockCompressPlugin::compress_count == 1);
  TEST_CHECK(Test::MockCompressPlugin::decompress_count == 0);
  TEST_CHECK(!r1.first);

  // 256 samples, should be able to compress.
  rawdata_t r2 = cfn(makeRawData<float>(256), index3_t{2,4,32});
  TEST_CHECK(Test::MockCompressPlugin::compress_count == 2);
  TEST_CHECK(Test::MockCompressPlugin::decompress_count == 0);
  TEST_CHECK(!!r2.first);

  // Decompress what we successfully compressed.
  CompressFactoryImpl::decompress(r2, BrickStatus::Compressed, index3_t{2,4,32});
  TEST_CHECK(Test::MockCompressPlugin::compress_count == 2);
  TEST_CHECK(Test::MockCompressPlugin::decompress_count == 1);
}

void
Test::TestCompression::test_zfp()
{
  auto cfn = CompressFactoryImpl::getCompressor
    ("ZFP", std::vector<std::string>{"snr", "35"});
  TEST_CHECK(bool(cfn));

  rawdata_t idata = makeRawData<float>(8*16*32);
  index3_t idata_size{8,16,32};
  float* idata_ptr =
    const_cast<float*>(static_cast<const float*>(idata.first.get()));
  idata_ptr[5] = 42;
  idata_ptr[6] = 40;
  rawdata_t cdata = cfn(idata, idata_size);
  if (verbose())
    std::cout << "ZFP compressed from " << idata.second
              << " to " << cdata.second << " bytes.\n";
  TEST_CHECK(!!cdata.first);
  TEST_CHECK(cdata.second == 64);

  // Shrink to fit the cdata buffer and add some trailing noise after it.
  // Otherwise the compressor might secretly be returning more data than
  // it says it claims, and the dcompressir might use them.
  auto tmp = std::shared_ptr<char>(new char[cdata.second + 4096], std::default_delete<char[]>());
  memcpy(tmp.get(), cdata.first.get(), cdata.second);
  memset(tmp.get()+cdata.second, 0xE5, 4096);
  cdata.first = tmp;

  rawdata_t ddata = CompressFactoryImpl::decompress
    (cdata, BrickStatus::Compressed, idata_size);

  const float *check = static_cast<const float*>(ddata.first.get());
  if (verbose()) {
    std::cout << "Round trip";
    for (int ii=0; ii<8; ++ii)
      std::cout << " " << check[ii];
    std::cout << "\n";
  }
  TEST_CHECK(!!ddata.first);
  TEST_CHECK(ddata.second == idata.second);
  TEST_CHECK(check[0] == 0);
  TEST_CHECK(check[1] == 0);
  TEST_CHECK(check[2] == 0);
  TEST_CHECK(check[3] == 0);
  TEST_CHECK(std::abs(check[4] -  0.0f) < 1.0f);
  TEST_CHECK(std::abs(check[5] - 42.0f) < 1.0f);
  TEST_CHECK(std::abs(check[6] - 40.0f) < 1.0f);
  TEST_CHECK(std::abs(check[7] -  0.0f) < 1.0f);
}

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("compress.registration",         Test::TestCompression::test_registration);
      register_test("compress.called",               Test::TestCompression::test_called);
#ifdef HAVE_ZFP
      register_test("compress.zfp",                  Test::TestCompression::test_zfp);
#endif
    }
  } dummy;
} // namespace for registration
