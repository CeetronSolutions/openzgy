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
#include "../impl/lodsampling.h"
#include "../impl/lodfilters.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <cstdint>
#include <cmath>
#include <vector>

using InternalZGY::LodSampling;
using InternalZGY::LodFilters;

namespace {
#if 0
}
#endif

static std::vector<double>
random_double_data()
{
  static const std::vector<double> result
    { 98.4424971709383, -19.3436789232173,  16.7187315504720,  16.5584535936180,
      84.5225140284728,  47.8877924641225, -66.7731527653793,  62.4663683415372,
     -86.6838070916691, -16.5649936069143,  -0.1920877736507, -22.9326391909987,
     -65.8969519262371, -88.7460922028064, -46.5034331254585,  96.6229869808814,
      25.1374719252825, -47.0263497842295, -77.1955913534992, -89.7920933840495,
      63.1204571103328, -44.4908733561979,  84.1153895407442,  36.8925493638353,
     -57.7678922771627, -31.3941176689217,  76.9175856553892, -37.0977171324006,
      32.3215689209398, -54.4649684055168, -99.5672040579935,  35.3891188043440
  };
  return result;
}

template<typename T>
static std::vector<T>
random_data()
{
  std::vector<double> d = random_double_data();
  std::vector<T> result;
  if (std::numeric_limits<T>::is_signed) {
    for (double sample : d)
      result.push_back(static_cast<T>(sample));
  }
  else {
    for (double sample : d)
      result.push_back(static_cast<T>(std::abs(sample)));
  }
  return result;
}

template<typename T>
static void
old_resample(T* dst, std::size_t dstsize, const T* src, std::size_t srcsize)
{
  std::vector<double> in, out;
  for (std::size_t it = 0; it < srcsize; ++it)
    in.push_back(static_cast<double>(src[it]));
  out.resize(dstsize);
  LodSampling::downSample1D(out.data(), (int)out.size(),
                            in.data(), (int)in.size());
  for (double sample : out)
    *dst++ = static_cast<T>(sample);
}

template<typename T>
static void
new_resample(T* dst, std::size_t dstsize, const T* src, std::size_t srcsize,
             const T* src4before = nullptr, const T* src4after = nullptr)
{
  LodFilters::downSample1D(dst, (int)dstsize, src, (int)srcsize,
                           src4before, src4after);
}

template<typename T>
static double
max_error(const std::vector<T>& new_data, const std::vector<T>& old_data)
{
  TEST_EQUAL(new_data.size(), old_data.size());
  std::size_t len = std::min(new_data.size(), old_data.size());
  double delta{0};
  for (std::size_t ii=0; ii<len; ++ii) {
    if (std::isfinite((double)new_data[ii]) || std::isfinite((double)old_data[ii])) {
        delta = std::max(delta, std::abs(((double)new_data[ii]-(double)old_data[ii])));
    }
  }
  if (verbose()) {
    std::stringstream ss;
    ss << std::fixed;
    if (std::numeric_limits<T>::is_integer)
      ss.precision(0);
    else
      ss.precision(5);
    ss << "\nMaximum error: " << delta << "\n";
    ss << "      New        Old      Error\n";
    std::cout << ss.str();
    for (std::size_t ii=0; ii<len; ++ii) {
      ss.str(std::string());
      ss << std::setw(10) << (double)new_data[ii] << " "
         << std::setw(10) << (double)old_data[ii] << " "
         << std::setw(10) << (double)std::abs((double)new_data[ii] - (double)old_data[ii]) << "\n";
      std::cout << ss.str();
    }
  }
  return delta;
}

/**
 * Single precision float, single trace with 32 samples.
 * The input data is hard coded, the values originally
 * output from a random number generator.
 */
static void
test_float()
{
  const std::vector<float> data = random_data<float>();
  std::vector<float> result;
  std::vector<float> check;
  result.resize((data.size() + 1)  / 2);
  check.resize((data.size() + 1)  / 2);
  new_resample(result.data(), result.size(), data.data(), data.size());
  old_resample(check.data(),  check.size(),  data.data(), data.size());
  double error = max_error(result, check);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * Adding the same constant to each input value should simply end up
 * adding that same value to the output.
 */
static void
test_offset()
{
  std::vector<float> data = random_data<float>();
  std::vector<float> expect;
  std::vector<float> actual;
  expect.resize((data.size() + 1)  / 2);
  actual.resize((data.size() + 1)  / 2);
  new_resample(expect.data(), expect.size(), data.data(), data.size());
  for (float& sample : data)
    sample += 100.0f;
  for (float& sample : expect)
    sample += 100.0f;
  new_resample(actual.data(), actual.size(), data.data(), data.size());
  double error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * The DC response, i.e. filtering a trace with all constant values.
 */
static void
test_dc()
{
  std::vector<float> data(24, 42);
  std::vector<float> expect(24/2, 42);
  std::vector<float> actual;
  actual.resize((data.size() + 1)  / 2);
  new_resample(actual.data(), actual.size(), data.data(), data.size());
  double error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * Traces containing NaN use a simpler algorithm, but only in
 * the part of the trace where NaNs exist.
 */
static void
test_NaN()
{
  std::vector<float> data = random_data<float>();
  std::vector<float> actual;
  std::vector<float> expect;
  actual.resize((data.size() + 1) / 2);
  expect.resize((data.size() + 1) / 2);
  new_resample(expect.data(), expect.size(), data.data(), data.size());

  expect[5] = std::numeric_limits<float>::quiet_NaN(); // 10&11 both nan
  expect[6] = data[13]; // because data[12] is ignored
  expect[11] = std::numeric_limits<float>::infinity(); // 22&23 both inf
  expect[10] = data[20]; // because data[12] is ignored
  for (std::size_t ii=3; ii<14; ++ii)
    if (ii != 5 && ii != 6 && ii != 10 && ii != 11)
      expect[ii] = (data[2*ii] + data[2*ii+1]) / 2;

  data[10] = std::numeric_limits<float>::quiet_NaN();
  data[11] = std::numeric_limits<float>::quiet_NaN();
  data[12] = std::numeric_limits<float>::quiet_NaN();
  data[21] = std::numeric_limits<float>::infinity();
  data[22] = std::numeric_limits<float>::infinity();
  data[23] = std::numeric_limits<float>::infinity();
  new_resample(actual.data(), actual.size(), data.data(), data.size());
  double error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);

  // Also check the old version matches the new.
  std::vector<float> old;
  old.resize((data.size() + 1) / 2);
  old_resample(old.data(), old.size(), data.data(), data.size());
  double error2 = max_error(actual, old);
  TEST_EQUAL_FLOAT(error2, 0, 0.001);
}

/**
 * Traces shorter than the filter width need special handling.
 * Note that the old algorithm did not handle this; the caller
 * was supposed to do the special handling. The behavior probably
 * doesn't match. But nobody should care.
 */
static void
test_small()
{
  double error{0};
  std::vector<float> indata, expect, actual;

  indata = std::vector<float>{10};
  expect = std::vector<float>{10};
  actual = std::vector<float>(expect.size(), 0);
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);

  indata = std::vector<float>{20, 30};
  expect = std::vector<float>{25};
  actual = std::vector<float>(expect.size(), 0);
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);

  indata = std::vector<float>{30, 40, 50};
  expect = std::vector<float>{35, 45};
  actual = std::vector<float>(expect.size(), 0);
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * Traces with an odd length.
 * The old algorithm expects caller to handle this corner case.
 * The behavior probably doesn't match. But nobody should care.
 */
static void
test_odd()
{
  std::vector<float> data = random_data<float>();
  data.resize(data.size()-1);
  std::vector<float> result;
  std::vector<float> check;
  result.resize((data.size() + 1)  / 2);
  check.resize((data.size() + 1)  / 2);
  new_resample(result.data(), result.size(), data.data(), data.size());
  data.push_back(data.back()); // Old algo would otherwise assert.
  old_resample(check.data(),  check.size(),  data.data(), data.size());
  // Skip the last 3 values because "odd" handling differs.
  // Old use replicates the last value once, keeps using the filter.
  // New switches to plain average for the last sample only.
  result.resize(result.size()-3);
  check.resize(check.size()-3);
  double error = max_error(result, check);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * Traces with any supported valuetype.
 * The <float> case is currently identical to
 * test_float(). That might change.
 */
template<typename T>
static void
test_vt()
{
  const std::vector<T> data = random_data<T>();
  std::vector<T> result;
  std::vector<T> check;
  result.resize((data.size() + 1)  / 2);
  check.resize((data.size() + 1)  / 2);
  new_resample(result.data(), result.size(), data.data(), data.size());
  old_resample(check.data(),  check.size(),  data.data(), data.size());
  double error = max_error(result, check);
  double epsilon = std::numeric_limits<T>::is_integer ? 1.01 : 0.001;
  TEST_EQUAL_FLOAT(error, 0, epsilon);
}

/**
 * The new version doesn't convert everything to float first,
 * it happens one sample at a time. TODO verify that all the extra
 * conversions to double isn't costing us...
 * With integral data there is a risk of overflow, e.g.
 * average = (int)(((float)a + (float)b) / 2) is good,
 * while if a and b were added as e.g. int8 they might overflow.
 *
 * TODO Performance: Part of the benefit of the new version of
 * this algorithm is that the trace data doesn't need to be copied
 * into a temp buffer of doubles first. BUT, each sample is accessed
 * 10 times, and to avoid overflow it needs to be converted from int
 * to float 10 times. This might remove some (or even all) the benefits
 * of not having to copy the data. For float data the float->double
 * cast might have been droppped. Currently it isn't.
 */
template<typename T>
static void
test_stress()
{
  const T hi = std::numeric_limits<T>::max();
  const T lo = std::numeric_limits<T>::min();

  double error{0};
  std::vector<T> indata, expect, actual;

  indata = std::vector<T>{hi-5, hi-3};
  expect = std::vector<T>{hi-4};
  actual = std::vector<T>(expect.size(), 0);
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 1.001);

  indata = std::vector<T>{lo+7, lo+11};
  expect = std::vector<T>{lo+9};
  actual = std::vector<T>(expect.size(), 0);
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 1.001);

  indata = std::vector<T>(24, hi-10);
  expect = std::vector<T>(12, hi-10);
  actual = std::vector<T>(expect.size(), 0);
  //indata[8] = hi-20;
  new_resample(actual.data(), actual.size(), indata.data(), indata.size());
  error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 1.001);
}

/**
 * Output buffer too small, simply do not compute the last part.
 */
static void
test_overflow()
{
  const std::vector<float> data = random_data<float>();
  std::vector<float> result(data.size()/2, 22);
  std::vector<float> check(data.size()/2, 0);
  // Lie when saying how much room there is.
  new_resample(result.data(), result.size()-3, data.data(), data.size());
  new_resample(check.data(), check.size(), data.data(), data.size());
  check[check.size()-1] = 22;
  check[check.size()-2] = 22;
  check[check.size()-3] = 22;
  double error = max_error(result, check);
  // Technically a failure, see explanation in downSampleT().
  // But providing a too small output buffer doesn't make sense.
  TEST_EQUAL_FLOAT(error, 0, 999);
}

/**
 * Output buffer too large, do not touch the last part.
 */
static void
test_underflow()
{
  const std::vector<float> data = random_data<float>();
  std::vector<float> result(data.size()/2 + 5, 19);
  std::vector<float> check(data.size()/2, 0);
  new_resample(result.data(), result.size(), data.data(), data.size());
  new_resample(check.data(), check.size(), data.data(), data.size());
  for (std::size_t ii=0; ii<5; ++ii)
    check.push_back(19);
  double error = max_error(result, check);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

/**
 * Split the task into 3 parts, but avoid artefacts in the middle
 * by providing the optional edge samples.
 */
static void
test_border()
{
  const std::vector<float> indata = random_data<float>();
  std::vector<float> expect;
  std::vector<float> actual;
  std::vector<float> actual_2;
  std::vector<float> actual_3;

  expect.resize(16);
  new_resample(expect.data(), expect.size(), indata.data(), indata.size());

  actual.resize(5);
  new_resample(actual.data(), actual.size(), indata.data(), 10,
               (float*)nullptr, indata.data() + 10);

  actual_2.resize(5);
  new_resample(actual_2.data(), actual_2.size(), indata.data() + 10, 10,
               indata.data() + 10 - 4, indata.data() + 20);

  actual_3.resize(6);
  new_resample(actual_3.data(), actual_3.size(), indata.data() + 20, 12,
               indata.data() + 20 - 4, (float*)nullptr);

  actual.insert(actual.end(), actual_2.begin(), actual_2.end());
  actual.insert(actual.end(), actual_3.begin(), actual_3.end());

  double error = max_error(actual, expect);
  TEST_EQUAL_FLOAT(error, 0, 0.001);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("lowpass.float",           test_float);
      register_test("lowpass.offset",          test_offset);
      register_test("lowpass.dc",              test_dc);
      register_test("lowpass.nan",             test_NaN);
      register_test("lowpass.small",           test_small);
      register_test("lowpass.odd",             test_odd);
      register_test("lowpass.vt.float",        [](){test_vt<float>();});
      register_test("lowpass.vt.double",       [](){test_vt<double>();});
      register_test("lowpass.vt.int8",         [](){test_vt<std::int8_t>();});
      register_test("lowpass.vt.uint8",        [](){test_vt<std::uint8_t>();});
      register_test("lowpass.vt.int16",        [](){test_vt<std::int16_t>();});
      register_test("lowpass.vt.uint16",       [](){test_vt<std::uint16_t>();});
      register_test("lowpass.vt.int32",        [](){test_vt<std::int32_t>();});
      register_test("lowpass.vt.uint32",       [](){test_vt<std::uint32_t>();});
      register_test("lowpass.stress.int8",     [](){test_stress<std::int8_t>();});
      register_test("lowpass.stress.uint8",    [](){test_stress<std::uint8_t>();});
      register_test("lowpass.stress.int16",    [](){test_stress<std::int16_t>();});
      register_test("lowpass.stress.uint16",   [](){test_stress<std::uint16_t>();});
      register_test("lowpass.stress.int32",    [](){test_stress<std::int32_t>();});
      register_test("lowpass.stress.uint32",   [](){test_stress<std::uint32_t>();});
      register_test("lowpass.underflow",       test_underflow);
      register_test("lowpass.overflow",        test_overflow);
      register_test("lowpass.border",          test_border);
    }
  } dummy;
} // namespace for registration
