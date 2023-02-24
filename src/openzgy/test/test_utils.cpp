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

// Based on Salmon/Shared/TestUtils/TempFileAutoDelete.cpp
// And wrapper/test_utils.py

#include "test_utils.h"
#include "test_all.h"
#include "../impl/environment.h"
#include "../impl/logger.h"

#include <sstream>
#include <chrono>
#include <algorithm>
#include <random>
#include <functional>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <thread>
#include <iostream>
#include <numeric>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#include <io.h>
#else
#include <unistd.h>
#endif

#include "../api.h"
#include "../iocontext.h"
#include "../exception.h"
#include "../impl/file.h"

using InternalZGY::Environment;

namespace Test_Utils {
#if 0
}
#endif

TempFileAutoDelete::TempFileAutoDelete(const std::string& name, const OpenZGY::IOContext* ctx)
  : armed_(true)
  , name_(name)
  , context_(ctx ? ctx->clone() : std::shared_ptr<OpenZGY::IOContext>())
{
}

TempFileAutoDelete::~TempFileAutoDelete()
{
  if (armed_ && !name_.empty()) {
    try {
      // Code smell, the reason is that I don't want virtuals here.
      if (name_.substr(0,5) == "sd://") {
        CloudFileAutoDelete::remove(name_, dynamic_cast<OpenZGY::SeismicStoreIOContext*>(savedcontext().get()));
      }
      else {
        LocalFileAutoDelete::remove(name_);
      }
      armed_ = false;
    }
    catch (const std::exception& ex) {
      std::cerr << "WARNING: Could not remove temporary file: " << ex.what() << "\n";
    }
  }
}

/**
 * Return an integer between 0 and 32767 inclusive.
 * Don't use this for anything important.
 * The number generator has very little entropy.
 */
std::uint32_t
TempFileAutoDelete::myrand()
{
  static std::uint32_t seed = static_cast<std::uint32_t>(time(0));
  seed = seed * 1103515245 + 12345;
  return((unsigned)(seed / 65536) % 32768);
}

/**
 * Return an integer between 0 and max inclusive.
 * Don't use this for anything important.
 * The number generator has very little entropy.
 */
std::uint32_t
TempFileAutoDelete::myrand(std::uint32_t max)
{
  return myrand() / (32767 / max);
}

bool
LocalFileAutoDelete::exists(const std::string &filename)
{
#ifdef WIN32
  return ::GetFileAttributesA(filename.c_str()) != INVALID_FILE_ATTRIBUTES;
#else
  return ::access(filename.c_str(), F_OK) >= 0;
#endif
}

void
LocalFileAutoDelete::remove(const std::string &filename)
{
  //std::cerr << "@ AUTO_DELETE \"" << filename << "\"\n";
  // TODO-Low, print warning on failure. Maybe support throw_on_error.
#ifdef WIN32
  ::DeleteFileA(filename.c_str());
#else
  ::remove(filename.c_str());
#endif
}

std::string
TempFileAutoDelete::join(const std::string& base, const std::string& file)
{
  if (file.empty())
    return file;
  else if (base.empty() || file.front() == '/' || file.front() == '\\')
    return file;
  else if (base.back() == '/' || base.back() == '\\')
    return base + file;
  else {
#ifdef WIN32
    return base + '\\' + file;
#else
    return base + '/' + file;
#endif
  }
}

std::string
TempFileAutoDelete::randomname()
{
  // The random number generator is low quality, but a collision
  // in temp file name clash isn't actually an earth shattering
  // bug. At worst some unit test might fail.
  std::uint32_t rnd = myrand();
  std::stringstream ss;
  ss << std::hex << std::setw(8) << std::setfill('0') << rnd << "-";
  return ss.str();
}

std::string
LocalFileAutoDelete::makePrefix()
{
  std::stringstream ss;
  ss << "tmp-"
     << std::hex << std::setw(8) << std::setfill('0')
     << std::uint32_t(time(nullptr))
     << "-";
  return join(Environment::getStringEnv("TESTRUNDIR", "."), ss.str());
}

std::string
LocalFileAutoDelete::getPrefix()
{
  // Every call will return the same value.
  static const std::string prefix = makePrefix();
  return prefix;
}

CloudFileAutoDelete::CloudFileAutoDelete(const std::string& suffix, const OpenZGY::SeismicStoreIOContext* ctx)
  : TempFileAutoDelete(getPrefix() + randomname() + suffix, ctx)
{
}

std::string
CloudFileAutoDelete::makePrefix()
{
  std::stringstream ss;
  ss << "/tmp-"
     << std::hex << std::setw(8) << std::setfill('0')
     << std::uint32_t(time(nullptr))
     << "-";
  return Environment::getStringEnv("OPENZGY_SDTESTSINK", "sd://sntc/testsink/d") + ss.str();
}

std::string
CloudFileAutoDelete::getPrefix()
{
  // Every call will return the same value.
  static const std::string prefix = makePrefix();
  return prefix;
}

void
CloudFileAutoDelete::remove(const std::string& name, const OpenZGY::SeismicStoreIOContext* ctx)
{
#ifdef HAVE_SD
  // std::cerr << "DELETE FROM CLOUD \"" << name << "\"\n";
  if (!name.empty()) {
    if (ctx) {
      // Using the low level delete functionality;
      // I could also have used OpenZGY::IZgyUtils.
      // The low level code needs a logger attached here,
      // unless I am ok with the default.
      OpenZGY::SeismicStoreIOContext localctx(*ctx);
      localctx.logger
        (InternalZGY::LoggerBase::standardCallback
         (InternalZGY::LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"),
          "openzgy-testutils: ", ""));
      std::shared_ptr<InternalZGY::IFileADT> fd = InternalZGY::FileFactory::instance().create
        (name, InternalZGY::OpenMode::Closed, &localctx);
      fd->deleteFile(name, /*missing_ok=*/true);
    }
    else {
      throw OpenZGY::Errors::ZgyInternalError("CloudFileAutoDelete::remove requires a context");
    }
  }
#endif
}

/**
 * \brief Create a vecor of random values.
 *
 * Is it just me, or is this horribly convoluted compared to the
 * equivalent function in Python?
 */
std::vector<float>
random_vector(std::size_t size)
{
  static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(static_cast<unsigned>(seed));
  static std::normal_distribution<float> distribution(0.0, 1.0);
  static std::function<float()> rnd = [&](){return distribution(generator);};
  std::vector<float> result(size);
  std::generate(result.begin(), result.end(), rnd);
  return result;
}

std::vector<double>
random_double_vector(std::size_t size)
{
  static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(static_cast<unsigned>(seed));
  static std::normal_distribution<double> distribution(0.0, 1.0);
  static std::function<double()> rnd = [&](){return distribution(generator);};
  std::vector<double> result(size);
  std::generate(result.begin(), result.end(), rnd);
  return result;
}

/**
 * Delay execution of the current thread by approximately the provided
 * number of milliseconds but varying between half and double the
 * specified time.
 */
int
random_delay(int ms)
{
  static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(static_cast<unsigned>(seed));
  std::uniform_int_distribution<int> distribution(ms/2, ms*2);
  int sleeptime = distribution(generator);
  std::this_thread::sleep_for(std::chrono::milliseconds(sleeptime));
  return sleeptime;
}

#ifdef HAVE_SD
/**
 * Convenience to hard code credentials for testing.
 */
const OpenZGY::SeismicStoreIOContext*
default_sd_context()
{
  // Picking up sdurl/sdapikey from the environment is redundant since
  // the library already does this as a fallback.
  using InternalZGY::Environment;
  static OpenZGY::SeismicStoreIOContext instance =
    OpenZGY::SeismicStoreIOContext()
    // Enable to prove that logging is configurable, for seismic storeat least.
#if 0
    .logger(InternalZGY::LoggerBase::standardCallback
            (InternalZGY::LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"),
             "openzgy-unittest: ", ""))
#endif
    .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
    .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
    .sdtoken(Environment::getStringEnv("OPENZGY_TOKEN"), "");
  return &instance;
}
#endif

const OpenZGY::IOContext*
default_context()
{
#ifdef HAVE_SD
  return default_sd_context();
#else
  return nullptr;
#endif
}

namespace {
  template<typename T, std::size_t  N>
  std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a)
  {
    os << "[";
    for (std::size_t ii=0; ii<N; ++ii)
      os << a[ii] << (ii == N-1 ? "" : ", ");
    os << "]";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const OpenZGY::SampleStatistics& in)
  {
    os << "cnt: " << in.cnt
       << " sum: " << in.sum
       << " ssq: " << in.ssq
       << " min: " << in.min
       << " max: " << in.max;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const OpenZGY::SampleHistogram& in)
  {
    os << "cnt: " << in.samplecount
       << " min: " << in.minvalue
       << " max: " << in.maxvalue
       << " bincount: " << in.bins.size();
    return os;
  }
}

static bool
similar(double a, double b, double eps)
{
  return std::abs(a - b)  <= eps * 0.5 * (std::abs(a) + std::abs(b));
}

static void
compare_stats(const OpenZGY::SampleStatistics& a, const OpenZGY::SampleStatistics& b, double eps)
{
  TEST_CHECK(a.cnt == b.cnt);
  TEST_CHECK(similar(a.sum, b.sum, eps));
  TEST_CHECK(similar(a.ssq, b.ssq, eps));
  TEST_CHECK(similar(a.min, b.min, eps));
  TEST_CHECK(similar(a.max, b.max, eps));
}

static void
compare_histo(const OpenZGY::SampleHistogram& a, const OpenZGY::SampleHistogram& b, double eps)
{
  // Consistency checks done on both histograms in turn.
  // The histogram on file stores sample_count. If there are no NaN, Inf,
  // or outside-range this should match the sum of all bins.
  const std::int64_t acount =
    std::accumulate(a.bins.begin(), a.bins.end(), std::int64_t(0));
  const std::int64_t bcount =
    std::accumulate(b.bins.begin(), b.bins.end(), std::int64_t(0));
  TEST_CHECK(a.samplecount == acount);
  TEST_CHECK(b.samplecount == bcount);
  // Now compare the two histograms.
  TEST_CHECK(a.bins.size() == b.bins.size());
  TEST_CHECK(a.samplecount == b.samplecount);
  TEST_CHECK(similar(a.minvalue, b.minvalue, eps));
  TEST_CHECK(similar(a.maxvalue, b.maxvalue, eps));
  for (std::size_t ii = 0; ii < std::min(a.bins.size(), b.bins.size()); ++ii) {
    TEST_CHECK(a.bins[ii] == b.bins[ii]);
    if (verbose())
      std::cout << "[" << ii << "] "
                << std::setw(8) << a.bins[ii] << " "
                << std::setw(8) << b.bins[ii] << " "
                << std::setw(8) << (a.bins[ii] - b.bins[ii]) << "\n";
  }
}

/**
 * Compare two files, possibly with differing value types and
 * bricksize. Bulk is compared after conversion to float. Meant to be
 * used for small files only. The entire file will be read into
 * memory. Any sample that is NaN will cause the compare to fail.
 */
void
compare_files(const std::string& a_name, const std::string& b_name, double epsilon, double lodepsilon)
{
  std::shared_ptr<OpenZGY::IZgyReader> a_reader =
    OpenZGY::IZgyReader::open(a_name, Test_Utils::default_context());
  std::shared_ptr<OpenZGY::IZgyReader> b_reader =
    OpenZGY::IZgyReader::open(b_name, Test_Utils::default_context());
  const OpenZGY::IZgyReader& a = *a_reader;
  const OpenZGY::IZgyReader& b = *b_reader;
  std::shared_ptr<const OpenZGY::FileStatistics>a_filestats = a.filestats();
  std::shared_ptr<const OpenZGY::FileStatistics>b_filestats = b.filestats();
  const OpenZGY::IZgyMeta::corners_t a_corners = a.corners();
  const OpenZGY::IZgyMeta::corners_t b_corners = b.corners();

  if (verbose()) {
    std::cout << "Statistics\n"
              << "a: " << a.statistics() << "\nb: " << b.statistics() << "\n";
    std::cout << "Histogram\n"
              << "a: " << a.histogram() << "\nb: " << b.histogram() << "\n";
  }

  TEST_CHECK(a.size() == b.size());
  TEST_CHECK(a.nlods() == b.nlods());
  TEST_CHECK(similar(a.datarange()[0], b.datarange()[0], epsilon));
  TEST_CHECK(similar(a.datarange()[1], b.datarange()[1], epsilon));
  TEST_CHECK(similar(a.annotstart()[0], b.annotstart()[0], 0.01));
  TEST_CHECK(similar(a.annotstart()[1], b.annotstart()[1], 0.01));
  TEST_CHECK(similar(a.annotinc()[0],   b.annotinc()[0], 0.01));
  TEST_CHECK(similar(a.annotinc()[1],b.annotinc()[1], 0.01));
  TEST_CHECK(similar(a.zstart(), b.zstart(), 0.01));
  TEST_CHECK(similar(a.zinc(), b.zinc(),  0.01));
  TEST_CHECK(a.hunitdim()    == b.hunitdim());
  TEST_CHECK(a.hunitname()   == b.hunitname());
  TEST_CHECK(a.hunitfactor() == b.hunitfactor());
  TEST_CHECK(a.zunitdim()    == b.zunitdim());
  TEST_CHECK(a.zunitname()   == b.zunitname());
  TEST_CHECK(a.zunitfactor() == b.zunitfactor());
  for (int corner = 0; corner < 4; ++corner) {
    for (int dim = 0; dim < 2; ++dim) {
      TEST_CHECK(a_corners[corner][dim] - b_corners[corner][dim] <= 3.0);
    }
  }
  if (a_filestats->fileVersion() >= 3 && a_filestats->fileVersion() >= 3) {
    compare_stats(a.statistics(), b.statistics(), epsilon);
    compare_histo(a.histogram(),  b.histogram(), epsilon);
    TEST_CHECK(a_name == b_name || a.verid() != b.verid());
    // Files with an integral type normally have histogram range set to
    // codingrange. I.e. we assume both the extreme minimum and the
    // extreme maximum are in use instead of checking. This might change
    // in the future for int16 files, so a failure here might not mean
    // there is anything wrong.
    if (a.datatype() != OpenZGY::SampleDataType::float32) {
      TEST_CHECK(similar(a.datarange()[0], a.histogram().minvalue, 1.0e-5));
      TEST_CHECK(similar(a.datarange()[1], a.histogram().maxvalue, 1.0e-5));
    }
    if (b.datatype() != OpenZGY::SampleDataType::float32) {
      TEST_CHECK(similar(b.datarange()[0], b.histogram().minvalue, 1.0e-5));
      TEST_CHECK(similar(b.datarange()[1], b.histogram().maxvalue, 1.0e-5));
    }
  }
  if (a.size() == b.size()) {
    const std::int64_t count = a.size()[0] * a.size()[1] * a.size()[2];
    const std::int32_t nlods = std::min(a.nlods(), b.nlods());
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const std::array<std::int64_t,3> orig{0,0,0};
    std::unique_ptr<float[]> a_data(new float[count]);
    std::unique_ptr<float[]> b_data(new float[count]);
    std::fill(a_data.get(), a_data.get() + count, nan);
    std::fill(b_data.get(), b_data.get() + count, nan);
    float *aptr = a_data.get();
    float *bptr = b_data.get();
    std::array<std::int64_t,3> size{
      std::min(a.size()[0], b.size()[0]),
      std::min(a.size()[1], b.size()[1]),
      std::min(a.size()[2], b.size()[2]),
    };
    for (int lod=0; lod<nlods; ++lod) {
      a.read(orig, size, aptr, lod);
      b.read(orig, size, bptr, lod);
      float worst = -1;
      bool reported_error = false;
      for (int ii=0; ii<count; ++ii) {
        worst = std::max(worst, std::abs(aptr[ii] - bptr[ii]));
        if (!reported_error) {
          TEST_CHECK(similar(aptr[ii], bptr[ii], (lod?lodepsilon:epsilon)));
          if (!similar(aptr[ii], bptr[ii], (lod?lodepsilon:epsilon))) {
            std::cout << "MISMATCH at offset " << ii
                      << " lod: " << lod
                      << ", a=" << aptr[ii]
                      << " != b=" << bptr[ii]
                      << "\n";
            // avoid flooding the log with failures.
            reported_error = true;
          }
        }
      }
      if (verbose()) {
            std::cout << "Most noise in lod " << lod
                      << " is " << worst
                      << "\n";
      }
      for (int ii=0; ii<3; ++ii)
        size[ii] = (size[ii]+1)/2;
    }
  }
}

bool
must_throw(const char *expect, const std::function<void()>& fn)
{
  try {
    fn();
    return TEST_CHECK_(false, "did not get expected exception \"%s\"\n", expect);
  }
  catch (const OpenZGY::Errors::ZgyError& ex) {
    if (verbose()) {
      if (strstr(ex.what(), expect) != nullptr)
        std::cout << "Got expected exception: " << ex.what() << std::endl;
    }
    return TEST_CHECK_(strstr(ex.what(), expect) != nullptr, "got wrong exception \"%s\" instead of \"%s\"\n", ex.what(), expect);
  }
}

}
