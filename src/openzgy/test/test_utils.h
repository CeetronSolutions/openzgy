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

#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <memory>
#include <functional>
#include <sstream>

namespace OpenZGY {
  class IOContext;
  class SeismicStoreIOContext;
}

#define TEST_EQUAL(actual,expect) TEST_CHECK_((actual)==(expect), "%s == %s: got %s, expect %s", #actual, #expect, Test_Utils::cutest_format((actual)).c_str(), Test_Utils::cutest_format((expect)).c_str())
#define TEST_EQUAL_FLOAT(actual,expect,epsilon) TEST_CHECK_(std::abs((actual)-(expect)) <= (epsilon), "%s == %s: got %s, expect %s", #actual, #expect, Test_Utils::cutest_format((actual)).c_str(), Test_Utils::cutest_format((expect)).c_str())
// Don't sue me if the next one doesn't work for even mildly complicated
// expressions. Instead, use must_throw directly. The drawback is that
// in that case you won't get the line number of the actual test.
#define MUST_THROW(expect, ...) TEST_CHECK(must_throw(expect, [&](){__VA_ARGS__}));

namespace Test_Utils {
#if 0
}
#endif

// Free functions
std::vector<float> random_vector(std::size_t size);
std::vector<double> random_double_vector(std::size_t size);
int random_delay(int ms);
#ifdef HAVE_SD
const OpenZGY::SeismicStoreIOContext* default_sd_context();
#endif
const OpenZGY::IOContext* default_context();
void compare_files(const std::string& a_name, const std::string& b_name, double epsilon, double lodepsilon);
bool must_throw(const char* expect, const std::function<void()>& fn);

template<typename T>
std::string
cutest_format(T value)
{
  std::stringstream ss;
  ss << std::boolalpha << value;
  return ss.str();
}

/**
 * Arrange for this explicitly named temporary file to be deleted when
 * the instance goes out of scope. FUTURE: both for seismic store and locally.
 *
 * Unlike the Python counterpart, no error is raised if the file does not
 * exist. The test code needs to explicitly delete the file if it wants
 * to control how errors are treated. In C++ it is not possible, or at least
 * not portable, to check whether the stack is being unwound due to a
 * pending exception. If it is then we definitely don't want to replace that
 * exception with a confusing "file not found". Also it is an exceptionally
 * bad idea to throw an error in the destructor.
 *
 * There is a fairly harmless race condition between the creation of the
 * instance and the actual creation of the file. If an exception is thrown
 * in that interval then a message is printed about not being able to delete
 * the file. You may alternatively disarm() the instance immediately on
 * creation and then arm() it again at the exact point you know the file
 * exists. Good luck. If you get an exception when creating a file on the
 * cloud, you might not know whether the file got created or not.
 *
 * The recommended naming convention is to use both the current time and
 * a per-process random number as a prefix. This makes it simpler to
 * clean up files that in spite of the aito delete got left behind.
 *
 * The class can also be used for files intended to be persistent, by
 * ensuring the file gets deleted if any error happened while creating
 * it. Call disarm() only when you are sure you want to keep it.
 * At that point your file is no longer considered temporary.
 */
class TempFileAutoDelete
{
public:
  TempFileAutoDelete(const std::string& name, const OpenZGY::IOContext* ctx = nullptr);
  ~TempFileAutoDelete();
  TempFileAutoDelete(const TempFileAutoDelete&) = delete;
  TempFileAutoDelete& operator=(const TempFileAutoDelete&) = delete;

  /**
   * The name of the file that will be deleted when the instance goes
   * out of scope.
   */
  const std::string& name() const { return name_; }

  /**
   * Tell the instance it is supposed to delete the file on exit.
   * Instances are armed by default.
   */
  void arm() { armed_ = true; }

  /**
   * Tell the instance to not try to delete the file. Either because we
   * deleted it ourselves, or because we for some reason decided to keep it,
   * or we want to guard agaist exceptions during create by temporarily
   * disarming and then arming the file.
   */
  void disarm() { armed_ = false; }

  /**
   * Combine a base folder and a file name or path into a single path.
   * If the second argument is empty this is an error and the result
   * will be empty as well.
   */
  static std::string join(const std::string& base, const std::string& file);

protected:
  static std::string randomname();
  static std::uint32_t myrand();
  static std::uint32_t myrand(std::uint32_t max);
  std::shared_ptr<OpenZGY::IOContext> savedcontext() const {return context_;}

private:
  bool armed_;
  std::string name_;
  std::shared_ptr<OpenZGY::IOContext> context_;
};

/**
 * As TempFileAutoDelete, but explicitly requesting a local file with
 * a random prefix in from of the name.
 */
class LocalFileAutoDelete : public TempFileAutoDelete
{
public:
  LocalFileAutoDelete(const std::string& suffix)
    : TempFileAutoDelete(getPrefix() + randomname() + suffix)
  {
  }

  /**
   * Check whether a file exists. Currently only supports on-prem files.
   * Declared public because it might be useful elsewhere.
   */
  static bool exists(const std::string& filename);

  /**
   * Delete file.
   * Declared public because it might be useful elsewhere.
   */
  static void remove(const std::string& filename);

private:
  static std::string getPrefix();
  static std::string makePrefix();
};

/**
 * As TempFileAutoDelete, but explicitly requesting a file on the
 * seismic store with a random prefix in from of the name.
 */
class CloudFileAutoDelete : public TempFileAutoDelete
{
public:
  CloudFileAutoDelete(const std::string& suffix, const OpenZGY::SeismicStoreIOContext* ctx);
  static void remove(const std::string& name, const OpenZGY::SeismicStoreIOContext* ctx);
private:
  static std::string getPrefix();
  static std::string makePrefix();
};

} // namespace
