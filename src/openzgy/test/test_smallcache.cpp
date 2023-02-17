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
#include "../api.h"
#include "../iocontext.h"
#include "../impl/environment.h"

#include <iostream>
#include <iomanip>
#include <memory>

using namespace OpenZGY;

namespace Test {
  class TestSmallCache {
  public:
    static void test_regular();
    static void test_fewbricks();
    static void test_manybricks();
  private:
    static std::int64_t open_with_trace(const std::string& filename);
  };
}

/**
 * Suitable for use as a _debug_trace callback.
 *
 * Can print the calls immediately for manual debuggig and also
 * collect them for later checking in unit tests. Currently in the C++
 * version only the count is checked.
 */
class TraceCallsToSD
{
  typedef std::function<void(const std::string&, std::int64_t, std::int64_t, std::int64_t, const std::vector<std::int64_t>&)> debugtrace_t;
  typedef struct {std::string what; std::int64_t need; std::int64_t want; std::int64_t parts; std::vector<std::int64_t> all;} entry_t;

  bool _verbose;
  std::vector<entry_t> _calls;

  TraceCallsToSD(const TraceCallsToSD&) = delete;
  TraceCallsToSD& operator=(const TraceCallsToSD&) = delete;

public:
  explicit TraceCallsToSD(bool verbose=false)
    : _verbose(verbose)
  {
  }
  void operator()(const std::string& what,
      std::int64_t need,
      std::int64_t want,
      std::int64_t parts,
      const std::vector<std::int64_t>& all)
  {
    _calls.push_back(entry_t{what, need, want, parts, all});
    if (_verbose)
      std::cerr << std::setw(9) << what
                << " need " << std::setw(10) << _pretty(need)
                << " want " << std::setw(10) << _pretty(want)
                << " parts " << parts
                << std::endl;
  }
  static std::string _pretty(std::int64_t n)
  {
    if ((n < 1024) || (n % (1024) != 0))
      return std::to_string(n) + " bytes";
    else if ((n < 1024*1024) || (n % (1024*1024) != 0))
      return std::to_string(n/1024) + " KB";
    else
      return std::to_string(n/(1024*102)) + " MB";
  }
  void reset()
  {
    _calls.clear();
  }
  std::int64_t count() const
  {
    return _calls.size();
  }
};

/*
 * Data preparation.
 * The --BS option to zgycopyc is not yet implemented, so you will need
 * to patch the tool calling .bricksize() on the ZgyWriterArgs.
 *
 * zgycopyc --size 100,100,100,1 -BS 8,16,64 -o sd://sntc/testdata/fewbricks.zgy
 * zgycopyc --size 200,200,200,1 -BS 4,4,4 -o sd://sntc/testdata/manybricks.zgy
 *
 * fewbricks.zgy has 221 bricks including lod, needs ~2Kb headers.
 * bricksize is 8 KB so all headers fit in one brick.
 * Size of segment 0 is thus 8 KB. Verify by enabling the log
 * statements near the end of the SeismicStoreFile constructor,
 * then use zgycopyc to read the file and discard the results.
 * zgy: Total 1818624 in 2 segments of size 8192, 1810432, 1810432
 *
 * manybricks.zgy has 143166 bricks including lod0, needs ~1200 KB.
 * bricksize is 64 bytes so all headers do not fit in one brick,
 * nor do they fit in the default cachesize of 256 KB.
 * zgy: Total 10337856 in 2 segments of size 1175232, 9162624, 9162624
 *
 * Synt2.zgy is an old test file, see elsewhere for instructions to create.
 */

/**
 * Without the meta cache, opening a file takes 6 reads. With the
 * cache there will be just one. This applies to both local and cloud
 * access, but currently only cloud access has the debug_trace hook
 * allowing the unit test to verify it.
 */
std::int64_t
Test::TestSmallCache::open_with_trace(const std::string& name)
{
#ifdef HAVE_SD
  std::string filename = InternalZGY::Environment::getStringEnv("OPENZGY_SDTESTDATA", "sd://sntc/testdata");
  if (filename.back() != '/')
    filename += "/";
  filename += name;
  TraceCallsToSD tracer(verbose());
  if (verbose())
    std::cerr << std::endl;
  auto context = SeismicStoreIOContext()
    .debug_trace(std::ref(tracer));
  std::shared_ptr<OpenZGY::IZgyReader> reader = OpenZGY::IZgyReader::open(filename, &context);
  return tracer.count();
#else
  return 0;
#endif
}

void
Test::TestSmallCache::test_regular()
{
  TEST_CHECK(open_with_trace("Synt2.zgy") == 1);
}

void
Test::TestSmallCache::test_fewbricks()
{
  TEST_CHECK(open_with_trace("fewbricks.zgy") == 1);
}

void
Test::TestSmallCache::test_manybricks()
{
  TEST_CHECK(open_with_trace("manybricks.zgy") == 2);
}

namespace {
  class Register
  {
  public:
    Register()
    {
      using Test::TestSmallCache;
#ifdef HAVE_SD
      register_sd_test("smallcache.regular",      TestSmallCache::test_regular);
      register_sd_test("smallcache.fewbricks",    TestSmallCache::test_fewbricks);
      register_sd_test("smallcache.manybricks",   TestSmallCache::test_manybricks);
#endif
    }
  } dummy;
} // namespace for registration
