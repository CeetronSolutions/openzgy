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
#include "../impl/meta.h"
#include "../impl/file.h"
#include "../impl/lookuptable.h"
#include "../impl/environment.h"

#include <iostream>
#include <memory>

using namespace InternalZGY;

namespace Test {
  class TestMeta {
  public:
    static void test_open_v1();
    static void test_open_v3();
    static void test_blup_v1();
    static void test_blup_v3();
  private:
    static void test_blup(const std::string&);
    static std::string get_testdata(const std::string& name);
  };
}

std::string
Test::TestMeta::get_testdata(const std::string& name)
{
#ifdef _WIN32
  std::string result = Environment::getStringEnv("OPENZGY_TESTDATA", "..\\..\\build\\testdata");
  if (result.back() != '\\')
    result += "\\";
  result += name;
  return result;
#else
  return
    Environment::getStringEnv("OPENZGY_TESTDATA", "../../build/testdata")
    + "/" + name;
#endif
}

void
Test::TestMeta::test_open_v1()
{
  const std::string filename = get_testdata("Empty-v1.zgy");
  std::shared_ptr<IFileADT> file =
    FileFactory::instance().create(filename, OpenMode::ReadOnly, nullptr);
  ZgyInternalMeta m(file);
  if (verbose()) {
    m.dump(std::cout, "  test:");
  }
}

void
Test::TestMeta::test_open_v3()
{
  const std::string filename = get_testdata("Empty-v3.zgy");
  std::shared_ptr<IFileADT> file =
    FileFactory::instance().create(filename, OpenMode::ReadOnly, nullptr);
  ZgyInternalMeta m(file);
  if (verbose()) {
    m.dump(std::cout, "  test:");
  }
}

#if 0
// Format like zgydump --offset does,
// so I can manually check the result with a diff.
static std::string
zgydump_format_pos(std::int64_t ii, std::int64_t jj, std::int64_t kk, std::int64_t lod)
{
  std::stringstream ss;
  ss << "[" << lod << "]"
     << "[" << ii << "]"
     << "[" << jj << "]"
     << "[" << kk << "]";
  return ss.str();
}

static std::string
zgydump_format(std::int64_t ii, std::int64_t jj, std::int64_t kk, std::int64_t lod, const LookupTable::LutInfo& info)
{
  std::stringstream ss;
  ss << std::setw(20) << std::left << zgydump_format_pos(ii, jj, kk, lod)
     << " = " << std::hex << std::setw(16) << std::right << info.offset_in_file
     << " Size " << std::setw(8) << info.size_in_file;
  return ss.str();
}
#endif

void
Test::TestMeta::test_blup(const std::string& filename)
{
#if 0 // Temporarily disabled, _lookup not accessible.
  std::shared_ptr<IFileADT> file =
    FileFactory::instance().create(filename, OpenMode::ReadOnly, nullptr);
  ZgyInternalMeta m(file);
  const auto& blup = dynamic_cast<const LookupTableV0Access&>(*m._blup);

  std::array<std::int64_t,3> size = m._ih->size();
  const std::array<std::int64_t,3> bs = m._ih->bricksize();
  for (std::int64_t lod= 0; lod< m._ih->nlods(); ++lod) {
    for (std::int64_t ii = 0; ii < size[0]; ii += bs[0]) {
      for (std::int64_t jj = 0; jj < size[1]; jj += bs[1]) {
        for (std::int64_t kk = 0; kk < size[2]; kk += bs[2]) {
          LookupTable::LutInfo info =
            LookupTable::getBrickFilePosition(ii/bs[0], jj/bs[1], kk/bs[2], lod,
                                 m._ih->lodsizes(), m._ih->brickoffsets(),
                                 blup._lookup,
                                 blup._lookend,
                                 64*64*64);
          if (verbose())
            std::cout << zgydump_format(ii/bs[0], jj/bs[1], kk/bs[2], lod, info) << "\n";
        }
      }
    }
    size[0] = (size[0]+1)/2;
    size[1] = (size[1]+1)/2;
    size[2] = (size[2]+1)/2;
  }
#endif
}

void
Test::TestMeta::test_blup_v1()
{
  test_blup(get_testdata("Empty-v1.zgy"));
}

void
Test::TestMeta::test_blup_v3()
{
  test_blup(get_testdata("Empty-v3.zgy"));
}

namespace {
  class Register
  {
  public:
    Register()
    {
      using Test::TestMeta;
      register_test("meta.open_v1",            TestMeta::test_open_v1);
      register_test("meta.open_v3",            TestMeta::test_open_v3);
      register_test("meta.blup_v1",            TestMeta::test_blup_v1);
      register_test("meta.blup_v3",            TestMeta::test_blup_v3);
    }
  } dummy;
} // namespace for registration
