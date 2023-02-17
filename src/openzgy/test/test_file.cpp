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
#include "../exception.h"
#include "../iocontext.h"
#include "../impl/file.h"
#include "../impl/file_smallcache.h"
#include "../impl/environment.h"

#include <iostream>
#include <memory>

// For fork/exec
#ifdef _WIN32
#else
#include <sys/unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif
#include <stdio.h>

using namespace OpenZGY;
using namespace OpenZGY::Errors;
using namespace InternalZGY;
using Test_Utils::LocalFileAutoDelete;
using Test_Utils::CloudFileAutoDelete;
using Test_Utils::must_throw;

namespace {
#if 0
}
#endif

static void
deleteFile(const std::shared_ptr<IFileADT>& file, const std::string& filename)
{
  file->deleteFile(filename, /*missing_ok=*/true);
}

static void
run_helloworld(const std::string& fn, const std::function<std::shared_ptr<IFileADT>(const std::string&, OpenMode)>& alloc)
{
  std::shared_ptr<IFileADT> file;
  file = alloc(fn, OpenMode::Truncate);
  file->xx_write("Hello, world, what else can I say?\n", 0, 35);
  file->xx_close();
  file.reset();

  file = alloc(fn, OpenMode::ReadOnly);
  char data[100]{0};
  file->xx_read(data, 7, 5);
  TEST_CHECK(std::string(data) == std::string("world"));
  auto sink = [&data](ReadRequest::data_t ptr, std::int64_t size) {
                memcpy(data, ptr.get(), size);
              };
  ReadList requests{ ReadRequest(14, 9, sink) };
  file->xx_readv(requests);
  TEST_CHECK(std::string(data) == std::string("what else"));

  // Also test a wrapper to cache the first 16(!) bytes
  file.reset(new FileWithSmallCache(file, 16));
  std::fill(&data[0], &data[100], 0);
  file->xx_read(data, 7, 5);
  TEST_CHECK(std::string(data) == std::string("world"));
  std::fill(&data[0], &data[100], 0);
  file->xx_read(data, 1, 4); // will be cached
  TEST_CHECK(std::string(data) == std::string("ello"));
  std::fill(&data[0], &data[100], 0);
  file->xx_readv(requests); // not cached, crosses end of the cache
  TEST_CHECK(std::string(data) == std::string("what else"));
  file->xx_close();
  file.reset();

  file = alloc(fn, OpenMode::Closed);
  deleteFile(file, fn);
  file.reset();

  // Attempt to delete non-existing file is not an error.
  file = alloc(fn, OpenMode::Closed);
  deleteFile(file, fn);
  file.reset();

  // Attempt to open a non-existing file is an error.
  try {
    file = alloc(fn, OpenMode::ReadOnly);
    TEST_CHECK(false && "Expected an exception about non-existing file.");
  }
  catch (const ZgyIoError& ex) {
    TEST_CHECK(ex.what() && std::string(ex.what()).find("No such file") != std::string::npos);
  }
  catch (const ZgyInternalError& ex) {
    // TODO-Low: A specific cloud-related error.
    // And/or a specific "file not found" error.
    TEST_CHECK(ex.what() && std::string(ex.what()).find("does not exist") != std::string::npos);
  }
}

static void
test_localfilefactory()
{
  LocalFileAutoDelete lad("testfile.zgy");
  run_helloworld(lad.name(),
                 [](const std::string& name, OpenMode mode) {
                   return FileFactory::instance().create(
                       name, mode, nullptr);
                 });
}

#ifdef HAVE_SD
static void
test_sdfilefactory()
{
  CloudFileAutoDelete cad("ozcpp-testfile.txt", Test_Utils::default_sd_context());
  run_helloworld(cad.name(),
                 [](const std::string& name, OpenMode mode) {
                   return FileFactory::instance().create(
                       name, mode, Test_Utils::default_sd_context());
                 });
  cad.disarm();
}

/**
 * This is primarily to verify that CloudFileAutoDelete works.
 */
static void
test_sdfiledelete()
{
  std::string filename;
  {
    CloudFileAutoDelete cad("ozcpp-testdelete.txt", Test_Utils::default_sd_context());
    filename = cad.name();
    std::shared_ptr<IFileADT> file = FileFactory::instance().create
      (filename, OpenMode::Truncate, Test_Utils::default_sd_context());
    file->xx_write("Hello, world, I will be auto-deleted\n", 0, 37);
    file->xx_close();
    file.reset();
  }

  // cad has gone out of scope, so the file ought to have been deleted.
  try {
    std::shared_ptr<IFileADT> file = FileFactory::instance().create
      (filename, OpenMode::ReadOnly, Test_Utils::default_sd_context());
    TEST_CHECK(false && "Expected an exception about non-existing file.");
  }
  catch (const ZgyIoError& ex) {
    if (verbose())
      std::cout << "Got expected I/O exception " << ex.what() << std::endl;
    TEST_CHECK(ex.what() && std::string(ex.what()).find("No such file") != std::string::npos);
  }
  catch (const ZgyInternalError& ex) {
    if (verbose())
      std::cout << "Got expected internal exception " << ex.what() << std::endl;
    // TODO-Low: A specific cloud-related error.
    // And/or a specific "file not found" error.
    TEST_CHECK(ex.what() && std::string(ex.what()).find("does not exist") != std::string::npos);
  }
}

/**
 * Python one-liner to show the expiry time of the OPENZGY_TOKEN
 * if one is present. No error handling whatsoever. may be useful to
 * debug problems in the build server. Note that this needs to be caller
 * before any test that will access the seismic store, so it should be
 * registered as "aaa.sdtoken" since tests are run alphabetically.
 *
 * It would of course be better to do this in C++ but I don't want to have
 * base64 and jsoncpp dependencies just because of this.
 */
static void
test_sdtoken()
{
  std::string token = Environment::getStringEnv("OPENZGY_TOKEN");
  auto pos = token.find(':');
  if (token.empty()) {
  }
  else if (pos != std::string::npos) {
    std::cout << "$OPENZGY_TOKEN is of type "
              << "\"" << token.substr(0, pos) << "\"" << std::endl;
  }
  else {
    std::cerr << std::flush;
    std::cout << std::flush;
    fflush(stderr);
    fflush(stdout);
#ifdef _WIN32
    // I could probably have used system() on Linux as well, but who cares?
    system("python -c \"import base64, json, os, datetime; print('\\nToken ending with', os.getenv('OPENZGY_TOKEN')[-5:], datetime.datetime.fromtimestamp(json.loads(base64.urlsafe_b64decode(os.getenv('OPENZGY_TOKEN').split('.')[1] + '====').decode())['exp']).strftime('expires %d/%m/%y %H:%M'))\"");
#else
    auto pid = fork();
    switch (pid) {
    case 0:
      execlp("python3", "python3", "-c", "import base64,json,os,datetime; print('\\nToken ending with', os.getenv('OPENZGY_TOKEN')[-5:],datetime.datetime.fromtimestamp(json.loads(base64.urlsafe_b64decode(os.getenv('OPENZGY_TOKEN').split('.')[1]+'====').decode())['exp']).strftime('expires %d/%m/%y %H:%M'))", (char*)0);
      _exit(0);
      break;
    case -1:
      break;
    default:
      waitpid(pid, nullptr, 0);
      break;
    }
#endif
  }
}

static void
update_sdreopen(const std::string& filename, std::int64_t start, std::int64_t nbytes, bool tweak)
{
  OpenZGY::SeismicStoreIOContext context =
    OpenZGY::SeismicStoreIOContext(*Test_Utils::default_sd_context())
    .setRoAfterWrite(false)
    .segsizebytes(2048)
    .segsplit(nbytes==10000 ? 8 : 1);
  std::shared_ptr<IFileADT> fd =
    FileFactory::instance().create
    (filename,
     start==0 && nbytes == 256 && !tweak ? OpenMode::Truncate : OpenMode::ReadWrite,
     &context);
  std::vector<std::int16_t> data(nbytes/sizeof(std::int16_t));
  for (std::int64_t ii=0; ii<(std::int64_t)data.size(); ++ii)
    data[ii] = static_cast<std::int16_t>((start / sizeof(std::int16_t)) + ii);
  if (tweak && data.size() > 42)
    data[42] *= -1;
  fd->xx_write(data.data(), start, nbytes, UsageHint::Unknown);
  fd->xx_close();
  fd.reset();
}

static bool
check_sdreopen(const std::string& filename, std::int64_t expect_bytes, bool tweak)
{
  bool ok = true;
  std::shared_ptr<IFileADT> fd =
    FileFactory::instance().create
    (filename, OpenMode::ReadOnly, Test_Utils::default_context());
  // TODO-@@@: Currently xx_segments() won't check all blocks.
  // Can change this but beware use of the debug callback gets expensive.
  const std::vector<std::int64_t> segs = fd->xx_segments(true);
  const std::int64_t expect_num_segs = 1 + (expect_bytes - 256 + 2047) / 2048;
  const std::int64_t expect_lastseg = ((expect_bytes + 1792 - 1) % 2048) + 1;
  if (TEST_EQUAL((std::int64_t)segs.size(), expect_num_segs)) {
    for (std::int64_t ii=0; ii<expect_num_segs; ++ii)
      if (!TEST_EQUAL(segs[ii], (ii==0 ? 256 :
                                 ii==expect_num_segs-1 ? expect_lastseg :
                                 2048)))
      {
        ok = false;
      }
  }
  else {
    ok = false;
  }
  if (TEST_EQUAL(fd->xx_eof(), expect_bytes)) {
    std::vector<std::int16_t> data(expect_bytes/sizeof(std::int16_t));
    fd->xx_read(data.data(), 0, expect_bytes, UsageHint::Unknown);
    for (std::int64_t ii=0; ii < (std::int64_t)data.size(); ++ii) {
      std::int16_t e = static_cast<std::int16_t>((tweak && ii==42) ? -ii : ii);
      if (data[ii] != e) {
        if (!TEST_EQUAL(data[ii], e)) {
          ok = false;
          break;
        }
      }
    }
  }
  else {
    ok = false;
  }
  fd->xx_close();
  fd.reset();
  return ok;
}

static void
test_sdreopen()
{
  // Segment size is 2048 bytes except seg0 which is 256 bytes.
  //
  // The test will make 6 writes to the file,
  //    Bytes written:    256, 2048+42,   12, 4096, 2048-42-12,  998, 10000
  //    Total size:       256,    2346, 2358, 6454,       8448, 9446, 19446
  //    Segments:           1,       3,    3,    5,          5,    6,    11
  //
  // Rationale for choosing those sizes is explained before each test.
  //
  // In (2) and (6) there is no need to re-open the last segment but in (6)
  // it will probably happen anyway to simplify the code.

  CloudFileAutoDelete cad("sdreopen.dat", Test_Utils::default_sd_context());

  // (1) Initial create and establish size of segment 0.
  update_sdreopen(cad.name(), 0, 256, false);
  TEST_CHECK(check_sdreopen(cad.name(), 256, false));

  // (2) Write a bit more that one segment, last segment will be odd size.
  update_sdreopen(cad.name(), 256, 2048+42, false);
  TEST_CHECK(check_sdreopen(cad.name(), 2346, false));

  // (3) Write less than one segment, last segment odd size before and after.
  update_sdreopen(cad.name(), 2346, 12, false);
  TEST_CHECK(check_sdreopen(cad.name(), 2358, false));

  // (4) Write more than one segment, last segment odd size before and after.
  update_sdreopen(cad.name(), 2358, 4096, false);
  TEST_CHECK(check_sdreopen(cad.name(), 6454, false));

  // (5) Write less than one segment, exactly filling up last segment.
  update_sdreopen(cad.name(), 6454, 2048-42-12, false);
  TEST_CHECK(check_sdreopen(cad.name(), 8448, false));

  // (6) Write less than one segment, last segment full before, odd after.
  update_sdreopen(cad.name(), 8448, 998, false);
  TEST_CHECK(check_sdreopen(cad.name(), 9446, false));

  // (7) Write lots of data with multi threaded upload.
  update_sdreopen(cad.name(), 9446, 10000, false);
  TEST_CHECK(check_sdreopen(cad.name(), 19446, false));

  // (8) Overwrite segment 0, no change in size.
  update_sdreopen(cad.name(), 0, 256, true);
  TEST_CHECK(check_sdreopen(cad.name(), 19446, true));

  // (*) Overwrite segment 0 with a smaller size.
  must_throw("Cannot change the size of block zero", [&](){
      update_sdreopen(cad.name(), 0, 200, false);});

  // (*) Overwrite segment 0 with a larger size.
  must_throw("Cannot write resized segment", [&](){
      update_sdreopen(cad.name(), 0, 2048, false);});

  // (*) Overwrite closed segment with no size change.
  // Support might be added for that if needed.
  must_throw("has already been flushed", [&](){
      update_sdreopen(cad.name(), 512, 2048, false);});
}

#endif // HAVE_SD

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("file.localfilefactory",   test_localfilefactory);
#ifdef HAVE_SD
      register_sd_test("aaa.sdtoken",             test_sdtoken);
      register_sd_test("file.sdfilefactory",      test_sdfilefactory);
      register_sd_test("file.sdfiledelete",       test_sdfiledelete);
      register_sd_test("file.sdreopen",           test_sdreopen);
#endif
    }
  } dummy;
} // namespace for registration
