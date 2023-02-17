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

#ifdef HAVE_SD // Rest of file

#include "test_sdutils.h"
#include "../impl/environment.h"
#define TEST_NO_MAIN
#define NOMINMAX // cutest includes Windows.h. Ugh!
#include "cutest.h"

#ifndef _WIN32
#include <SDManager.h>
#include <SDGenericDataset.h>
#else
#include <SDAPI/SDManager.h>
#include <SDAPI/SDGenericDataset.h>
#endif

#include <stdexcept>

namespace InternalZGY {
  extern void OPENZGY_TEST_API hack_setAuthProvider(seismicdrive::SDManager*,const std::string&);
}

namespace Test_Utils {
#if 0
}
#endif

static std::shared_ptr<seismicdrive::SDManager>
get_manager()
{
  using InternalZGY::Environment;
  const std::string sdurl{Environment::getStringEnv("OPENZGY_SDURL")};
  const std::string sdapikey{Environment::getStringEnv("OPENZGY_SDAPIKEY")};
  const std::string sdtoken{Environment::getStringEnv("OPENZGY_TOKEN")};
  if (sdurl.empty() || sdapikey.empty() || sdtoken.empty())
    throw std::runtime_error("Missing $OPENZGY_{SDURL,SDAPIKEY,TOKEN}");
  std::shared_ptr<seismicdrive::SDManager> sdmanager
    (new seismicdrive::SDManager(sdurl, sdapikey));
  InternalZGY::hack_setAuthProvider(sdmanager.get(), sdtoken);
  return sdmanager;
}

/**
 * Return a functor suitable for use in ctx.sdtokencb().
 * New tokens are generated from a SDManager created in the default
 * manner with credentials in $OPENZGY_TOKEN. If those credentials
 * support refreshing then the callback does as as well.
 *
 * Using one manager to provide credentials from another might not
 * be technically legal, because the callback should not make
 * calls back to SDAPI. In this case it will probably work,
 * but don't use this trick in production.
 *
 * getIDToken() is thread safe so I shouldn't need additional
 * locking here.
 */
std::function<std::string()>
get_token_callback()
{
  std::shared_ptr<seismicdrive::SDManager> mgr = get_manager();
  return [mgr]() {
           std::string token = mgr->getIDToken();
           static std::string bearer("Bearer ");
           return token.rfind(bearer) ? token : token.substr(bearer.size());
         };
}

/**
 * Copy a file on seismic store to another file on seismic store.
 * This uses direct calls to SDAPI only; there is no OpenZGY involved.
 * To be used for setting up unit test data. Caveat: No legaltag
 * is provided, so the destination subproject must have a default set.
 */
void
copy_sd_to_sd(const std::string& srcname, const std::string& dstname)
{
  std::shared_ptr<seismicdrive::SDManager> manager = get_manager();
  seismicdrive::SDGenericDataset src(manager.get(), srcname);
  src.open(seismicdrive::SDDatasetDisposition::READ_ONLY);
  const std::uint64_t nblocks = src.getBlockNum();
  std::vector<std::string> names;
  for (std::uint64_t ii = 0; ii < nblocks; ++ii)
    names.push_back(std::to_string(ii));
  const std::vector<long long> sizearray = src.getBlocksSize(names);
  for (long long segsize : sizearray)
    if (segsize < 0)
      throw std::runtime_error("Segment size must be > 0");

  seismicdrive::SDGenericDataset dst(manager.get(), dstname);
  dst.open(seismicdrive::SDDatasetDisposition::OVERWRITE);

  for (std::uint64_t segnum = 0; segnum < nblocks; ++segnum) {
    //std::cout << "copy " << segnum << " size " << sizearray[segnum] << "\n";
    std::unique_ptr<char[]> data(new char[sizearray[segnum]]);
    src.readBlock((int)segnum, data.get(), 0, (size_t)sizearray[segnum]);
    dst.writeBlock((int)segnum, data.get(), (size_t)sizearray[segnum]);
  }
  dst.close();
  src.close();
}

std::vector<std::int64_t>
get_segsizes(const std::string& name)
{
  std::shared_ptr<seismicdrive::SDManager> manager = get_manager();
  seismicdrive::SDGenericDataset src(manager.get(), name);
  src.open(seismicdrive::SDDatasetDisposition::READ_ONLY);
  const std::uint64_t nblocks = src.getBlockNum();
  std::vector<std::string> names;
  for (std::uint64_t ii = 0; ii < nblocks; ++ii)
    names.push_back(std::to_string(ii));
  const std::vector<long long> sizearray = src.getBlocksSize(names);
  src.close();
  std::vector<std::int64_t> result;
  for (auto it : sizearray)
    result.push_back(it);
  return result;
}

} // namespace Test_Utils

#endif // HAVE_SD
