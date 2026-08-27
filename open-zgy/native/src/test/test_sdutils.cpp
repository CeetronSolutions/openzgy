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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#ifndef _WIN32
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#endif

namespace InternalZGY {
  extern std::shared_ptr<seismicdrive::SDManager>
  OPENZGY_TEST_API
  hack_createAuthorizedManager(
     const std::string& sdurl,
     const std::string& sdkey,
     const std::string& token);
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
  return InternalZGY::hack_createAuthorizedManager(sdurl, sdapikey, sdtoken);
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

static std::string
replchar(const std::string& in, char from, char to)
{
  std::string str(in);
  std::transform(str.begin(), str.end(), str.begin(),
                 [from,to](char ch){ return ch == from ? to : ch; }
                 );
  return str;
}

std::vector<std::pair<std::string,std::int64_t>>
get_segments(const std::string& name)
{
  std::vector<std::pair<std::string,std::int64_t>> result;
  std::shared_ptr<seismicdrive::SDManager> manager = get_manager();
  seismicdrive::SDGenericDataset src(manager.get(), name);
  src.open(seismicdrive::SDDatasetDisposition::READ_ONLY);
  std::int64_t tally_count{0}, tally_size{0};
  auto it = src.getIterator(true);
  while (it.hasNext()) {
    auto info = it.next();
    tally_count += 1;
    tally_size += info.getSize();
    result.push_back(std::make_pair(info.getName(), info.getSize()));
  }
  std::int64_t expect_count = (std::int64_t)src.getBlockNum();
  std::int64_t expect_size = (std::int64_t)src.getSize();
  std::cout << "expect " << expect_count << " " << expect_size
            << " actual " << tally_count << " " << tally_size
            << std::endl;
  src.close();
  return result;
}

#ifndef _WIN32
/**
 * Naive directory listing, including file sizes.
 */
std::vector<std::pair<std::string,std::int64_t>>
ls_local_files(const std::string& dir, const std::string& prefix)
{
  std::vector<std::pair<std::string,std::int64_t>> segments;
  DIR* handle = opendir(dir.c_str());
  if (!handle)
    throw std::runtime_error("Cannot open folder \"" + dir + "\"");
  for (;;) {
    errno = 0;
    const struct dirent *entry = readdir(handle);
    if (entry == nullptr) {
      if (errno != 0)
        throw std::runtime_error("Cannot read folder \"" + dir + "\"");
      break;
    }
    std::string basename(entry->d_name);
    std::string fullpath(dir + "/" + basename);
    if (basename.empty() || basename.substr(0, 1) == ".")
      continue;
    if (!prefix.empty() && basename.substr(0, prefix.size()) != prefix)
      continue;
    struct stat statbuf{};
    if (stat(fullpath.c_str(), &statbuf) < 0)
      throw std::runtime_error("Cannot stat \"" + fullpath + "\"");
    segments.push_back(std::make_pair(basename, statbuf.st_size));
  }
  if (closedir(handle) < 0)
    throw std::runtime_error("Cannot close folder \"" + dir + "\"");
  std::sort(segments.begin(), segments.end());
  return segments;
}
#else
std::vector<std::pair<std::string, std::int64_t>>
ls_local_files(const std::string&, const std::string&)
{
  return std::vector<std::pair<std::string, std::int64_t>>();
}
#endif

/**
 * Copy a data set on seismic store to a local folder, with one file
 * per segment. To be used for setting up unit test data. Caveat: No
 * legaltag is provided, so the destination subproject must have a
 * default set.
 */
void
copy_sd_to_folder(const std::string& srcname, const std::string& dstname, int verbose)
{
  if (verbose)
    std::cout << "\nCopy \"" << srcname
              << "\" to \"" << dstname << "\"" << std::endl;
  if (!ls_local_files(dstname, "seg.").empty())
    throw std::runtime_error("Target \"" + dstname + "\" is not empty");
  std::vector<std::pair<std::string,std::int64_t>> segments;
  std::shared_ptr<seismicdrive::SDManager> manager = get_manager();
  seismicdrive::SDGenericDataset src(manager.get(), srcname);
  src.open(seismicdrive::SDDatasetDisposition::READ_ONLY);
  {
    auto it = src.getIterator(true);
    while (it.hasNext()) {
      auto info = it.next();
      segments.push_back(std::make_pair(info.getName(), info.getSize()));
    }
  }
  //std::sort(segments.begin(), segments.end());
  // TODO check for duplicates?
  for (const auto& it : segments) {
    // If sehment name contains '%', we are toast.
    std::string quoted = replchar(it.first, '/', '%');
    std::string fname = dstname + "/seg." + quoted;
    if (verbose)
      std::cout << "Copy " << std::setw(10) << it.second
                << " \"" << it.first << "\"" << std::endl;
    // read
    if (it.second < 0)
      throw std::runtime_error("Bad size in segment " + it.first + "\"");
    std::unique_ptr<char[]> data(new char[it.second]);
    src.readBlock(it.first, data.get(), 0, (size_t)it.second);
    // write
    std::ofstream fhandle(fname, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
    if (!fhandle.good())
      throw std::runtime_error("Error creating local file \"" + fname + "\"");
    fhandle.write(data.get(), it.second);
    fhandle.close();
  }
  src.close();
}

void
copy_folder_to_sd(const std::string& srcname, const std::string& dstname, int verbose)
{
  if (verbose)
    std::cout << "\nCopy \"" << srcname
              << "\" to \"" << dstname << "\"" << std::endl;
  const std::string prefix("seg.");
  const auto segments = ls_local_files(srcname, prefix);
  if (segments.empty())
    throw std::runtime_error("No \"" + prefix + "*\" files found in \"" + srcname + "\"");
  std::shared_ptr<seismicdrive::SDManager> manager = get_manager();
  seismicdrive::SDGenericDataset dst(manager.get(), dstname);
  dst.open(seismicdrive::SDDatasetDisposition::CREATE);
  for (const auto& it : segments) {
    std::string unquoted = replchar(it.first, '%', '/');
    const std::string segname = unquoted.substr(prefix.size());
    const std::string fullpath(srcname + "/" + it.first);
    const std::int64_t segsize = it.second;
    if (verbose)
      std::cout << "Copy " << std::setw(10) << segsize
                << " \"" << segname << "\"" << std::endl;
    // read
    if (segsize < 0)
      throw std::runtime_error("Bad size in segment " + it.first + "\"");
    std::unique_ptr<char[]> data(new char[segsize]);
    std::ifstream fhandle(fullpath, std::ios_base::in|std::ios_base::binary);
    if (!fhandle.good())
      throw std::runtime_error("Error reading local file \"" + fullpath + "\"");
    fhandle.read(data.get(), segsize);
    if (!fhandle.good())
      throw std::runtime_error("Error reading local \"" + fullpath + "\"");
    // write
    dst.writeBlock(segname, data.get(), (std::size_t)segsize, false);
    fhandle.close();
  }
  dst.setReadonlyMode(true);
  dst.close();
}

} // namespace Test_Utils

#endif // HAVE_SD
