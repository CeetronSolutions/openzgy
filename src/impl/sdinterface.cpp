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

#include "sdinterface.h"

#ifdef HAVE_SD // Most of the file

#ifndef _WIN32
#include <SDManager.h>
#include <SDGenericDataset.h>
#else
#include <SDAPI/SDManager.h>
#include <SDAPI/SDGenericDataset.h>
#endif

namespace InternalZGY {
#if 0
}
#endif

class SDGenericDatasetPlain : public ISDGenericDataset
{
public:
  SDGenericDatasetPlain(std::shared_ptr<seismicdrive::SDGenericDataset> sgds)
    : relay_(sgds)
  {
  }

  virtual ~SDGenericDatasetPlain();

  virtual void open(
       seismicdrive::SDDatasetDisposition disposition,
       const std::unordered_map<std::string, std::string> &args) override
  {
    relay().open(disposition, args);
  }

  virtual void close() override
  {
    return relay().close();
  }

  virtual std::string getConsistencyID() const override
  {
    return mutableRelay().getConsistencyID();
  }

  virtual std::string getSerializedContext() const override
  {
    return mutableRelay().getSerializedContext();
  }

  virtual bool getReadonlyMode() const override
  {
    return relay().getReadonlyMode();
  }

  virtual void setReadonlyMode(bool readonly) override
  {
    return relay().setReadonlyMode(readonly);
  }

  virtual void setExponentialRetryBackoffPolicy(
       const seismicdrive::ExponentialRetryBackoffPolicy *policy,
       const seismicdrive::HttpConnectionLink link) override
  {
    return relay().setExponentialRetryBackoffPolicy(policy, link);
  }

  virtual void readBlock(int blocknum, char *data, std::size_t offset, std::size_t numBytes) const override
  {
    return mutableRelay().readBlock(blocknum, data, offset, numBytes);
  }

  virtual void writeBlock(int blocknum, const char *data, std::size_t len, bool check_and_overwrite) override
  {
    return relay().writeBlock(blocknum, data, len, check_and_overwrite);
  }

  virtual void deleteBlock(const std::string &blockName) override
  {
    return relay().deleteBlock(blockName);
  }

  virtual long long getSize() const override
  {
    return mutableRelay().getSize();
  }

  virtual std::uint64_t getBlockNum() const override
  {
    return mutableRelay().getBlockNum();
  }

  virtual long long getBlockSize(int blocknum) const override
  {
    return mutableRelay().getBlockSize(blocknum);
  }

  virtual std::vector<long long>
  getBlockSizes(const std::vector<std::string> &blockNames) const override
  {
    // TODO use getBlockSizes() once SDAPI has been upgraded for windows.
    return mutableRelay().getBlocksSize(blockNames);
  }

private:
  seismicdrive::SDGenericDataset& relay() { return *relay_; }
  const seismicdrive::SDGenericDataset& relay() const { return *relay_; }
  // to be used only when the official SDAPI has declared a method
  // as non-const when it ought to have been const.
  seismicdrive::SDGenericDataset& mutableRelay() const {
    return const_cast<seismicdrive::SDGenericDataset&>(*relay_);
  }

private:
  std::shared_ptr<seismicdrive::SDGenericDataset> relay_;
};

/**
 * The destructor is a no-op, but the compiler likes to habe at least
 * one virtual non-inlined method so it knows where to put the vtbl.
 */
SDGenericDatasetPlain::~SDGenericDatasetPlain()
{
}

/**
 * The destructor is a no-op, but the compiler likes to habe at least
 * one virtual non-inlined method so it knows where to put the vtbl.
 */
ISDGenericDataset::~ISDGenericDataset()
{
}

std::shared_ptr<ISDGenericDataset>
ISDGenericDataset::createPlainInstance(seismicdrive::SDManager *manager, const std::string &filename, const bool log)
{
  auto relay =
    std::make_shared<seismicdrive::SDGenericDataset>(manager, filename, log);
  return std::make_shared<SDGenericDatasetPlain>(relay);
}

} // namespace

#endif // HAVE_SD
