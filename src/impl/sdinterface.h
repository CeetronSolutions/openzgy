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

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>

namespace seismicdrive {
  class SDGenericDataset;
  class SDManager;
  struct ExponentialRetryBackoffPolicy;
  enum class HttpConnectionLink;
  enum class SDDatasetDisposition;
}

namespace InternalZGY {
#if 0
}
#endif

/**
 * Wrap seismicdrive::SDGenericDataset in a pure virtual interface.
 * This is done to make it simpler to insert optional modules.
 * It can also be helpful if SDAPI needs to be mocked.
 *
 * Currently only the methods needed by OpenZGY are included.
 *
 * It is possible to add a relay() back door that allows callers to
 * invoke any method, even if not wrapped, by bypassing one layer. Or
 * use instance.relay().relay() to bypass two layers, etc. That would
 * however be dangerous and make this interface class harder to use
 * correctly. The application would need to know which methods are
 * safe to call. If the interface is used for mocking the answer might
 * simply be "none".
 *
 * The interface might later be extended to hold all methods not
 * already deprecated. This makes it more generally useful.
 *
 * The code might also be moved inside SDAPI. BUT, note the
 * caveat that a class with virtual methods is more difficult
 * to keep binary compatible.
 *
 * FUTURE: The existing SDGenericDataset has inconsistent usage of signed
 * vs unsigned, long long vs. std::int64_t, std::int32_t vs std::int64_t,
 * and const vs. mutable. This interface doesn't technically need to match
 * what SDGenericDataset uses. So, very soon it won't. I prefer const for
 * anything semantically const even when SDAPI does caching behind the
 * scenes. Because applications don't need to know that. I also prefer
 * using std::int64_t for all integral types unless there is a really
 * good reason not to.
 *
 * Original SDAPI rules for integer:
 *
 *  - Block numbers, as a convenience instead of strings, are signed int.
 *    Negative numbers are technically legal but when numerical block
 *    names are used the convention is to name them consecutively from 0
 *    up to getBlockNum(), which can actually be larger then an int.
 *    SUGGESTION: Could be fixed by declaring the block numbers as uint64_t
 *    and then have the wrapper use std::to_string() and call the
 *    overload that expects a name.
 *
 *  - Block offset/size arguments are size_t which is normally unsigned 64-bit.
 *
 *  - BUT, block sizes returned from functions are signed long long.
 *    A negative result might be used to signal an error.
 *    SUGGESTION: It would have better to at least declare these as ssize_t,
 *    to not mix old fashioned int, long, etc. with higher level types.
 *    There is still a signed/unsigned inconsistency between creating
 *    a block and getting its size (or size of all blocks).
 *
 *  - The total number of blocks is returned as uint64_t defined in stdint.h
 */
class ISDGenericDataset
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;
  virtual ~ISDGenericDataset();
  virtual void          open(seismicdrive::SDDatasetDisposition disposition, const std::unordered_map<std::string, std::string> &args) = 0;
  virtual void          close() = 0;
  virtual std::string   getConsistencyID() const = 0;
  virtual std::string   getSerializedContext() const = 0;
  virtual bool          getReadonlyMode() const = 0;
  virtual void          setReadonlyMode(bool readonly) = 0;
  virtual void          setExponentialRetryBackoffPolicy(const seismicdrive::ExponentialRetryBackoffPolicy *policy, const seismicdrive::HttpConnectionLink link) = 0;
  virtual void          readBlock(int blocknum, char *data, std::size_t offset, std::size_t numBytes) const = 0;
  virtual void          writeBlock(int blocknum, const char *data, std::size_t len, bool check_and_overwrite = false) = 0;
  virtual void          deleteBlock(const std::string &blockName) = 0;
  virtual long long     getSize() const = 0;
  virtual std::uint64_t getBlockNum() const = 0;
  virtual long long     getBlockSize(int blocknum) const = 0;
  virtual std::vector<long long>
    getBlockSizes(const std::vector<std::string> &blockNames) const = 0;

  // Dangerous! Only enable if really needed. See class documentation.
  //virtual seismicdrive::SDGenericDataset& relay() = 0; //{ return *relay_; }
  //virtual const seismicdrive::SDGenericDataset& relay() const = 0; //{ return *relay_; }

  // Doesn't belong in the pure interface but is more convenient.
  static std::shared_ptr<ISDGenericDataset> createPlainInstance(seismicdrive::SDManager *manager, const std::string &filename, const bool log = false);
  static std::shared_ptr<ISDGenericDataset> injectAsyncInstance(const std::shared_ptr<ISDGenericDataset>& relay, std::int64_t highwater, const LoggerFn& logger);
};

} // namespace
