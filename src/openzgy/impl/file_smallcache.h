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

/**
 * \file: file_smallcache.h
 * \brief IFileADT wrapper adding a naive caching of the ZGY header area.
 */

#include <cstdint>
#include <vector>
#include <string>
#include <memory>

#include "../declspec.h"
#include "file.h"
#include "file_relay.h"

namespace InternalZGY {
#if 0
}
#endif

/**
 * \brief IFileADT wrapper adding a naive caching of the ZGY header area.
 *
 * \details
 *
 * Problem: Header information is read from the file back-end one header
 * at a time. This is because the size of each of the headers might not
 * be known until the previous headers have been read.
 *
 * Solution: Use a short lived cache, guessing how large it should be.
 *
 * The cache should only be used while opening a file for read. It
 * makes no sense to keep it afterwards since OpenZGY reads all the
 * metadata up front. This also means the code doesn't need to worry
 * about the cache going stale.
 *
 * The cache knows neither the brick size nor the total size of the
 * header area due to a chicken and egg situation. This means it
 * doesn't know where the boundary between segment 0 and segment 1 is.
 * Segment 0 doesn't need to be larger than what is needed for the
 * headers.
 *
 * Choosing to cache 256 KB, which is the size of a single brick with
 * default size and int8 storage. Even if the code later realizes that
 * the header was larger the cache won't expand. Note that ZgyPublic
 * and OpenZGY both pad the header area to a multiple of the brick
 * size. So most files will have at least 256 KB.
 *
 * If the actual segment 0 size is less (meaning a small file with
 * small brick size) then the cache will needlessly read part of
 * segment 1 as well. Because it doesn't realize that it is crossing a
 * segment boundary. I *could* avoid that by querying the file backend
 * but that complicates the code to account for something that is very
 * unlikely to happen.
 *
 * If the actual segment 0 size is larger than 256 KB this means the
 * file is huge, containing at least 30,000 bricks. The code will then
 * get one or more cache misses and may end up issuing up to 6 reads
 * to get all the headers.
 *
 * Uncompressed files v1 will also see some cache misses because some
 * of the headers are stored at the end of the file.
 *
 * The cache size also needs to be chosen so that reading too much
 * won't waste time. For cloud access anything less than 2 MB will
 * probably take the same time. So that is not a concern here.
 *
 * Thread safety: Not safe for writes. Writes will currently raise an
 * exception because users of this wrapper won't need them. Read
 * methods do not have race conditions.
 */
class OPENZGY_TEST_API FileWithSmallCache : public FileRelayBase
{
private:
  std::shared_ptr<IFileADT> _relay;
  std::unique_ptr<char[]> _cache;
  std::int64_t _cachesize;
  FileWithSmallCache(const FileWithSmallCache&) = delete;
  FileWithSmallCache& operator=(const FileWithSmallCache&) = delete;

public:
  explicit FileWithSmallCache(std::shared_ptr<IFileADT> relay, std::int64_t size);
  virtual ~FileWithSmallCache();
  void         xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint) override;
  void         xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint) override;
  void         xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint) override;
  void         xx_close()         override;
  std::int64_t xx_eof()     const override;
  std::vector<std::int64_t> xx_segments(bool complete) const override;
  bool         xx_iscloud() const override;
};

} // namespace
