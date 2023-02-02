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

/**
 * \file: file_smallcache.cpp
 * \brief IFileADT wrapper adding a naive caching of the ZGY header area.
 */
#include "file_smallcache.h"
#include "exception.h"
#include <string.h>

namespace InternalZGY {
#if 0
}
#endif

FileWithSmallCache::FileWithSmallCache(std::shared_ptr<IFileADT> relay, std::int64_t size)
  : FileRelayBase(relay)
  , _relay(relay)
  , _cache()
  , _cachesize(size)
{
}

FileWithSmallCache::~FileWithSmallCache()
{
}

void
FileWithSmallCache::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  if (offset + size > _cachesize || _cachesize <= 0) {
    //printf("Cache miss offset %llx size %llx\n", offset, size);
    _relay->xx_read(data, offset, size, usagehint);
    return;
  }
  if (!_cache) {
    std::int64_t filesize = _relay->xx_eof();
    if (filesize < _cachesize) {
      _cachesize = filesize;
      //printf("Cache entire file %lld bytes (0x%llx)\n", _cachesize, _cachesize);
    }
    if (_cachesize > 0) {
      //printf("Cache load offset %llx size %llx\n", 0, _cachesize);
      std::unique_ptr<char[]> newcache(new char[_cachesize]);
      _relay->xx_read(newcache.get(), 0, _cachesize, UsageHint::Unknown);
      _cache.swap(newcache);
    }
  }
  // Repeat the test in case _cachesize changed.
  if (offset + size > _cachesize || _cachesize <= 0) {
    //printf("Cache miss #2 offset %llx size %llx\n", offset, size);
    _relay->xx_read(data, offset, size, usagehint);
    return;
  }
  //printf("Cache hit offset %llx size %llx\n", offset, size);
  memcpy(data, _cache.get() + offset, size);
}

void
FileWithSmallCache::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
  // When OpenZGY reads header data it currently doesn't use xx_readv.
  // So just treat this as a cache miss.
  _relay->xx_readv(requests, parallel_ok, immutable_ok, transient_ok, usagehint);
}

void
FileWithSmallCache::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown)
{
  throw OpenZGY::Errors::ZgyInternalError("The small cache doesn't allow writes");
}

void
FileWithSmallCache::xx_close()
{
  _relay->xx_close();
}

std::int64_t
FileWithSmallCache::xx_eof() const
{
  return _relay->xx_eof();
}

std::vector<std::int64_t>
FileWithSmallCache::xx_segments(bool complete) const
{
  return _relay->xx_segments(complete);
}

bool
FileWithSmallCache::xx_iscloud() const
{
  return _relay->xx_iscloud();
}

} // namespace
