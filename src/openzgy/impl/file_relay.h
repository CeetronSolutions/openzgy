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

#include "file.h"

namespace InternalZGY {
#if 0
}
#endif

/**
 * Static class for (in some situation) speeding up read from cloud.
 */

/**
 * \brief A class that forwards all IFileBase requests to another instance.
 *
 * \details
 * This class is meant to be extended to add special handling of some
 * of IFileADT members. Used by itself it is a no-op. The class is
 * typically used when just a few methods are to be intercepted and
 * it makes sense to have a default that just passes on the call.
 */
class OPENZGY_TEST_API FileRelayBase : public FileADT
{
private:
  std::shared_ptr<IFileBase> _relay;
  FileRelayBase(const FileRelayBase&) = delete;
  FileRelayBase& operator=(const FileRelayBase&) = delete;

protected:
  const IFileBase& relay() const { return *_relay; };
  IFileBase& relay() { return *_relay; };

public:
  explicit FileRelayBase(std::shared_ptr<IFileBase> relay);
  virtual ~FileRelayBase();
  void deleteFile(const std::string& name, bool missing_ok) const override
  {
    _relay->deleteFile(name, missing_ok);
  }
  std::string  altUrl(const std::string& name) const override
  {
    return _relay->altUrl(name);
  }
  std::string  idToken() const override
  {
    return _relay->idToken();
  }
};

/**
 * \brief An IFileADT that forwards all requests to another instance.
 *
 * \details
 * This class is meant to be extended to add special handling of some
 * of IFileADT members. Used by itself it is a no-op. The class is
 * typically used when just a few methods are to be intercepted and
 * it makes sense to have a default that just passes on the call.
 */
class FileRelay : public FileRelayBase
{
private:
  std::shared_ptr<IFileADT> _relay;
  FileRelay(const FileRelay&) = delete;
  FileRelay& operator=(const FileRelay&) = delete;

protected:
  const IFileADT& relay() const { return *_relay; };
  IFileADT& relay() { return *_relay; };

public:
  explicit FileRelay(std::shared_ptr<IFileADT> relay);
  virtual ~FileRelay();
  void         xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint hint) override {
    _relay->xx_read(data, offset, size, hint);
  }
  void         xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint hint) override {
    _relay->xx_readv(requests, parallel_ok, immutable_ok, transient_ok, hint);
  }
  void         xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint hint) override {
    _relay->xx_write(data, offset, size, hint);
  }
  void         xx_close()         override {
    _relay->xx_close();
  }
  std::int64_t xx_eof()     const override {
    return _relay->xx_eof();
  }
  std::vector<std::int64_t> xx_segments(bool complete) const override
  {
    return _relay->xx_segments(complete);
  }
  bool         xx_iscloud() const override {
    return _relay->xx_iscloud();
  }
};

} // namespace
