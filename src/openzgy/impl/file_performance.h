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
 * \file: file_performance.h
 * \brief IFileADT wrapper logging performance statistics.
 */

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include <ostream>

#include "file.h"
#include "file_relay.h"

namespace InternalZGY {
#if 0
}
#endif

class PerformanceLogger;

/**
 * \brief IFileADT wrapper for logging performance statistics.
 */
class FileWithPerformanceLogger : public FileRelay
{
private:
  std::shared_ptr<PerformanceLogger> _recorder;

  FileWithPerformanceLogger(const FileWithPerformanceLogger&) = delete;
  FileWithPerformanceLogger(FileWithPerformanceLogger&&) = delete;
  FileWithPerformanceLogger& operator=(const FileWithPerformanceLogger&) = delete;
  FileWithPerformanceLogger& operator=(FileWithPerformanceLogger&&) = delete;

public:
  explicit FileWithPerformanceLogger(std::shared_ptr<IFileADT> relay, const std::string& outname, std::int64_t chunksize, int hist_bincount, double hist_min, double hist_max, int interval, const std::string& srcname);
  virtual ~FileWithPerformanceLogger();

  // Intercept
  void         xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint) override;
  void         xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint) override;
  void         xx_close() override;
  static std::shared_ptr<IFileADT> inject(std::shared_ptr<IFileADT> file, const std::string& srcname);
};

} // namespace
