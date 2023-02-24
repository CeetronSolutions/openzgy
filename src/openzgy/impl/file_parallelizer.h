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

#include "file_relay.h"

namespace InternalZGY {
#if 0
}
#endif

/**
 * Static class for (in some situation) speeding up read from cloud.
 * See SeismicStoreFile::xx_readv() for some explanation.
 */
class FileParallelizer : public FileRelay
{
  std::int64_t _cputhreads;
public:
  FileParallelizer(std::shared_ptr<IFileADT> relay, std::int64_t cputhreads);
  FileParallelizer(const FileParallelizer&)            = delete;
  FileParallelizer& operator=(const FileParallelizer&) = delete;
  virtual ~FileParallelizer();
  void xx_readv(const ReadList& requests,
                        bool parallel_ok,
                        bool immutable_ok,
                        bool transient_ok,
                        UsageHint hint) override;
  static std::shared_ptr<IFileADT> inject(std::shared_ptr<IFileADT> file, std::int64_t cputhreads);
};

} // namespace
