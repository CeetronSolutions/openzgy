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

#include "file_performance.h"
#include "perflogger.h"
#include "timer.h"
#include "environment.h"
#include "exception.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <atomic>
#include <algorithm>

namespace InternalZGY {
#if 0
}
#endif

FileWithPerformanceLogger::FileWithPerformanceLogger(std::shared_ptr<IFileADT> relay, const std::string& outname, std::int64_t chunksize, int hist_bincount, double hist_min, double hist_max, int interval, const std::string& srcname)
  : FileRelay(relay)
  , _recorder(new PerformanceLogger
              (outname, chunksize, hist_bincount, hist_min, hist_max, interval, srcname))

{
}

FileWithPerformanceLogger::~FileWithPerformanceLogger()
{
  _recorder->dumpToFile("destructed");
}

void
FileWithPerformanceLogger::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  if (_recorder->logThisSize(size)) {
    Timer timer;
    relay().xx_read(data, offset, size, usagehint);
    timer.stop();
    _recorder->add(timer, size);
  }
  else {
    relay().xx_read(data, offset, size, usagehint);
  }
}

void
FileWithPerformanceLogger::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
  // This is only effective in the fully random access case,
  // To include the "consolidate" case I might need to implenent it
  // inside SeismicStoreFile:. Or make the entire consolidate logic
  // into a separate module that can be chained. Either way I may need
  // multiple histograms to cover different brick sizes.
  if ((requests.size() == 1 && _recorder->logThisSize(requests.front().size)) ||
      _recorder->logThisSize(-1)) {
    Timer timer;
    relay().xx_readv(requests, parallel_ok, immutable_ok, transient_ok, usagehint);
    timer.stop();
    std::int64_t size = 0;
    for (const ReadRequest& it : requests)
      size += it.size;
    _recorder->add(timer, size);
  }
  else {
    relay().xx_readv(requests, parallel_ok, immutable_ok, transient_ok, usagehint);
  }
}

void
FileWithPerformanceLogger::xx_close()
{
  relay().xx_close();
  _recorder->dumpToFile("closed");
}

/**
 * Inject a telemetry module if enabled by environment variables.
 * If not enabled the telemetry code has zero impact on the system.
 *
 *    OPENZGY_MEASURE_KB       = brick size to monitor, or -1 for all reads.
 *    OPENZGY_MEASURE_LOGFILE  = optionally write to this file
 *    OPENZGY_MEASURE_BINS     = bins in histogram (default 251)
 *    OPENZGY_MEASURE_TIME     = highest latency in ms (default 500)
 *    OPENZGY_MEASURE_INTERVAL = periodic report of throughput and latency
 *
 * There are three distinct reports being output:
 * - Periodic latency and throughput, enable by MEASURE_INTERVAL, tags CSV0,CSV8
 * - Latency statistics for the entire file, tags CSV1,CSV2
 * - Latency histogram for the entire file, tags CSV3,CSV4,CSV5,CSV6
 *   and can be fine tuned by MEASURE_BINS and MEASURE_TIME.
 */
std::shared_ptr<IFileADT>
FileWithPerformanceLogger::inject(std::shared_ptr<IFileADT> file, const std::string& srcname)
{
  int target = Environment::getNumericEnv("OPENZGY_MEASURE_KB", 0);
  if (target != 0) {
    int bincount = Environment::getNumericEnv("OPENZGY_MEASURE_BINS", 251);
    int maxtime = Environment::getNumericEnv("OPENZGY_MEASURE_TIME", 500);
    int interval = Environment::getNumericEnv("OPENZGY_MEASURE_INTERVAL", 0);
    std::string filename = Environment::getStringEnv("OPENZGY_MEASURE_LOGFILE");
    std::shared_ptr<std::ostream> out;
    file = std::shared_ptr<IFileADT>(new FileWithPerformanceLogger(file, filename, target*1024, bincount, 0.0, maxtime, interval, srcname));
  }
  return file;
}

} // namespace
