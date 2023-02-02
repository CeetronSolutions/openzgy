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
 * \file: perflogger.h
 * \brief Telemetry,
 */

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include <ostream>
#include <atomic>

namespace InternalZGY {
#if 0
}
#endif

class Timer;
class SummaryTimer;

/**
 * Inject code to log throughput, latency, etc.
 *
 * Note the separation of concerns.
 *
 * - PerformanceLogger knows how record the information but not how to
 *   get injected into the production code. The same class can be used
 *   to inject the logger at different levels ane the code can be
 *   shared among projects. Albeit typically using copy/paste.
 *
 * - Timer, PrintingTimer, etc. are low level classes used to store
 *   the results.
 *
 * - FileRelay implements the entire interface where the injection
 *   should happen, making each method just forward to the layer below.
 *   This is a convenient base class for code that only needs to
 *   intercept a few of the methods in the interface. Deriving from
 *   this class makes for cleaner code in the actual interceptor.
 *   The class only knows about the interface that will be used for the
 *   intercept. It may be shared if there are multiple loggers written
 *   for the same interface, In ZGY-Cloud it could have been used by
 *   class TelemetryFileWrapper. But if it works, don't fix it.
 *
 * - FileWithPerformanceLogger glues together the two classes above
 *   and contains the actual code to intercept use of the interface.
 */
class PerformanceLogger
{
private:
  std::shared_ptr<std::ostream> _outfile;
  const std::int64_t _chunksize;
  mutable std::mutex _mutex;
  // Latency a.k.a. round trip time: reported one for each thread.
  std::int64_t _nsamples;
  const double _histmin, _histmax;
  double _statmin, _statmax, _statsum, _statssq;
  std::vector<std::int64_t> _histbins;
  // Throughput: Report transfered bytes vs. elapsed time all threads.
  std::shared_ptr<SummaryTimer> _sumtimer;
  double _sumtimerbeg;
  double _sumtimerend;
  const double _suminterval;
  std::int64_t _sumbytes;
  bool _first;
  int _id;
  const std::string _srcname;
  static std::atomic<int> _last_id;

  PerformanceLogger(const PerformanceLogger&) = delete;
  PerformanceLogger(PerformanceLogger&&) = delete;
  PerformanceLogger& operator=(const PerformanceLogger&) = delete;
  PerformanceLogger& operator=(PerformanceLogger&&) = delete;

public:
  explicit PerformanceLogger(const std::string& outname, std::int64_t chunksize, int hist_bincount, double hist_min, double hist_max, int interval, const std::string& srcname);
  virtual ~PerformanceLogger();
  bool logThisSize(std::int64_t size);
 public:
  void add(const Timer& timer, std::int64_t blocksize);
  // Nonvirtual; might be called from destructor.
  std::string dumpLatency(bool clear);
  std::string dumpThroughput(bool clear);
  void dumpToFile(const std::string& comment);
};

} // namespace
