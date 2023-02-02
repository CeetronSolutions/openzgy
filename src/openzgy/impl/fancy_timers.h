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

#include "timer.h"
#include <atomic>
#include <cstdint>
#include <string>
#include <ostream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file fancy-timers.h
 * \brief Easy to use performance measurement.
 * \details Adds bells and whistles specifically for this project.
 * The base class is deliberately kept simple and with few
 * dependencies so it can be copy/pasted between projects.
 */


/**
 * \brief SummaryTimer that prints its result when going out of scope.
 * \details This is a private extension that also knows about total
 * number of bytes transfered, so it can report the actuar throughput.
 */
class OPENZGY_TEST_API SummaryPrintingTimerEx : public SummaryPrintingTimer
{
  std::atomic<std::int64_t> bytes_read_;
  std::atomic<std::int64_t> bytes_written_;
 public:
  explicit SummaryPrintingTimerEx(const char *name);
  virtual ~SummaryPrintingTimerEx();
  static std::string niceSize(std::int64_t n);
  virtual void print();
  void printToFile(std::ostream& os, bool csv, bool clear);
  void addBytesRead(std::int64_t nbytes);
  void addBytesWritten(std::int64_t nbytes);
  static bool isCSVEnabled();
};

/**
 * The only difference between a SimpleTimer and a SimpleTimerEx
 * is that the latter is enabled or disabled using the environment
 * variable OPENZGY_TIMERS. So it becomes specific to this library.
 */
class OPENZGY_TEST_API SimpleTimerEx : public SimpleTimer
{
 public:
  explicit SimpleTimerEx(SummaryTimer& owner);
  static bool isTimerEnabled();
};

} // namespace
