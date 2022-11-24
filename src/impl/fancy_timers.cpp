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

#include "fancy_timers.h"
#include "environment.h"
#include <iostream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Create a new instance. The base class is constructed with the default
 * csv=false but that makes no difference for us because print() is
 * redefined here. So we decide for ourselves whether to output csv
 * or not.
 *
 * In the simplest case, instances are statically constructed and
 * print their result on program exit.
 */
SummaryPrintingTimerEx::SummaryPrintingTimerEx(const char *name)
  : SummaryPrintingTimer(name)
  , bytes_read_(0)
  , bytes_written_(0)
{
}

SummaryPrintingTimerEx::~SummaryPrintingTimerEx()
{
  print();
}

/**
 * \brief Convert a long integer to a human readable string.
 * \details
 * Use the largest suffix (TB, MB, etc.) that still allows the number
 * to be displayed without any decimal point. So, 2*(1024^3) will be
 * disokated as  "2 GB" while that number plus 1024 will display
 * returns the size in kilobytes. I.e. "2097153 KB".
 */
std::string
SummaryPrintingTimerEx::niceSize(std::int64_t n)
{
  static struct {std::int64_t factor; const char *unit;} lookup[]{
    {1024LL*1024LL*1024LL*1024LL, " TB"},
    {1024*1024*1024, " GB"},
    {1024*1024, " MB"},
    {1024, " KB"},
    {1, " bytes"},
  };
  std::string neg(n<0?"-":"");
  n = std::abs(n);
  for (const auto& it : lookup)
    if (n >= it.factor && (n % it.factor) == 0)
      return neg + std::to_string(n / it.factor) + std::string(it.unit);
  return neg + std::to_string(n) + " bytes";
}

void
SummaryPrintingTimerEx::printToFile(std::ostream& outstream, bool csv, bool clear)
{
  if (getCount() != 0) {
    std::string msg(csv ? getCSV() : getValue(true, true));
    if (!msg.empty() && msg.back() == '\n')
      msg = msg.substr(0, msg.size()-1);
    if (csv)
      outstream << msg
                << "," << bytes_read_.load()
                << "," << bytes_written_.load()
                << std::endl;
    else
      outstream << msg
                << ", R: " << niceSize(bytes_read_.load())
                << ", W: " << niceSize(bytes_written_.load())
                << std::endl;
  }
  if (clear)
    reset();
}

void
SummaryPrintingTimerEx::print()
{
  printToFile(std::cerr, isCSVEnabled(), true);
}

void
SummaryPrintingTimerEx::addBytesRead(std::int64_t nbytes)
{
  bytes_read_.fetch_add(nbytes);
}

void SummaryPrintingTimerEx::addBytesWritten(std::int64_t nbytes) {
  bytes_written_.fetch_add(nbytes);
}

bool
SummaryPrintingTimerEx::isCSVEnabled()
{
  static int enable = Environment::getNumericEnv("OPENZGY_TIMERS", 0);
  return enable > 1;
}

/**
 * The only difference between a SimpleTimer and a SimpleTimerEx
 * is that the latter is enabled or disabled using the environment
 * variable OPENZGY_TIMERS. So it becomes specific to this library.
 */
SimpleTimerEx::SimpleTimerEx(SummaryTimer& owner)
  : SimpleTimer(owner, isTimerEnabled())
{
};

bool
SimpleTimerEx::isTimerEnabled()
{
  static int enable = Environment::getNumericEnv("OPENZGY_TIMERS", 0);
  return enable > 0;
}

} // namespace
