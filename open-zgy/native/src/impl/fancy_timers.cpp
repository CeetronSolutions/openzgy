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
#include <fstream>

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
    else if (!bytes_read_.load() && !bytes_written_.load())
      outstream << msg << std::endl;
    else
      outstream << msg
                << ", R: " << niceSize(bytes_read_.load())
                << ", W: " << niceSize(bytes_written_.load())
                << std::endl;
  }
  if (clear)
    reset();
}

/**
 * Output a single line of timer information, by default to std::cerr
 * but might be redirected to a file.
 *
 * CAVEATS when writing to a file:
 *
 * The output file will be opened and closed once per line. That can
 * easily be changed by making "os" a static variable. But in that
 * case, closing the file happens in a static destructor. And some
 * timers are printed from static destructors. Be very careful with
 * teardown order if changing this.
 *
 * Even with the mitigation of open/close for each line, there is some
 * risk of crash on exit. Possibly the std::ostream mechanism itself
 * might be shut down, making it illegal to open the log file. The
 * code should perhaps use lower level I/O than what STL provides.
 *
 * Arguably it might have been better to redirect timer output to the
 * standard loggers. But keeping the timer mechanism separate from
 * other logging might be cleaner. Since loggers are typically owned
 * by some data structure, not global. Also, the issues with ordering
 * of static destructors would likely get worse.
 *
 * The frequent open/close could skew the timing results when there
 * are nested timers and the outer timer includes the time to print
 * the inner timers.
 *
 * The mechanism to include an id in the file name, to output one log
 * file for each file handle to be measured, is not in use. The reason
 * is that several of the timers don't know which file is being
 * acessed.
 *
 * In the future, the code might offer output to ::OutputDebugStringA()
 * on Windows. Yagni until somebody says they need it. If they do,
 * the logger mechanism probsably needs that feature as well.
 */
void
SummaryPrintingTimerEx::print()
{
  std::shared_ptr<std::ostream> os = getVerboseFileFromEnv(0);
  printToFile(*os, isCSVEnabled(), true);
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

/**
 * Copy-paste from LoggerBaseImpl::getVerboseFileFromEnv() and
 * ApiRecorder::getVerboseFileFromEnv(). Arguably it might have been
 * better to redirect timer output to the standard loggers. But
 * keeping the timer mechanism separate from other logging might be
 * cleaner. Since loggers are typically owned by some data structure,
 * not global.
 */
std::shared_ptr<std::ostream>
SummaryPrintingTimerEx::getVerboseFileFromEnv(int id)
{
  std::string name = Environment::getStringEnv("OPENZGY_TIMERS_LOGFILE");
  if (name.empty() || name == "cerr" || name == "stderr") {
    return std::shared_ptr<std::ostream>(&std::cerr, [](std::ostream*){});
    // Never closed. The smart pointer has a dummy destructor.
  }
  else if (name == "cout" || name == "stdout" || name == "con:") {
    return std::shared_ptr<std::ostream>(&std::cout, [](std::ostream*){});
  }
  else {
    std::size_t pos = name.find("{}");
    if (pos != std::string::npos)
      name = name.substr(0, pos) + std::to_string(id) + name.substr(pos+2);
    return std::make_shared<std::ofstream>(name, std::ofstream::app);
  }
}

} // namespace
