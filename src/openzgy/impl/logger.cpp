// Copyright 2017-2020, Schlumberger
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

#include "logger.h"
#include "environment.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <atomic>
#include <memory>
#include <chrono>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Static class so it doesn't actually need to be attached as a pimpl.
 * The methods and data is used by standardCallback().
 */
class LoggerBaseImpl
{
public:
  static std::atomic<int> next_id_;
  static std::mutex mutex_;
public:
  static std::shared_ptr<std::ostream> getVerboseFileFromEnv(int id);
  static int getVerboseDetailsFromEnv();
  static double timestamp();
};

std::atomic<int> LoggerBaseImpl::next_id_{1};
std::mutex LoggerBaseImpl::mutex_{};

int
LoggerBase::getVerboseFromEnv(const char *envname)
{
  return Environment::getNumericEnv(envname, 0);
}

std::shared_ptr<std::ostream>
LoggerBaseImpl::getVerboseFileFromEnv(int id)
{
  std::string name = Environment::getStringEnv("OPENZGY_VERBOSE_LOGFILE");
  if (!name.empty()) {
    std::size_t pos = name.find("{}");
    if (pos != std::string::npos)
      name = name.substr(0, pos) + std::to_string(id) + name.substr(pos+2);
    return std::make_shared<std::ofstream>(name, std::ofstream::app);
  }
  else {
    return std::shared_ptr<std::ostream>(&std::cerr, [](std::ostream*){});
  }
}

int
LoggerBaseImpl::getVerboseDetailsFromEnv()
{
  return Environment::getNumericEnv("OPENZGY_VERBOSE_DETAILS", 0);
}

double
LoggerBaseImpl::timestamp()
{
  // Remember the good old days with ::time(nullptr) ?
  typedef std::chrono::high_resolution_clock clock;
  static constexpr double factor =
    static_cast<double>(clock::duration::period::num) /
    static_cast<double>(clock::duration::period::den);
  return clock::now().time_since_epoch().count() * factor;
}

/**
 * Invoke the specified logger function with a string argument.
 * This isn't much different from invoking the callback directly.
 * But it makes debugging slightly simpler to have an easy place
 * to set a breakpoint. It also adds more symmetry with respect
 * to the stringstream version, which does add value.
 */
bool
LoggerBase::logger(const LoggerFn& callback, int priority, const std::string& str)
{
  return callback(priority, str);
}

/**
 * Invoke the specified logger function with a stringstream argument.
 *
 * Example usage:
 *
 *    if (logger(priority))
 *      logger(priority, std::stringstream() << "Hello" << ", world.");
 */
bool
LoggerBase::logger(const LoggerFn& callback, int priority, const std::ios& ss)
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return callback(priority, sstream ? sstream->str() : std::string());
}

/**
 * To be used instead of an empty functor, so the code knows it is
 * always safe to call the logger. All output is discarded.
 */
LoggerBase::LoggerFn
LoggerBase::emptyCallback()
{
  static const LoggerFn null = [](int, const std::string&) { return false; };
  return null;
}

/**
 * Convert an old style integer level to a logger that outputs to
 * std::cerr all messages of the specified priority level and below.
 *
 * Level 42 is treated specially. It causes all logging code to be
 * executed but nothing will actually be output.
 *
 * Level < 0 will never show anything. Not even messages with pri < 0.
 *
 * Suggested use of priority levels:
 * - -1 => Serious errors, always show unless messages disbled completely.
 * -  0 => Shown during normal execution. Also in production mode.
 * - 1+ => For debugging and testing.
 *
 * By default the following environment variables are recognized:
 * - OPENZGY_VERBOSE (choose how much is dumped)
 * - OPENZGY_VERBOSE_LOGFILE (redirect output to a file)
 * - OPENZGY_VERBOSE_DETAILS (bitmask to turn on optional features)
 *
 * It is the caller of standardCallback that looks up OPENZGY_VERBOSE
 * and it can choose to use some other mechanism to enable or disable.
 * The other two are currently local to this file and can only be set
 * by these environment variables.
 *
 * Thread safety: May be used concurrently from different threads.
 * Uses a global lock. The lock might be skipped is std::ios is thread
 * safe, but that would technically trigger undefined behavior.
 */
LoggerBase::LoggerFn
LoggerBase::standardCallback(int currentlevel, const std::string& prefix_in, const std::string& suffix_in)
{
  static const LoggerFn null = [](int, const std::string&) { return false; };
  static const LoggerFn test = [](int, const std::string&) { return true; };
  if (currentlevel < 0) {
    return null;
  }
  else if (currentlevel == 42) {
    return test;
  }
  else {
    std::string prefix(prefix_in);
    std::string suffix(suffix_in + "\n");
    int id = LoggerBaseImpl::next_id_++;
    int details = LoggerBaseImpl::getVerboseDetailsFromEnv();
    std::shared_ptr<std::ostream> os =
      LoggerBaseImpl::getVerboseFileFromEnv(id);
    if (details & 2)
      prefix = prefix + "<" + std::to_string(id) + "> ";
    return [currentlevel, prefix, suffix, os, details](int pri, const std::string& msg) {
      if (!msg.empty() && pri <= currentlevel) {
        std::istringstream f(msg);
        std::string s;
        std::lock_guard<std::mutex> lk(LoggerBaseImpl::mutex_);
        if (os->good()) {
          while (std::getline(f, s, '\n')) {
            if (!s.empty()) {
              if (details & 1) {
                std::stringstream ss;
                ss << std::setprecision(3) << std::fixed
                   << LoggerBaseImpl::timestamp() << ": ";
                ss << prefix << s << suffix;
                *os << ss.str();
              }
              else {
                *os << (prefix + s + suffix);
              }
            }
          }
          // This will slow down logging to file but is safer if there
          // is a crash or if logging is done from a static destructor.
          *os << std::flush;
        }
      }
      return pri <= currentlevel;
    };
  }
}

} // namespace
