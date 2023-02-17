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

#pragma once

#include "../declspec.h"

#include <string>
#include <sstream>
#include <functional>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file logger.h
 * \brief Logging framework.
 *
 * See InternalZGY::LoggerBase for details.
 */

/**
 * Logging framework.
 *
 * The actual logger is a simple functor that will be invoked whenever
 * the code wants to output a debug message. The functor may or may
 * not output the message somewhere. The return value indicates
 * whether the message was or should have been output. If called with
 * an empty message don't output anything but still return the bool
 * should_be_output.
 *
 * The logger can be stored in a global singleton or be passed around
 * as function arguments or class members. OpenZGY currently only
 * does the latter. This file used to contain a class Logger to be
 * used as a singleton. Yagni, or fetch from source control if you
 * need it.
 *
 * The callback should be guaranteed to not be empty. Attempts to set
 * an empty callback should cause a dummy callback with no output to be set.
 *
 * Simple use:
 * \code
 *   callback(priority, message);
 * \endcode
 *
 * Advanced use if message may be expensive to compute:
 * \code
 *   if (callback(priority, ""))
 *     callback(priority, getMyMessage());
 * \endcode
 *
 * The functor is allowed to lie. Specifically, when the second form
 * is used there will be code that is never executed unless the user
 * turns on a huge amount of logging. It is useful to run some of the
 * tests with a logger that discards all the output but claims it was
 * or would be shown. Then run the same test with no logging, and
 * verify there was no observed difference. This strategy also does
 * wonders for the code coverage.
 *
 * LoggerBase contains static convenience functions only. It is
 * possible to avoid this class and use the logger functor directly if
 * it is desirable to limit dependencies. What LoggerBase provides
 * is a few simple loggers that write to std::cerr or never writes.
 * There is also a wrapper to simplify logging of formatted text.
 *
 * Need to be exported because it is used by the Python wrapper.
 *
 * Thread safety: Safe because it contains only static methods.
 *
 * Any LoggerFn implementation must be thread safe.
 */
class OPENZGY_API LoggerBase
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;
  static int getVerboseFromEnv(const char *envname);
  static bool logger(const LoggerFn& logger, int priority, const std::string& str = std::string());
  static bool logger(const LoggerFn& logger, int priority, const std::ios& ss);
  static LoggerFn emptyCallback();
  static LoggerFn standardCallback(int level, const std::string& prefix, const std::string& suffix);
};

} // namespace
