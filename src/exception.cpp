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

#include "exception.h"
#include <errno.h>
#include <string.h>
#include <string>

namespace OpenZGY { namespace Errors {
#if 0
}}
#endif

ZgyError::ZgyError(const std::string& arg) : std::runtime_error(arg) {}
ZgyFormatError::ZgyFormatError(const std::string& arg) : ZgyError(arg) {}
ZgyNeedOldLibrary::ZgyNeedOldLibrary(const std::string& arg) : ZgyFormatError(arg) {}
ZgyUpdateRules::ZgyUpdateRules(const std::string& arg) : ZgyFormatError(arg) {}
ZgyCorruptedFile::ZgyCorruptedFile(const std::string& arg) : ZgyError(arg) {}
ZgyUserError::ZgyUserError(const std::string& arg) : ZgyError(arg) {}
ZgyInternalError::ZgyInternalError(const std::string& arg) : ZgyError(arg) {}
ZgyEndOfFile::ZgyEndOfFile(const std::string& arg) : ZgyError(arg) {}
ZgySegmentIsClosed::ZgySegmentIsClosed(const std::string& arg) : ZgyError(arg) {}
ZgyAborted::ZgyAborted(const std::string& arg) : ZgyError(arg) {}
ZgyMissingFeature::ZgyMissingFeature(const std::string& arg) : ZgyError(arg) {}
ZgyNotReadOnlyError::ZgyNotReadOnlyError(const std::string& arg) : ZgyError(arg) {}

namespace {
  static std::string get_error(const std::string& filename, int system_errno)
  {
    std::string errstring;
    char buf[1024]{0};
#if defined(_WIN32)
    // crt gives us thread_local buffer, no need for strerror_s.
    // So why is the compiler still complaining?
    //errstring = std::string(strerror(system_errno));
    if (0 == strerror_s(buf, sizeof(buf) - 1, system_errno))
      errstring = std::string(buf);
#elif (_POSIX_C_SOURCE >= 200112L) && ! _GNU_SOURCE
    int rc = strerror_r(system_errno, buf, sizeof(buf)-1);
    buf[sizeof(buf)-1] = '\0';
    errstring = std::string(buf);
#else
    char *e = strerror_r(system_errno, buf, sizeof(buf)-1);
    errstring = std::string(e ? e : "");
#endif
    if (errstring.empty()) // should never happen; strerror should do it.
      errstring = "Unknown errno " + std::to_string(system_errno);
    return filename + ": " + errstring;
  }
}

ZgyIoError::ZgyIoError(const std::string& filename, int system_errno)
  : ZgyError(get_error(filename, system_errno))
{
}

}} // namespace
