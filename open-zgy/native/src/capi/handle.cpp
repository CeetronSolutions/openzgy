// Copyright 2017-2022, Schlumberger
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

#include "handle.h"
#include "../api.h"
#include "../exception.h"
#include "../iocontext.h"

#include <iostream>
#include <sstream>
#include <iomanip>

namespace InternalZGY { namespace CAPI {
#if 0
}}
#endif

ZgyHandleBase::ZgyHandleBase(HandleType magic)
  : magic_(magic)
{
}

ZgyHandleBase::~ZgyHandleBase()
{
  // This is to help debugging, making it more likely to catch use of
  // a recently freed handle. The code should not test for this particular
  // value anywhere. It should be handled just like any other unknown tag.
  magic_ = ZGY_FREE_HANDLE;
}

bool
ZgyHandleBase::magicIsKnown(ZgyHandleBase::HandleType magic)
{
  switch (magic) {
  case ZGY_SUCCESS_HANDLE:
  case ZGY_ERROR_HANDLE:
  case ZGY_CLEANUP_HANDLE:
  case ZGY_STRING_HANDLE:
  case ZGY_READER_HANDLE:
  case ZGY_WRITER_HANDLE:
  case ZGY_UTILS_HANDLE:
  case ZGY_STATISTICS_HANDLE:
  case ZGY_HISTOGRAM_HANDLE:
  case ZGY_FILESTATS_HANDLE:
  case ZGY_IOCONTEXT_HANDLE:
  case ZGY_WRITERARGS_HANDLE:
    return true;
  case ZGY_FREE_HANDLE:
  default:
    return false;
  }
}

std::string
ZgyHandleBase::magicToString(ZgyHandleBase::HandleType magic)
{
  switch (magic) {
  case ZGY_SUCCESS_HANDLE:    return "SUCCESS";
  case ZGY_ERROR_HANDLE:      return "ERROR";
  case ZGY_CLEANUP_HANDLE:    return "CLEANUP";
  case ZGY_STRING_HANDLE:     return "STRING";
  case ZGY_READER_HANDLE:     return "READER";
  case ZGY_WRITER_HANDLE:     return "WRITER";
  case ZGY_UTILS_HANDLE:      return "UTILS";
  case ZGY_STATISTICS_HANDLE: return "STATISTICS";
  case ZGY_HISTOGRAM_HANDLE:  return "HISTOGRAM";
  case ZGY_FILESTATS_HANDLE:  return "FILESTATS";
  case ZGY_IOCONTEXT_HANDLE:  return "IOCONTEXT";
  case ZGY_WRITERARGS_HANDLE: return "WRITERARGS";
  default:
    {
      std::stringstream ss;
      ss << "0x" << std::hex << std::setfill('0') << std::setw(8) << magic;
      return ss.str();
    }
  }
}

/**
 * Report a garbage handle by forcing the application to terminate.
 * A message is printed first. If the library expected a particular
 * kind of handle then this can be passed in, to be used as part
 * of the error message.
 *
 * At this level a null handle is treated as a fatal error.
 */
void
ZgyHandleBase::checkHandle(const void *handle, HandleType expected) noexcept
{
  const HandleType magic =
    handle ? static_cast<const ZgyHandleBase*>(handle)->magic_ : (HandleType)0;
  if (!handle || !magicIsKnown(magic)) {
    std::string message =
      ("ZGY C API: Fatal error. " +
       std::string(handle ? "Garbage" : "Null") +
       std::string(" handle. Expected ") +
       (expected==(HandleType)0 ? std::string("any") : magicToString(expected))+
       std::string(".\n"));
    std::cerr << message << std::flush;
    abort();
  }
}

/**
 * Output a message to the console and crash if the handle is garbage.
 * Throw an exception if the handle is good but of the wrong type.
 * Null handles at this level are treated as exceptions, not fatal.
 * A non-null handle with an empty pointer inside is not even
 * detected here. Because those pointers belong in subclasses.
 */
void
ZgyHandleBase::checkMagic(const void* handle, HandleType expected)
 {
   if (handle)
     checkHandle(handle, expected);
   const HandleType magic =
     handle ? static_cast<const ZgyHandleBase*>(handle)->magic_ : (HandleType)0;
   if (!handle || magic != expected) {
    std::stringstream ss;
    ss << "ZGY C API: Wrong handle type "
       << (handle ? magicToString(magic) : std::string("null"))
       << ", expected "
       << magicToString(expected) << ".";
    //std::cerr << "ZgyHandleBase::checkMagic: " << ss.str() << std::endl;
    throw OpenZGY::Errors::ZgyInternalError(ss.str());
  }
}

/**
 * Output a message to the console and crash if the handle is garbage.
 * Return false if the handle is of the wrong type. To be used from nothrow
 * methods such as oz_resultGetString()
 */
bool
ZgyHandleBase::checkMagicNoThrow(const void* handle, HandleType expected) noexcept
{
  if (handle) {
    checkHandle(handle, expected);
    const HandleType magic = static_cast<const ZgyHandleBase*>(handle)->magic_;
    return magic == expected;
  }
  else {
    return false;
  }
}

/**
 * Return false if the handle represents an error, else true.
 * Crash if the header is garbage.
 *
 * This method will be exposed as extern "C" oz_resultIsSuccess(ZgyHandle).
 */
bool
ZgyHandleBase::resultIsSuccess(const void *handle) noexcept
{
  if (!handle)
    return true; // nullptr may be used instead of SUCCESS.
  checkHandle(handle);
  const int magic = static_cast<const ZgyHandleBase*>(handle)->magic_;
  return magic != ZGY_ERROR_HANDLE;
}

/**
 * Destroy the handle and release all its resources.
 * This method will be exposed as extern "C" oz_freeHandle(ZgyHandle).
 */
void
ZgyHandleBase::freeHandle(void* handle) noexcept
{
  if (handle) {
    checkHandle(handle);

    // This check could probably have beem more elegant.
    // Some types MUST be disposed, not just finalized.
    // Print a warning in shutdown.
    std::string name;
    std::shared_ptr<void> victim;
    try {
      if (checkMagicNoThrow(handle, ZGY_READER_HANDLE)) {
        name = "ZgyReader";
        victim = ZgyReaderHandle::getsmart(handle, false);
      }
      else if (checkMagicNoThrow(handle, ZGY_WRITER_HANDLE)) {
        name = "ZgyWriter";
        victim = ZgyWriterHandle::getsmart(handle, false);
      }
      else if (checkMagicNoThrow(handle, ZGY_UTILS_HANDLE)) {
        name = "ZgyUtils";
        victim = ZgyUtilsHandle::getsmart(handle, false);
      }
    }
    catch(...)
    {
      // Cannot happen. getsmart() can technically throw, but not after
      // the magic number has already been checked
      std::cerr << "ERROR: Problems disposing a ZgyBase." << std::endl;
      // I choose to continue with the regular free, not the leak.
    }
    //if (!name.empty() && !victim) {
    //  std::cerr << "Debug: Free " + name + " that was already closed." << std::endl;
    //}
    if (!name.empty() && victim)
    {
      std::cerr << "ERROR: Failed to dispose " + name << "." << std::endl;
      // It is too dangerous to release the smart pointer and call
      // close() from a finalizer. Leak the handle instead. Could have
      // tried to free "safe" resources only but there is not much point.
      static_cast<ZgyHandleBase*>(handle)->magic_ = ZGY_FREE_HANDLE;
    }
    else {
      delete static_cast<const ZgyHandleBase*>(handle);
    }
  }
}

OpenZGY::IZgyTools*
ZgyMetaHandle::get(void* handle) {
  if (checkMagicNoThrow(handle, ZGY_READER_HANDLE)) {
    return ZgyReaderHandle::get(handle);
  }
  else if (checkMagicNoThrow(handle, ZGY_WRITER_HANDLE)) {
    return ZgyWriterHandle::get(handle);
  }
  else {
    checkMagic(handle, ZGY_READER_HANDLE);
    return nullptr; // actually never reached.
  }
}

void test_instanciate_handles()
{
  // Hand written
  ZgySuccessHandle success;
  ZgyErrorHandle error("error", "oops!");
  ZgyCleanupHandle cleanup;
  ZgyStringHandle string("hello!");
  // From a single template
  ZgyReaderHandle reader(nullptr);
  ZgyWriterHandle writer(nullptr);
  ZgyUtilsHandle utils(nullptr);
  ZgyStatisticsHandle stats(nullptr);
  ZgyHistogramHandle histo(nullptr);
  ZgyFileStatsHandle fstat(nullptr);
  ZgyIOContextHandle ctxt(nullptr);
  ZgyWriterArgsHandle args(nullptr);
}

}} // namespace

