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

#pragma once

#include <string>
#include <vector>
#include <memory>

#include "../exception.h"

namespace OpenZGY {
  class IZgyTools;
  class IZgyReader;
  class IZgyWriter;
  class IZgyUtils;
  class SampleStatistics;
  class SampleHistogram;
  class FileStatistics;
  class SeismicStoreIOContext;
  class ZgyWriterArgs;
}

namespace InternalZGY { namespace CAPI {
#if 0
}}
#endif

/**
 * Base class of all handles. Instanciated as just the base class it
 * represents a successful method invocation.
 */
class ZgyHandleBase
{
public:
  enum HandleType {
    ZGY_FREE_HANDLE       = 0xDEADBEEF,
    // Ephemeral
    ZGY_SUCCESS_HANDLE     = 0x1A386235,
    ZGY_ERROR_HANDLE      = 0x7C6FFD03,
    ZGY_CLEANUP_HANDLE    = 0x59B54AE0,
    ZGY_STRING_HANDLE     = 0x024A6C25,
    // Long lived
    ZGY_READER_HANDLE     = 0x37E38FAC,
    ZGY_WRITER_HANDLE     = 0xD8B03197,
    ZGY_UTILS_HANDLE      = 0xE260EE12,
    ZGY_STATISTICS_HANDLE = 0x74B464BD,
    ZGY_HISTOGRAM_HANDLE  = 0x57F156ED,
    ZGY_FILESTATS_HANDLE  = 0xE93E9590,
    ZGY_IOCONTEXT_HANDLE  = 0x935F8DDD,
    ZGY_WRITERARGS_HANDLE = 0xCD100F46,
  };

protected:
  explicit ZgyHandleBase(HandleType magic);

public:
  virtual ~ZgyHandleBase();

private:
  static bool magicIsKnown(HandleType magic);
  static std::string magicToString(HandleType magic);

protected:
  static void checkHandle(const void* handle, HandleType expected = (HandleType)0) noexcept;
  static void checkMagic(const void* handle, HandleType expected);
  static bool checkMagicNoThrow(const void* handle, HandleType expected) noexcept;

public:
  static bool resultIsSuccess(const void *header) noexcept;
  static void freeHandle(void* handle) noexcept;

private:
  HandleType magic_;
};

/**
 * Class representing a successful result.
 * Might not be needed, perhaps just use nullptr for that?
 */
class ZgySuccessHandle : public ZgyHandleBase
{
public:
  explicit ZgySuccessHandle()
    : ZgyHandleBase(ZGY_SUCCESS_HANDLE)
  {
  }
};

/**
 * Class representing an exception from a method invocation.
 */
class ZgyErrorHandle : public ZgyHandleBase
{
  std::string type_;
  std::string message_;
public:
  explicit ZgyErrorHandle(const char *type, const char *message)
    : ZgyHandleBase(ZGY_ERROR_HANDLE)
    , type_(type)
    , message_(message)
  {
  }
  // To be exported as extern "C".
  static const char *resultGetExceptionType(const void* handle) noexcept
  {
    if (checkMagicNoThrow(handle, ZGY_ERROR_HANDLE))
      return static_cast<const ZgyErrorHandle*>(handle)->type_.c_str();
    else
      return "";
  }
  // To be exported as extern "C".
  static const char *resultGetExceptionMessage(const void* handle) noexcept
  {
    if (checkMagicNoThrow(handle, ZGY_ERROR_HANDLE))
      return static_cast<const ZgyErrorHandle*>(handle)->message_.c_str();
    else
      return "";
  }
};

/**
 * Use like ZgySuccessHandle but manages some private resources,
 */
class ZgyCleanupHandle : public ZgyHandleBase
{
  std::vector<std::shared_ptr<void>> cleanup_;
public:
  explicit ZgyCleanupHandle()
    : ZgyHandleBase(ZGY_CLEANUP_HANDLE)
  {
  }
  void cleanup(std::shared_ptr<void> ptr)
  {
    cleanup_.push_back(ptr);
  }
};

/**
 * Class representing the result from a function returning a string.
 */
class ZgyStringHandle : public ZgyHandleBase
{
  std::string result_;
public:
  explicit ZgyStringHandle(const char *result)
    : ZgyHandleBase(ZGY_STRING_HANDLE)
    , result_(result)
  {
  }
  explicit ZgyStringHandle(const std::string& result)
    : ZgyHandleBase(ZGY_STRING_HANDLE)
    , result_(result)
  {
  }
  // To be exported as extern "C".
  static const char *resultGetString(const void* handle) noexcept
  {
    if (checkMagicNoThrow(handle, ZGY_STRING_HANDLE))
      return static_cast<const ZgyStringHandle*>(handle)->result_.c_str();
    else
      return "";
  }
};

/**
 * Generic handle that wraps a single smart pointer e.g. to ZgyReader.
 * The get() method returns the unmanaged pointer after doing
 * consistency checks and casts on the handle. Throws on error.
 * The get() method is internal to the wrapper library. In contrast
 * to e.g. ZgyStringHandle where the contents can be retrieved
 * by the client code.
 *
 * Release will remove and eventually free the unmanaged ZgyXXX.
 * Only the few bytes in the ZgyHandleT that the managed SafeHandle
 * knows about will remain, to eventually get freed by oz_freeHandle().
 * The smart pointer is returned so the caller can invoke close() etc.
 */
template<typename T, ZgyHandleBase::HandleType MAGIC>
class ZgyHandleT : public ZgyHandleBase
{
  std::shared_ptr<T> ptr_;
public:
  typedef ZgyHandleT<T, MAGIC> self_type;
  explicit ZgyHandleT(const std::shared_ptr<T>& ptr)
    : ZgyHandleBase(MAGIC), ptr_(ptr) {}

  static self_type* cast(void* handle, bool throw_if_empty) {
    checkMagic(handle, MAGIC);
    auto self = static_cast<ZgyHandleT<T, MAGIC>*>(handle);
    if (throw_if_empty && !self->ptr_)
      throw OpenZGY::Errors::ZgyInternalError("ZGY C API: Handle was disposed.");
    return self;
  }
  static T* get(void* handle) {
    return cast(handle, true)->ptr_.get();
  }
  static std::shared_ptr<T> getsmart(void* handle, bool throw_if_empty) {
    return cast(handle, throw_if_empty)->ptr_;
  }
  static std::shared_ptr<T> release(void* handle)
  {
    self_type* self = cast(handle, false);
    auto victim = self->ptr_;
    self->ptr_.reset();
    return victim;
  }
};

typedef ZgyHandleT<OpenZGY::IZgyReader,       ZgyHandleBase::ZGY_READER_HANDLE> ZgyReaderHandle;
typedef ZgyHandleT<OpenZGY::IZgyWriter,       ZgyHandleBase::ZGY_WRITER_HANDLE> ZgyWriterHandle;
typedef ZgyHandleT<OpenZGY::IZgyUtils,        ZgyHandleBase::ZGY_UTILS_HANDLE>  ZgyUtilsHandle;
typedef ZgyHandleT<OpenZGY::SampleStatistics, ZgyHandleBase::ZGY_STATISTICS_HANDLE> ZgyStatisticsHandle;
typedef ZgyHandleT<OpenZGY::SampleHistogram,  ZgyHandleBase::ZGY_HISTOGRAM_HANDLE> ZgyHistogramHandle;
typedef ZgyHandleT<const OpenZGY::FileStatistics,   ZgyHandleBase::ZGY_FILESTATS_HANDLE> ZgyFileStatsHandle;
typedef ZgyHandleT<OpenZGY::SeismicStoreIOContext, ZgyHandleBase::ZGY_IOCONTEXT_HANDLE> ZgyIOContextHandle;
typedef ZgyHandleT<OpenZGY::ZgyWriterArgs,    ZgyHandleBase::ZGY_WRITERARGS_HANDLE> ZgyWriterArgsHandle;

/**
 * Access a handle that is either a ZgyReader or a ZgyWriter.
 * There is no handle specifically of the META type.
 */
class ZgyMetaHandle : public ZgyHandleBase
{
  ZgyMetaHandle() = delete;
  ZgyMetaHandle(const ZgyMetaHandle&) = delete;
public:
  static OpenZGY::IZgyTools* get(void* handle);
};

}} // namespace
