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

#pragma region includes and local typedefs

#include "capi.h"
#include "handle.h"
#include "../api.h"
#include "../iocontext.h"
#include "../exception.h"
#include "../declspec.h"

#include <functional>
#include <cstdint>
#include <tuple>
#include <string.h>

using namespace InternalZGY::CAPI;
typedef void* ZgyHandle;
typedef void* CALLBACK; // placeholder
typedef int ENUM;

namespace CallbackTypes
{
  // Callbacks in the C++ API.
  typedef std::function<bool(int, const std::string&)> cxx_logger_t;
  typedef std::function<bool(std::int64_t, std::int64_t)> cxx_progress_t;
  typedef std::function<std::string()> cxx_tokencb_t;

  // class XXXForwarder in C++ (currently in this file) takes a C callback in
  // and returns a C++ std::function suitable for passing to a method in the
  // C++ API. It will be used in the implementation of the oz_xxx function
  // (also currently in this file) i.e. the external "C" linkage function.

  // Callbacks in the C api, note int instead of bool.
  typedef int(*c_logger_t)(int, const char*);
  typedef int(*c_progress_t)(std::int64_t, std::int64_t);
  typedef const char* (*c_tokencb_v1_t)();
  typedef std::size_t(*c_tokencb_v2_t)(char*, std::size_t);
  typedef const char* (*c_tokencb_v3_t)();

  // Callbacks in the C# API, marshaling to the C API.
  // These are used in the [DllExport] declarations.
  //[UnmanagedFunctionPointer(CallingConvention.Cdecl)]public delegate int InteropLoggerCallback(int level, string message);
  //[UnmanagedFunctionPointer(CallingConvention.Cdecl)]public delegate int InteropProgressCallback(long done, long total);
  //[UnmanagedFunctionPointer(CallingConvention.Cdecl)]public delegate string InteropTokenCallbackV1();
  //[UnmanagedFunctionPointer(CallingConvention.Cdecl)]public delegate long InteropTokenCallbackV2(long buffersize, string message);
  //[UnmanagedFunctionPointer(CallingConvention.Cdecl)]public delegate System.IntPtr InteropTokenCallbackV3();

  // Callbacks::forwardXXX takes an application XXXDelegate instance and will
  // return an InteropXXXCallback. The forwarder is unsed in the C# funcion
  // that invokes the oz_xxx function i.e. the external "C" linkage function.

  // Callbacks in the C# API, application visible.
  //public delegate bool ProgressDelegate(long done, long total);
  //public delegate bool LoggerDelegate(int level, string message);
  //public delegate string TokenDelegate();

  // Beware object lifetimes for the logger and the tokencb.
  // The functor returned by the forwarder will be held onto in C++ code,
  // and inside the functor is a pointer to the C callback to be invoked.
  // On the C# end the code there is responsible for keeping the C# delegate
  // alive. As a minimum store it in the ZgyMeta instance. Not the IOContext
  // which is ephemeral.
  // The lifetime of the progress callback is only the call in which it is
  // used. So that is less of a problem.
}

#pragma endregion

namespace {

#pragma region "Internal functions with C++ linkage"
#if 0
}
#endif

/**
 * Convenience function for wrapping any method.
 * An exception is converted to a return of a ZgyErrorHandle.
 * Otherwise the functor's return value becomes the end result.
 * If the functor returned null or didn't have a return type
 * then a ZgyResultHandle, i.e. a plain success, is returned.
 *
 * Might I actually treat "null handle" as success
 * and not have to allocate all those handles?
 * Yes, but if it works (and runs fast enough) don't fix it.
 *
 * Note that most of the ~130 methods return just a short lived
 * pass/fail handle. Those should use protect() instead of
 * protectReturn().
 */
static ZgyHandleBase* protectReturn(const std::function<ZgyHandleBase* ()>& fn)
{
  try {
    ZgyHandleBase* result = fn();
    return result ? result : new ZgySuccessHandle();
  }
  // The list of recognized exceptions should be kept in sync with
  // ../exception.h. Be very careful about ordering. Specialized
  // errors must be listed before the base class.
  catch (const OpenZGY::Errors::ZgyNeedOldLibrary& ex) { return new ZgyErrorHandle("ZgyNeedOldLibrary", ex.what()); }
  catch (const OpenZGY::Errors::ZgyUpdateRules& ex) { return new ZgyErrorHandle("ZgyUpdateRules", ex.what()); }
  catch (const OpenZGY::Errors::ZgyFormatError& ex) { return new ZgyErrorHandle("ZgyFormatError", ex.what()); }
  catch (const OpenZGY::Errors::ZgyCorruptedFile& ex) { return new ZgyErrorHandle("ZgyCorruptedFile", ex.what()); }
  catch (const OpenZGY::Errors::ZgyUserError& ex) { return new ZgyErrorHandle("ZgyUserError", ex.what()); }
  catch (const OpenZGY::Errors::ZgyInternalError& ex) { return new ZgyErrorHandle("ZgyInternalError", ex.what()); }
  catch (const OpenZGY::Errors::ZgyEndOfFile& ex) { return new ZgyErrorHandle("ZgyEndOfFile", ex.what()); }
  catch (const OpenZGY::Errors::ZgySegmentIsClosed& ex) { return new ZgyErrorHandle("ZgySegmentIsClosed", ex.what()); }
  catch (const OpenZGY::Errors::ZgyAborted& ex) { return new ZgyErrorHandle("ZgyAborted", ex.what()); }
  catch (const OpenZGY::Errors::ZgyMissingFeature& ex) { return new ZgyErrorHandle("ZgyMissingFeature", ex.what()); }
  catch (const OpenZGY::Errors::ZgyIoError& ex) { return new ZgyErrorHandle("ZgyIoError", ex.what()); }
  catch (const OpenZGY::Errors::ZgyNotReadOnlyError& ex) { return new ZgyErrorHandle("ZgyNotReadOnlyError", ex.what()); }
  catch (const OpenZGY::Errors::ZgyError& ex) { return new ZgyErrorHandle("ZgyError", ex.what()); }
  catch (const std::exception& ex) { return new ZgyErrorHandle("std::exception", ex.what()); }
  catch (...) { return new ZgyErrorHandle("unknown", "unknown exception"); }
}

/**
 * Convenience function for wrapping a method that only returns a
 * simple pass/fail result handle and not a longer lived handle with
 * more information. Most wrappers will use this. The functor itself
 * does not return a value, it is deemed a success if it did not throw
 * an exception.
 *
 * Using protect() instead of protectReturn() saves one line of
 * "return null" in ~100 methods.
 */
static ZgyHandleBase* protect(const std::function<void()>& fn)
{
  return protectReturn([&]() { fn(); return nullptr; });
}

/**
 * Convenience function for wrapping a method that returns int64[3]
 * a.k.a. size3i_t into a C function with 3 out parameters.
 */
static std::tuple<std::int64_t, std::int64_t, std::int64_t>
sizeToTuple(const std::array<std::int64_t, 3>& s)
{
  return std::make_tuple(s[0], s[1], s[2]);
}

/**
 * Convenience function for wrapping a method that returns float[3]
 * into a C function with 2 out parameters.
 */
static std::tuple<float, float>
rangeToTuple(const std::array<float, 2>& s)
{
  return std::make_tuple(s[0], s[1]);
}

/**
 * Convenience function for wrapping a method that returns float[3]
 * into a C function with 2 out parameters.
 */
static std::tuple<double, double>
doubleToTuple(const std::array<double, 2>& s)
{
  return std::make_tuple(s[0], s[1]);
}

/**
 * Convert a std::function callback instance to something that can be
 * passed down to a function expecting an unsafe C function pointer.
 *
 * The documentation comments for this class also apply to
 * ProgressForwarder and TokenCallbackForwarder, so those are not
 * documented themselves.
 */
class LoggerForwarder
{
private:
  int(*callback_)(int, const char*);

public:
  /**
   * The constructor's argument  is what we got from the client,
   * marshaled as an unsafe C-style callback.
   */
  explicit LoggerForwarder(int(*callback)(int, const char*))
    : callback_(callback)
  {
    //std::cerr << "@@ new LoggerForwarder " << std::hex << (std::int64_t)callback << std::dec << std::endl;
  }

  /**
   * Invoke code in the client. Possibly on a foreign thread. Knows how
   * to handle the case where the client doesn't want to be bothered.
   * If get() was used to get the functor we won't actually get here
   * with a null callback_, and the library decides what the default is.
   */
  bool operator()(int level, const std::string& message)
  {
    return callback_ ? (*callback_)(level, message.c_str()) != 0 : false;
  }

  /**
   * Get a std::function that can be passed to a method in this library.
   * A null unsafe callback becomes an empty std::function, telling the
   * function in the C++ API that it should not invoke the callback.
   */
  std::function<bool(int, const std::string&)> get()
  {
    if (callback_)
      return *this;
    return nullptr;
  }
};

/**
 * See LoggerForwarder for an explanation of how this works.
 */
class ProgressForwarder
{
private:
  int(*callback_)(std::int64_t, std::int64_t);

public:
  explicit ProgressForwarder(int(*callback)(std::int64_t, std::int64_t))
    : callback_(callback)
  {
    //std::cerr << "@@ new ProgressForwarder " << std::hex << (std::int64_t)callback << std::dec << std::endl;
  }

  bool operator()(std::int64_t done, std::int64_t total)
  {
    return callback_ ? (*callback_)(done, total) != 0 : 1;
  }

  std::function<bool(std::int64_t, std::int64_t)> get()
  {
    if (callback_)
      return *this;
    return nullptr;
  }
};

/**
 * See LoggerForwarder for an explanation of how this works.
 *
 * This version assumes a C-style callback with no parameters and the
 * result being returned as a char*. CAVEAT: The string might be leaked
 * because the C++ code invoking the callback doesn't know how to free it
 * and the code that implements it doesn't know when it is safe to do so.
 * Or worse, the returned string might be freed too early.
 */
class TokenCallback1Forwarder
{
private:
  const char* (*callback_)();

public:
  TokenCallback1Forwarder(const char* (*cb)())
    : callback_(cb)
  {
    //std::cerr << "@@ new TokenCallback1Forwarder " << std::hex << (std::int64_t)callback_ << std::dec << std::endl;
  }

  std::string operator()()
  {
    const char* ret = callback_ ? (*callback_)() : "";
    return std::string(ret ? ret : "");
  }

  std::function<std::string()> get()
  {
    if (callback_)
      return *this;
    return nullptr;
  }
};

/**
 * See LoggerForwarder for an explanation of how this works.
 *
 * This version is nearly identical to version 1, except the
 * ownership of the string's memory gets transferred to us.
 * The memory was allocated by malloc so we know how to free it,
 */
class TokenCallback3Forwarder
{
private:
  const char* (*callback_)();

public:
  TokenCallback3Forwarder(const char* (*cb)())
    : callback_(cb)
  {
    //std::cerr << "@@ new TokenCallback3Forwarder " << std::hex << (std::int64_t)callback_ << std::dec << std::endl;
  }

  std::string operator()()
  {
    std::string result;
    if (callback_) {
      const char* ret = (*callback_)();
      if (ret) {
        result = std::string(ret);
        //std::cerr << "@@ Free " << std::hex << (std::int64_t)ret << std::dec << std::endl;
        // Cast is ok because we allocated the buffer.
        ::free(const_cast<char*>(ret));
      }
    }
    return result;
  }

  std::function<std::string()> get()
  {
    if (callback_)
      return *this;
    return nullptr;
  }
};

/**
 * See LoggerForwarder for an explanation of how this works.
 *
 * This version expects a C-style callback that will copy out the
 * result to a buffer provided by the code that implements the callback.
 * CAVEAT: In some languages it might be difficult or impossible
 * to handle return by reference.
 */
class TokenCallback2Forwarder
{
private:
  std::size_t(*callback_)(char*, std::size_t);

public:
  TokenCallback2Forwarder(std::size_t(*cb)(char*, std::size_t))
    : callback_(cb)
  {
  }

  std::string operator()()
  {
    if (callback_) {
      std::vector<char> buffer(10); // For testing, make it too small.
      std::size_t len = (*callback_)(buffer.data(), buffer.size());
      while (len >= buffer.size()) {
        // Even if len is exactly the string length there would
        // not have been space for the terminating null.
        // Note that there is technically a race condition here,
        // because the result migh change to a larger string between
        // calls. But the code is protected against buffer overflow.
        // And assuming that the result doesn't keep growing
        // indefinitely the loop will eventually get it right..
        buffer.resize(len + 1);
        len = (*callback_)(buffer.data(), buffer.size());
      }
      return std::string(buffer.data());
    }
    else {
      return std::string();
    }
  }

  std::function<std::string()> get()
  {
    if (callback_)
      return *this;
    return nullptr;
  }
};

#pragma endregion

} // end anonymous namespace

#pragma region Tests and experiments

class CallbackTest
{
public:
  // Callbacks in the C++ API
  typedef std::function<bool(int, const std::string&)> Logger_t;
  typedef std::function<bool(std::int64_t, std::int64_t)> progress_t;
  typedef std::function<std::string()> tokencb_t;

  // Callbacks in the C api, note int instead of bool.
  typedef int(*c_logger)(int, const char*);
  typedef int(*c_progress)(std::int64_t, std::int64_t);
  typedef const char* (*c_tokencb_v1)();
  typedef std::size_t(*c_tokencb_v2)(char*, std::size_t);
  typedef const char* (*c_tokencb_v3)();

  // Callbacks in the C# API, matching the C api.
  //public delegate int InteropLoggerCallback(int level, string message);
  //public delegate int InteropProgressCallback(long done, long total);
  //public delegate string InteropTokenCallbackV1();
  //public delegate long InteropTokenCallbackV2(long buffersize, string message);
  //public delegate System.IntPtr InteropTokenCallbackV3();

  // Callbacks in the C# API, application visible.
  //public delegate bool ProgressDelegate(long done, long total);
  //public delegate bool LoggerDelegate(int level, string message);
  //public delegate string TokenDelegate();

  /**
   * Test method that takes all three callback functions as input.
   * * Then call all three callbacks and print the result.
   */
  static void TestExampleV1(Logger_t logger, progress_t progress, tokencb_t tokencb)
  {
    if (logger) {
      /*bool ok_logger = */logger(2, "Testing from example_1");
      //std::cerr << "CAPI: Logger delegate returned " << ok_logger << std::endl;
    }
    else
    {
      std::cerr << "No logger delegate provided." << std::endl;
    }
    if (progress)
    {
      /*bool ok_progress = */progress(5, 42);
      //std::cerr << "CAPI: progress delegate returned " << ok_progress << std::endl;
    }
    else
    {
      std::cerr << "No progress delegate provided." << std::endl;
    }
    if (tokencb)
    {
      std::string my_token = tokencb();
      //std::cerr << "CAPI: token delegate returned \"" << my_token << "\"" << std::endl;
    }
    else
    {
      std::cerr << "No progress delegate provided." << std::endl;
    }
  }
};

#pragma endregion

extern "C" {

#pragma region Methods not returning a normal handle

/////////////////////////////////////////////////////////////////////////////
/// Functions that do not follow the normal pattern of returning a handle ///
/// and do not use a try/catch block. Be careful to not call anything     ///
/// that might throw. Non fatal errors (wrong handle type) will return    ///
/// the default value e.g. 0, false, or "".                               ///
/////////////////////////////////////////////////////////////////////////////

#if 0
}
#endif

#define LIB_MAJOR 0
#define LIB_MINOR 1
#define LIB_PATCH 1

OPENZGY_API ZgyHandle oz_checkLibraryVersion(unsigned int version)
{
  unsigned int major = (version >> 16) & 0xFFFF;
  unsigned int minor = (version >> 8) & 0xFF;
  //unsigned int patch = version & 0xFF;
  // major, minor, patch is from the C# wrapper
  // LIB_ is compiled into this code, the C wrapper.
  // If the C# wrapper is seriously older (different major) it will fail.
  // If the C# wrapper is even slightly newer this will fail.
  if (major < LIB_MAJOR)
    return new ZgyErrorHandle("LibraryVersionMismatch", "The C wrapper binary is too new");
  else if (major > LIB_MAJOR || minor > LIB_MINOR)
    return new ZgyErrorHandle("LibraryVersionMismatch", "The C wrapper binary is too new");
  else
    return new ZgySuccessHandle();
}

OPENZGY_API int oz_resultIsSuccess(ZgyHandle handle)
{
  // this method is noexcept
  return ZgyHandleBase::resultIsSuccess(handle) ? 1 : 0;
}

OPENZGY_API const char* oz_resultGetExceptionType(ZgyHandle handle)
{
  // this method is noexcept
  return ZgyErrorHandle::resultGetExceptionType(handle);
}

OPENZGY_API const char* oz_resultGetExceptionMessage(ZgyHandle handle)
{
  // this method is noexcept
  return ZgyErrorHandle::resultGetExceptionMessage(handle);
}

OPENZGY_API const char* oz_resultGetString(ZgyHandle handle)
{
  // this method is noexcept
  return ZgyStringHandle::resultGetString(handle);
}

OPENZGY_API void oz_freeHandle(ZgyHandle handle)
{
  // this method is noexcept
  ZgyHandleBase::freeHandle(handle);
}

OPENZGY_API void* oz_malloc(std::int64_t size)
{
  void* data = ::calloc(size, 1);
  //std::cerr << "@@ Allocate " << size << " bytes  at " << std::hex << (std::int64_t)data << std::dec << std::endl;
  return data;
}

#pragma endregion

#pragma region Create instances

/////////////////////////////////////////////////////////////////////////////
/// Functions that create new instances from scratch, normally mapped to  ///
/// a static method in the C++ API and returning a nontrivial handle.     ///
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_reader_open(const char *filename, ZgyHandle iocontext)
{
  return protectReturn([&]() {
    return new ZgyReaderHandle
      (OpenZGY::IZgyReader::open
       (std::string(filename ? filename : "(null)"),
        iocontext ? ZgyIOContextHandle::get(iocontext) : nullptr));});

  /* why make it readable?
  auto fn = [&]() -> ZgyReaderHandle* {
     auto ctxt = iocontext ? ZgyIOContextHandle::get(iocontext) : nullptr;
     auto name = std::string(filename ? filename : "(null)");
     auto reader = OpenZGY::IZgyReader::open(name, ctxt);
     return new ZgyReaderHandle(reader);
  };
  return protectReturn(fn);
  */
}

OPENZGY_API ZgyHandle oz_writer_open(ZgyHandle args)
{
  return protectReturn([&]() {
    return new ZgyWriterHandle
      (OpenZGY::IZgyWriter::open
       (*ZgyWriterArgsHandle::get(args)));});
}

OPENZGY_API ZgyHandle oz_writer_reopen(ZgyHandle args)
{
  return protectReturn([&]() {
    return new ZgyWriterHandle
      (OpenZGY::IZgyWriter::reopen
       (*ZgyWriterArgsHandle::get(args)));});
}

OPENZGY_API ZgyHandle oz_utils_utils(const char *filename, ZgyHandle iocontext)
{
  return protectReturn([&]() {
    return new ZgyUtilsHandle
      (OpenZGY::IZgyUtils::utils
       (std::string(filename ? filename : "(null)"),
        iocontext ? ZgyIOContextHandle::get(iocontext) : nullptr));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_create()
{
  // Actually this one cannot throw, but why not protect it anyway?
  return protectReturn([&]() {
    return new ZgyIOContextHandle
      (std::make_shared<OpenZGY::SeismicStoreIOContext>());});
}

OPENZGY_API ZgyHandle oz_writerargs_create()
{
  // Actually this one cannot throw, but why not protect it anyway?
  return protectReturn([&]() {
    return new ZgyWriterArgsHandle
      (std::make_shared<OpenZGY::ZgyWriterArgs>());});
}

#pragma endregion

#pragma region Obtain instances

/////////////////////////////////////////////////////////////////////////////
/// Functions that obtain new instances from some member function.        ///
/// Returning a nontrivial handle.                                        ///
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_meta_statistics(ZgyHandle handle)
{
  // The accessor returns an actual instance, need to make a smart pointer.
  return protectReturn([&]() {
    return new ZgyStatisticsHandle
      (std::make_shared<OpenZGY::SampleStatistics>
       (ZgyMetaHandle::get(handle)->statistics()));});
}

OPENZGY_API ZgyHandle oz_meta_histogram(ZgyHandle handle)
{
  // The accessor returns an actual instance, need to make a smart pointer.
  return protectReturn([&]() {
    return new ZgyHistogramHandle
      (std::make_shared<OpenZGY::SampleHistogram>
       (ZgyMetaHandle::get(handle)->histogram()));});
}

OPENZGY_API ZgyHandle oz_meta_filestats(ZgyHandle handle)
{
  return protectReturn([&]() {
    return new ZgyFileStatsHandle
      (ZgyMetaHandle::get(handle)->filestats());});
}

OPENZGY_API ZgyHandle oz_ssiocontext_clone(ZgyHandle handle)
{
  // This is messy because the C api doesn't match the polymorphic IOContext.
  return protectReturn([&]() {
    if (!handle)
      throw std::runtime_error("Unable to clone a null IOContext.");
    std::shared_ptr<OpenZGY::IOContext> cloned = ZgyIOContextHandle::get(handle)->clone();
    if (!cloned)
      throw std::runtime_error("Unable to clone the IOContext.");
    auto cloned_ss = std::dynamic_pointer_cast<OpenZGY::SeismicStoreIOContext>(cloned);
    if (!cloned_ss)
      throw std::runtime_error("IOContext type not recognized by the C API.");
    return new ZgyIOContextHandle(cloned_ss);
    });
}

#pragma endregion

#pragma region Important functions

/////////////////////////////////////////////////////////////////////////////
/// Important functions that are not plain getters or setters,            ///
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_reader_read_float(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, float* data, int lod)
{
  return protect([&]() {ZgyReaderHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data, lod);});
}

OPENZGY_API ZgyHandle oz_reader_read_int16(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, std::int16_t* data, int lod)
{
  return protect([&]() {ZgyReaderHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data, lod);});
}

OPENZGY_API ZgyHandle oz_reader_read_int8(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, std::int8_t* data, int lod)
{
  return protect([&]() {ZgyReaderHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data, lod);});
}

OPENZGY_API ZgyHandle oz_reader_readconst(ZgyHandle handle, int* is_const_intflag, double* value, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, int lod, int as_float_intflag)
{
  return protect([&]() {
                   *is_const_intflag = 0;
                   *value = 0;
                   std::pair<bool, double> ret = ZgyReaderHandle::get(handle)->readconst(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, lod, as_float_intflag != 0);
                   *is_const_intflag = ret.first ? 1 : 0;
                   *value = ret.second;
                 });
}

OPENZGY_API ZgyHandle oz_writer_read_float(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, float* data)
{
  // Without the lod argument, and expects a different handle type.
  // Applies to all 4 "read" methods available in the writer.
  return protect([&]() {ZgyWriterHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_read_int16(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, std::int16_t* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_read_int8(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, std::int8_t* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->read(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_readconst(ZgyHandle handle, int* is_const_intflag, double* value, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, int as_float_intflag)
{
  return protect([&]() {
                   *is_const_intflag = 0;
                   *value = 0;
                   std::pair<bool, double> ret = ZgyWriterHandle::get(handle)->readconst(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, as_float_intflag != 0);
                   *is_const_intflag = ret.first ? 1 : 0;
                   *value = ret.second;
                 });
}

OPENZGY_API ZgyHandle oz_writer_write_float(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const float* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->write(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_write_int16(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const std::int16_t* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->write(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_write_int8(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const std::int8_t* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->write(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_writeconst_float(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const float* data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->writeconst(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_writeconst_int16(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const std::int16_t *data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->writeconst(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_writer_writeconst_int8(ZgyHandle handle, std::int64_t i0, std::int64_t j0, std::int64_t k0, std::int64_t ni, std::int64_t nj, std::int64_t nk, const std::int8_t *data)
{
  return protect([&]() {ZgyWriterHandle::get(handle)->writeconst(std::array<std::int64_t,3>{i0,j0,k0}, std::array<std::int64_t,3>{ni,nj,nk}, data);});
}

OPENZGY_API ZgyHandle oz_reader_close(ZgyHandle handle)
{
  return protect([&]() {
    auto ptr = ZgyReaderHandle::release(handle);
    if (ptr)
      ptr->close();
  });
}

OPENZGY_API ZgyHandle oz_writer_finalize(ZgyHandle handle, ENUM* decimation, int num_decimation, CallbackTypes::c_progress_t progress, ENUM action,int force_intflag)
{
  return protect([&]() {
    std::vector<OpenZGY::DecimationType> cxx_decimation;
    if (decimation != nullptr)
      for (int ii = 0; ii < num_decimation; ++ii)
        cxx_decimation.push_back((OpenZGY::DecimationType)decimation[ii]);
    ProgressForwarder cxx_progress(progress);
    auto cxx_action = (OpenZGY::FinalizeAction)action;
    ZgyWriterHandle::get(handle)->finalize(cxx_decimation, cxx_progress, cxx_action, force_intflag != 0);
  });
}

OPENZGY_API ZgyHandle oz_writer_close(ZgyHandle handle)
{
  return protect([&]() {
    auto ptr = ZgyWriterHandle::release(handle);
    if (ptr)
      ptr->close();
  });
}

OPENZGY_API ZgyHandle oz_writer_close_incomplete(ZgyHandle handle)
{
  return protect([&]() {
    auto ptr = ZgyWriterHandle::release(handle);
    if (ptr)
      ptr->close_incomplete();
  });
}

OPENZGY_API ZgyHandle oz_utils_deletefile(ZgyHandle handle, const char *name, int missing_ok_intflag)
{
  return protect([&]() {ZgyUtilsHandle::get(handle)->deletefile(name, missing_ok_intflag != 0); });
}

OPENZGY_API ZgyHandle oz_utils_alturl(ZgyHandle handle, const char* name)
{
  // One of the very few methods that return a string.
  return protectReturn([&]() {
    return new ZgyStringHandle
    (ZgyUtilsHandle::get(handle)->alturl(name)); });
}

OPENZGY_API ZgyHandle oz_utils_idtoken(ZgyHandle handle)
{
  // One of the very few methods that return a string.
  return protectReturn([&]() {
    return new ZgyStringHandle
    (ZgyUtilsHandle::get(handle)->idtoken()); });
}

OPENZGY_API ZgyHandle oz_utils_close(ZgyHandle handle)
{
  return protect([&]() {
    auto ptr = ZgyUtilsHandle::release(handle);
    // There is no close() for this type.
    });
}

#pragma endregion

#pragma region Functions with callbacks

/////////////////////////////////////////////////////////////////////////////
/// Functions that involve callbacks. Ouch.                               ///
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_test_example_v1(
  int(*logger)(int, const char*),
  int(*progress)(std::int64_t, std::int64_t),
  const char* (*tokencb)())
{
  return protect([&]() {
    // Note that protectXXX in C#, responsible for invoking the actual
    // application callback, will also need a try/catch/swallow to make
    // sure exceptions aren't propagated C# to C.
    LoggerForwarder loggerf(logger);
    ProgressForwarder progressf(progress);
    TokenCallback1Forwarder tokencb1f(tokencb);
    CallbackTest::TestExampleV1(loggerf.get(), progressf.get(), tokencb1f.get());
    });
}

typedef bool(*c_logger)(int, const char*);
typedef bool(*c_progress)(std::int64_t, std::int64_t);
typedef const char* (*c_tokencb_v1)();
typedef std::size_t(*c_tokencb_v2)(char*, std::size_t);
typedef const char* (*c_tokencb_v3)();

OPENZGY_API ZgyHandle oz_test_example_v2(
     int(*logger)(int, char*),
     int(*progress)(std::int64_t, std::int64_t),
     std::size_t(*tokencb)(char*, std::size_t))
{
  // Left as an exercise for the reader ;-)
  return new ZgySuccessHandle();
}

OPENZGY_API ZgyHandle oz_test_example_v3(
  int(*logger)(int, const char*),
  int(*progress)(std::int64_t, std::int64_t),
  const char*(*tokencb)())
{
  return protect([&]() {
    LoggerForwarder loggerf(logger);
    ProgressForwarder progressf(progress);
    TokenCallback3Forwarder tokencb3f(tokencb);
    // Using TokenCallback3Forwarder ensures the result is freed.
    // On the C# side, forwardTokenV3 ensures it is malloced.
    CallbackTest::TestExampleV1(loggerf.get(), progressf.get(), tokencb3f.get());
    });
}

#pragma endregion

#pragma region ZgyMeta getters

/////////////////////////////////////////////////////////////////////////////
/// Vanilla wrappers for member function getters and setters.             ///
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
///   ZgyMeta   /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/* Simplest possible wrap, with a single out parameters.
OPENZGY_API ZgyHandle oz_prefix_EXAMPLE(ZgyHandle handle, TYPE *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->EXAMPLE();});
}
*/

OPENZGY_API ZgyHandle oz_meta_size(ZgyHandle handle, std::int64_t* ni_return, std::int64_t* nj_return, std::int64_t *nk_return)
{
  return protect([&]() {std::tie(*ni_return, *nj_return, *nk_return) = sizeToTuple(ZgyMetaHandle::get(handle)->size());});
}

OPENZGY_API ZgyHandle oz_meta_datatype(ZgyHandle handle, ENUM *value)
{
  return protect([&]() {*value = (ENUM)ZgyMetaHandle::get(handle)->datatype();});
}

OPENZGY_API ZgyHandle oz_meta_datarange(ZgyHandle handle, float* lo_return, float* hi_return)
{
  return protect([&]() {std::tie(*lo_return, *hi_return) = rangeToTuple(ZgyMetaHandle::get(handle)->datarange());});
}

OPENZGY_API ZgyHandle oz_meta_raw_datarange(ZgyHandle handle, float* lo_return, float* hi_return)
{
  return protect([&]() {std::tie(*lo_return, *hi_return) = rangeToTuple(ZgyMetaHandle::get(handle)->raw_datarange());});
}

OPENZGY_API ZgyHandle oz_meta_zunitdim(ZgyHandle handle, ENUM *value)
{
  return protect([&]() {*value = (ENUM)ZgyMetaHandle::get(handle)->zunitdim();});
}

OPENZGY_API ZgyHandle oz_meta_hunitdim(ZgyHandle handle, ENUM *value)
{
  return protect([&]() {*value = (ENUM)ZgyMetaHandle::get(handle)->hunitdim();});
}

OPENZGY_API ZgyHandle oz_meta_zunitname(ZgyHandle handle)
{
  // One of the very few methods that return a string.
  return protectReturn([&]() {
    return new ZgyStringHandle
      (ZgyMetaHandle::get(handle)->zunitname());});
}

OPENZGY_API ZgyHandle oz_meta_hunitname(ZgyHandle handle)
{
  // One of the very few methods that return a string.
  return protectReturn([&]() {
    return new ZgyStringHandle
      (ZgyMetaHandle::get(handle)->hunitname());});
}

OPENZGY_API ZgyHandle oz_meta_zunitfactor(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->zunitfactor();});
}

OPENZGY_API ZgyHandle oz_meta_hunitfactor(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->hunitfactor();});
}

OPENZGY_API ZgyHandle oz_meta_zstart(ZgyHandle handle, float *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->zstart();});
}

OPENZGY_API ZgyHandle oz_meta_zinc(ZgyHandle handle, float *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->zinc();});
}

OPENZGY_API ZgyHandle oz_meta_annotstart(ZgyHandle handle, float* il_return, float* xl_return)
{
  return protect([&]() {std::tie(*il_return, *xl_return) = rangeToTuple(ZgyMetaHandle::get(handle)->annotstart());});
}

OPENZGY_API ZgyHandle oz_meta_annotinc(ZgyHandle handle, float* il_return, float* xl_return)
{
  return protect([&]() {std::tie(*il_return, *xl_return) = rangeToTuple(ZgyMetaHandle::get(handle)->annotinc());});
}

OPENZGY_API ZgyHandle oz_meta_corners(ZgyHandle handle, double *X00, double *Y00, double* XN0, double* YN0, double* X0M, double* Y0M, double* XNM, double* YNM)
{
  return protect([&]() {
    OpenZGY::IZgyMeta::corners_t ret = ZgyMetaHandle::get(handle)->corners();
    *X00 = ret[0][0];
    *Y00 = ret[0][1];
    *XN0 = ret[1][0];
    *YN0 = ret[1][1];
    *X0M = ret[2][0];
    *Y0M = ret[2][1];
    *XNM = ret[3][0];
    *YNM = ret[3][1];
  });
}

OPENZGY_API ZgyHandle oz_meta_indexcorners(ZgyHandle handle, double *X00, double *Y00, double* XN0, double* YN0, double* X0M, double* Y0M, double* XNM, double* YNM)
{
  return protect([&]() {
    OpenZGY::IZgyMeta::corners_t ret = ZgyMetaHandle::get(handle)->indexcorners();
    *X00 = ret[0][0];
    *Y00 = ret[0][1];
    *XN0 = ret[1][0];
    *YN0 = ret[1][1];
    *X0M = ret[2][0];
    *Y0M = ret[2][1];
    *XNM = ret[3][0];
    *YNM = ret[3][1];
  });
}

OPENZGY_API ZgyHandle oz_meta_annotcorners(ZgyHandle handle, double *X00, double *Y00, double* XN0, double* YN0, double* X0M, double* Y0M, double* XNM, double* YNM)
{
  return protect([&]() {
    OpenZGY::IZgyMeta::corners_t ret = ZgyMetaHandle::get(handle)->annotcorners();
    *X00 = ret[0][0];
    *Y00 = ret[0][1];
    *XN0 = ret[1][0];
    *YN0 = ret[1][1];
    *X0M = ret[2][0];
    *Y0M = ret[2][1];
    *XNM = ret[3][0];
    *YNM = ret[3][1];
  });
}

OPENZGY_API ZgyHandle oz_meta_bricksize(ZgyHandle handle, std::int64_t* ni_return, std::int64_t* nj_return, std::int64_t *nk_return)
{
  return protect([&]() {std::tie(*ni_return, *nj_return, *nk_return) = sizeToTuple(ZgyMetaHandle::get(handle)->bricksize());});
}

//  NOT WRAPPED. RARELY USEFUL. std::vector<size3i_t> brickcount()

OPENZGY_API ZgyHandle oz_meta_nlods(ZgyHandle handle, std::int32_t *value)
{
  return protect([&]() {*value = ZgyMetaHandle::get(handle)->nlods();});
}

OPENZGY_API ZgyHandle oz_meta_verid(ZgyHandle handle)
{
  // One of the very few methods that return a string.
  return protectReturn([&]() {
    return new ZgyStringHandle
      (ZgyMetaHandle::get(handle)->verid());});
}

#pragma endregion

#pragma region ZgyTools functions

/////////////////////////////////////////////////////////////////////////////
///   ZgyTools   ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_tools_annotToIndex(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->annotToIndex(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_tools_annotToWorld(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->annotToWorld(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_tools_indexToAnnot(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->indexToAnnot(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_tools_indexToWorld(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->indexToWorld(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_tools_worldToAnnot(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->worldToAnnot(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_tools_worldToIndex(ZgyHandle handle, double* x_return, double* y_return, double x, double y)
{
  return protect([&]() {std::tie(*x_return, *y_return) = doubleToTuple(ZgyMetaHandle::get(handle)->worldToIndex(std::array<double,2>{x, y}));});
}

OPENZGY_API ZgyHandle oz_meta_toString(ZgyHandle handle)
{
  // In C++ this writes to a stream. C# returns a string.
  return protectReturn([&]() {
    std::stringstream ss;
    auto hh = ZgyMetaHandle::get(handle);
    hh->dump(ss);
    return new ZgyStringHandle(ss.str());
  });
}

//  INTERNAL, NOT WRAPPED: oz_tools_transform1()

#pragma endregion

#pragma region SampleStatistics getters

/////////////////////////////////////////////////////////////////////////////
///   SampleStatistics   ////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_statistics_cnt(ZgyHandle handle, std::int64_t *value)
{
  // This class exposes data, not accessors.
  return protect([&]() {*value = ZgyStatisticsHandle::get(handle)->cnt;});
}

OPENZGY_API ZgyHandle oz_statistics_sum(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyStatisticsHandle::get(handle)->sum;});
}

OPENZGY_API ZgyHandle oz_statistics_ssq(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyStatisticsHandle::get(handle)->ssq;});
}

OPENZGY_API ZgyHandle oz_statistics_min(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyStatisticsHandle::get(handle)->min;});
}

OPENZGY_API ZgyHandle oz_statistics_max(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyStatisticsHandle::get(handle)->max;});
}

#pragma endregion

#pragma region SampleHistogram getters

/////////////////////////////////////////////////////////////////////////////
///   SampleHistogram   /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_histogram_samplecount(ZgyHandle handle, std::int64_t *value)
{
  // This class exposes data, not accessors.
  return protect([&]() {*value = ZgyHistogramHandle::get(handle)->samplecount;});
}

OPENZGY_API ZgyHandle oz_histogram_minvalue(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyHistogramHandle::get(handle)->minvalue;});
}

OPENZGY_API ZgyHandle oz_histogram_maxvalue(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyHistogramHandle::get(handle)->maxvalue;});
}

OPENZGY_API ZgyHandle oz_histogram_bins(ZgyHandle handle, std::int64_t *data, std::int32_t *actual_size_return, std::int32_t allocated_size)
{
  // The data is variable length.
  return protect([&]() {
    auto hist = ZgyHistogramHandle::get(handle);
    std::int32_t size = static_cast<std::int32_t>(hist->bins.size());
    *actual_size_return = size;
    if (allocated_size > 0) {
      memcpy(data, hist->bins.data(), sizeof(std::int64_t) * std::min(size, allocated_size));
    }
  });
}

#pragma endregion

#pragma region FileStatistics getters

/////////////////////////////////////////////////////////////////////////////
///   FileStatistics   //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Save some typing for the 18 getters returning int64.
#define LONG_GETTER(name) \
OPENZGY_API ZgyHandle oz_filestats_##name(ZgyHandle handle, std::int64_t *value) { \
return protect([&]() {*value = ZgyFileStatsHandle::get(handle)->name();});}

LONG_GETTER(fileVersion)
LONG_GETTER(fileSize)
LONG_GETTER(headerSize)
LONG_GETTER(dataStart)
LONG_GETTER(alphaNormalCount)
LONG_GETTER(alphaNormalSizePerEntry)
LONG_GETTER(alphaCompressedCount)
LONG_GETTER(alphaCcompressedSize)
LONG_GETTER(alphaMissingCount)
LONG_GETTER(alphaConstantCount)
LONG_GETTER(brickNormalCount)
LONG_GETTER(brickNormalSizePerEntry)
LONG_GETTER(brickCompressedCount)
LONG_GETTER(brickCompressedSize)
LONG_GETTER(brickMissingCount)
LONG_GETTER(brickConstantCount)
LONG_GETTER(usedSize)
LONG_GETTER(usedIfUncompressed)

OPENZGY_API ZgyHandle oz_filestats_compressionFactor(ZgyHandle handle, double *value)
{
  return protect([&]() {*value = ZgyFileStatsHandle::get(handle)->compressionFactor();});
}

OPENZGY_API ZgyHandle oz_filestats_isCompressed(ZgyHandle handle, int *value_intflag)
{
  return protect([&]() {*value_intflag = ZgyFileStatsHandle::get(handle)->isCompressed() ? 1 : 0;});
}

OPENZGY_API ZgyHandle oz_filestats_toString(ZgyHandle handle)
{
  // In C++ this writes to a stream. C# returns a string.
  // Also, remove the "prefix" second arg.
  return protectReturn([&]() {
    std::stringstream ss;
    ZgyFileStatsHandle::get(handle)->dump(ss, "");
    return new ZgyStringHandle(ss.str());
  });
}

#pragma endregion

#pragma region ZgyWriterArgs setters

/////////////////////////////////////////////////////////////////////////////
///   ZgyWriterArgs   ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_writerargs_filename(ZgyHandle handle, const char* name)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->filename(std::string(name ? name : ""));});
}

OPENZGY_API ZgyHandle oz_writerargs_iocontext(ZgyHandle handle, ZgyHandle iocontext)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->iocontext(ZgyIOContextHandle::get(iocontext));});
}

OPENZGY_API ZgyHandle oz_writerargs_zfp_compressor(ZgyHandle handle, float sqnr)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->zfp_compressor(sqnr);});
}

OPENZGY_API ZgyHandle oz_writerargs_zfp_lodcompressor(ZgyHandle handle, float sqnr)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->zfp_lodcompressor(sqnr);});
}

OPENZGY_API ZgyHandle oz_writerargs_size(ZgyHandle handle, std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->size(ni, nj, nk);});
}

OPENZGY_API ZgyHandle oz_writerargs_bricksize(ZgyHandle handle, std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->bricksize(ni, nj, nk);});
}

OPENZGY_API ZgyHandle oz_writerargs_datatype(ZgyHandle handle, ENUM value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->datatype((OpenZGY::SampleDataType)value);});
}

OPENZGY_API ZgyHandle oz_writerargs_datarange(ZgyHandle handle, float lo, float hi)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->datarange(lo, hi);});
}

OPENZGY_API ZgyHandle oz_writerargs_zunit(ZgyHandle handle, ENUM dimension, const char* name, double factor)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->zunit((OpenZGY::UnitDimension)dimension, std::string(name ? name : ""), factor);});
}

OPENZGY_API ZgyHandle oz_writerargs_hunit(ZgyHandle handle, ENUM dimension, const char* name, double factor)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->hunit((OpenZGY::UnitDimension)dimension, std::string(name ? name : ""), factor);});
}

OPENZGY_API ZgyHandle oz_writerargs_ilstart(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->ilstart(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_ilinc(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->ilinc(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_xlstart(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->xlstart(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_xlinc(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->xlinc(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_zstart(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->zstart(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_zinc(ZgyHandle handle, float value)
{
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->zinc(value);});
}

OPENZGY_API ZgyHandle oz_writerargs_corners(ZgyHandle handle, double X00, double Y00, double XN0, double YN0, double X0M, double Y0M, double XNM, double YNM)
{
  OpenZGY::ZgyWriterArgs::corners_t corners {{{ X00,Y00 }, { XN0,YN0 }, { X0M,Y0M }, { XNM,YNM }}};
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->corners(corners);});
}

OPENZGY_API ZgyHandle oz_writerargs_metafrom(ZgyHandle handle, ZgyHandle other)
{
  // Note, metafrom() expects a smartptr. get() would return a raw one.
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->metafrom(ZgyReaderHandle::getsmart(other, true));});
}

OPENZGY_API ZgyHandle oz_writerargs_merge(ZgyHandle handle, ZgyHandle other)
{
  // Note, merge() expects a reference, not a pointer.
  return protect([&]() {ZgyWriterArgsHandle::get(handle)->merge(*ZgyWriterArgsHandle::get(other));});
}

OPENZGY_API ZgyHandle oz_writerargs_toString(ZgyHandle handle)
{
  // In C++ this writes to a stream. C# returns a string.
  return protectReturn([&]() {
    std::stringstream ss;
    auto hh = ZgyWriterArgsHandle::get(handle);
    hh->dump(ss);
    return new ZgyStringHandle(ss.str());
  });
}

#pragma endregion

#pragma region SeismicStoreIOContext getters

/////////////////////////////////////////////////////////////////////////////
///   SeismicStoreIOContext   ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

OPENZGY_API ZgyHandle oz_ssiocontext_sdurl(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->sdurl(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_sdapikey(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->sdapikey(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_sdtoken(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->sdtoken(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_maxsize(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->maxsize(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_maxhole(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->maxhole(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_aligned(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->aligned(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_segsize(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->segsize(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_segsplit(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->segsplit(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_iothreads(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->iothreads(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_cputhreads(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->cputhreads(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_writethreads(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->writethreads(value);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_legaltag(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->legaltag(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_writeid(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->writeid(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_seismicmeta(ZgyHandle handle, const char* value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->seismicmeta(std::string(value ? value : ""));});
}

OPENZGY_API ZgyHandle oz_ssiocontext_setRoAfterWrite(ZgyHandle handle, int value_intflag)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->setRoAfterWrite(value_intflag != 0);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_forceRoBeforeRead(ZgyHandle handle, int value_intflag)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->forceRoBeforeRead(value_intflag != 0);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_forceRwBeforeWrite(ZgyHandle handle, int value_intflag)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->forceRwBeforeWrite(value_intflag != 0);});
}

OPENZGY_API ZgyHandle oz_ssiocontext_retryCount(ZgyHandle handle, int value)
{
  return protect([&]() {ZgyIOContextHandle::get(handle)->retryCount(value);});
}

/**
 * Beware object lifetimes.
 * See comments in namespace CallbackTypes.
 */
OPENZGY_API ZgyHandle oz_ssiocontext_logger(ZgyHandle handle, CallbackTypes::c_logger_t logger)
{
  return protect([&]() {
    LoggerForwarder loggerf(logger);
    ZgyIOContextHandle::get(handle)->logger(loggerf);
    });
}

OPENZGY_API ZgyHandle oz_ssiocontext_sdtokencb(ZgyHandle handle, CallbackTypes::c_tokencb_v3_t callback)
{
  return protect([&]() {
    TokenCallback3Forwarder callbackf(callback);
    ZgyIOContextHandle::get(handle)->sdtokencb(callbackf);
    });
}

OPENZGY_API ZgyHandle oz_ssiocontext_credentialsFrom(ZgyHandle handle, ZgyHandle utils_handle)
{
  return protect([&]() {
    ZgyIOContextHandle::get(handle)->credentialsFrom
      (ZgyUtilsHandle::getsmart(utils_handle, true));
    });
}

OPENZGY_API ZgyHandle oz_ssiocontext_toString(ZgyHandle handle)
{
  return protectReturn([&]() {return new ZgyStringHandle(ZgyIOContextHandle::get(handle)->toString()); });
}
#pragma endregion

} // extern "C"
