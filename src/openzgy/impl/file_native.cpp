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

#ifdef _WIN32 // Entire file is windows only

// TODO-Low: Get rid of this temporary kludge on Windows.
#define _CRT_SECURE_NO_WARNINGS 1

#include "file.h"
#include "../exception.h"
#include "fancy_timers.h"
#include "environment.h"
#include "file_parallelizer.h"
#include "file_performance.h"

#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <io.h>
#include <mutex>

#ifdef _WIN32 // Entire file is Windows only
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#endif

#define XXX_DISABLE_ALL_IO 0

using OpenZGY::IOContext;
namespace InternalZGY {
#if 0
}
#endif

/**
 * \file: file_local.cpp
 * \brief Low level I/O, regular files.
 */

/**
 * For debugging, ensure all local disk buffers are written to disk
 * before closing. This makes timing results more reproducible. N/A
 * for cloud I/O.
 *
 * ZgyTool on Linux currently does its own sync (of the entire system).
 * In that situation it is better to enable fsync here, and disregard
 * the sync reported by the tool. Because that will now report only the
 * time needed to sync unrelated buffers.
 *
 * If the sync time is suspiciously slow then this suggests that fsync
 * is called or implied somewhere else.
 */
static int enable_fsync()
{
  static int enable = Environment::getNumericEnv("OPENZGY_ENABLE_FSYNC", 0);
  return enable;
}

/**
 * Choose between posix style _open, _read, etc. and Windows
 * CreateFile, ReadFile, etc.
 */
static int use_windows_fileaccess()
{
  static int enable = Environment::getNumericEnv("OPENZGY_WINDOWS_FILEACCESS", 0);
  return enable > 0;
}

/**
 * Thread safety when used for reading:
 * Designed to be thread safe as no internal data structures should change
 * after the file has been opened. Any lazy-evaluated information needs
 * to do appropriate locking.
 *
 * Thread safety when used for writing:
 * Yes, but caller won't mak use of this because the high level design
 * states that writes (i.e. calls  from ZgyWriter) are not thread safe.
 * Besides, other file backends might not be thread safe and those might
 * not be so easy to change.
 *
 * Thread safety when closing a file:
 * Not thread safe.
 *
 * The class is noncopyable with copy and assign method deleted.
 * Not that the users could easily copy it anyway, as the class should
 * always be used via the IFileADT interface.
 */
class LocalFileNativeWindows : public FileCommon
{
  LocalFileNativeWindows(const LocalFileNativeWindows&) = delete;
  LocalFileNativeWindows& operator=(const LocalFileNativeWindows&) = delete;
public:
  LocalFileNativeWindows(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  virtual ~LocalFileNativeWindows();
  static std::shared_ptr<IFileADT> xx_make_instance(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  // Methods from IFileBase
  void deleteFile(const std::string& name, bool missing_ok) const override;
  std::string altUrl(const std::string& name) const override;
  std::string idToken() const override;
  // Methods from IFileADT:
  void xx_close() override;
  std::int64_t xx_eof() const override;
  std::vector<std::int64_t> xx_segments(bool complete) const override;
  bool xx_iscloud() const override;
  void xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  void xx_readv(const ReadList& requests, bool parallel_ok=false, bool immutable_ok=false, bool transient_ok=false, UsageHint usagehint=UsageHint::Unknown) override;
  void xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  std::int64_t _real_eof() const override;
private:
  std::int64_t _get_filesize() const;
private:
  HANDLE _handle;
  mutable std::mutex _mutex;
};

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> FileCommon -> LocalFileNativeWindows   ///////////////////////////
/////////////////////////////////////////////////////////////////////////////

LocalFileNativeWindows::LocalFileNativeWindows(const std::string& filename, OpenMode mode, const IOContext*)
  : FileCommon(filename, mode)
  , _handle(INVALID_HANDLE_VALUE)
  , _mutex()
{
  // Function specific to native vs. posix version.
  ::SetLastError(ERROR_SUCCESS);
  _handle = INVALID_HANDLE_VALUE;
  switch (mode) {
  case OpenMode::ReadOnly:
    // If opening for read, we don't mind that some other process
    // opens it for write later. Thus the FILE_SHARE_WRITE in that case.
    _handle = ::CreateFileA(filename.c_str(),
      GENERIC_READ,
      FILE_SHARE_READ | FILE_SHARE_WRITE,
      NULL,
      OPEN_EXISTING,
      FILE_ATTRIBUTE_NORMAL,
      NULL);
    break;

  case OpenMode::ReadWrite:
    // If we are opening it for write, we specify FILE_SHARE_READ
    // meaning we allow others to open it for reading but not writing.
    _handle = ::CreateFileA(filename.c_str(),
      GENERIC_READ | GENERIC_WRITE,
      FILE_SHARE_READ,
      NULL,
      OPEN_EXISTING,
      FILE_ATTRIBUTE_NORMAL,
      NULL);
    break;

  case OpenMode::Truncate:
    _handle = ::CreateFileA(filename.c_str(),
      GENERIC_READ | GENERIC_WRITE,
      FILE_SHARE_READ,
      NULL,
      CREATE_ALWAYS,
      FILE_ATTRIBUTE_NORMAL,
      NULL);
    break;

  case OpenMode::Closed:
  default:
    _handle = INVALID_HANDLE_VALUE;
    break;
  }

  DWORD err = ::GetLastError();
  // This result from CREATE_ALWAYS or OPEN_ALWAYS
  // is informational, not an actual error.
  //if (err == ERROR_ALREADY_EXISTS)
  //  std::cerr << "Opened existing ZGY file\n";

  if (mode != OpenMode::Closed && _handle == INVALID_HANDLE_VALUE)
    throw OpenZGY::Errors::ZgyWindowsError(filename, err);

  _eof = _get_filesize();

  if (false)
    std::cerr << "Opened "
    << (err == ERROR_ALREADY_EXISTS ? "existing" : "new")
    << " file \"" << filename
    << "\" size " << std::hex << _eof << std::dec << "\n";
}

LocalFileNativeWindows::~LocalFileNativeWindows()
{
  if (_mode != OpenMode::Closed) {
    try {
      xx_close();
    }
    catch (const std::exception& ex) {
      // The calling layer is supposed to do an explicit xx_close()
      // so it can catch and handle exceptions. This blind catch is
      // just a desperate attempt to avoid an application crash.
      std::cerr << "EXCEPTION closing file: " << ex.what() << std::endl;
    }
  }
}

std::shared_ptr<IFileADT>
LocalFileNativeWindows::xx_make_instance(const std::string& filename, OpenMode mode, const IOContext *iocontext)
{
  if (filename.find("://") == std::string::npos && use_windows_fileaccess()) {
    std::cerr << "Using LocalFileNativeWindows to access \"" << filename << "\"\n";
    auto file = std::shared_ptr<IFileADT>(new LocalFileNativeWindows(filename, mode, iocontext));

    // This is a no-op unless enabled by enviroment variables
    file = FileWithPerformanceLogger::inject(file, filename);

    // The following is to parallelize only decompression and copy-out.
    // Unlike the Linux case there is no parallel read at the lowest level.
    // Which would have automatically parallelized those two.
    // Unlike SeismicStoreFile the thread count isn't configuarble.
    // Just set a very high count. It gets clipped to what OpenMP suggests.
    // Caveat: If application makes a lot of threads itself then even the
    // OpenMP default might be to much.
    file = FileParallelizer::inject(file, 256);

    return file;
  }
  else {
    return std::shared_ptr<IFileADT>();
  }
}

void
LocalFileNativeWindows::deleteFile(const std::string& name, bool missing_ok) const
{
  // Function specific to native vs. posix version.
  if (0 == ::DeleteFileA(name.c_str()))
    throw OpenZGY::Errors::ZgyWindowsError(name, ::GetLastError());
}

std::string
LocalFileNativeWindows::altUrl(const std::string& name) const
{
  return name;
}

std::string
LocalFileNativeWindows::idToken() const
{
  return std::string();
}

/**
 * \details: Thread safety: No. All other operations must be completed first.
 */
void
LocalFileNativeWindows::xx_close()
{
  if (_mode == OpenMode::Closed) {
    // Note: I might "be nice" to the application and simply ignore a duplicate
    // close or a close on a file that was never open in the first place.
    // In that case I should probably check using an atomic_flag.
    // But if the application issues extraneous xx_close(), let alone multiple
    // concurrent calls to xx_close(), this is a bug. That may indicate there
    // is something else wrong as well.
    throw OpenZGY::Errors::ZgyUserError("Attemping to close a file twice.");
  }

  OpenMode mode = _mode;
  _mode = OpenMode::Closed;     // In case we throw.

  // Begin specific to native vs. posix version.
  switch (mode) {

  default:
  case OpenMode::Closed:
    break;

  case OpenMode::ReadOnly:
  case OpenMode::ReadWrite:
  case OpenMode::Truncate:
    if (mode != OpenMode::ReadOnly && enable_fsync() > 0) {
      SimpleTimerEx mm(*_synctimer);
      (void)::FlushFileBuffers(_handle); // errors are not fatal.
    }
    if (0 == ::CloseHandle(_handle)) {
      _handle = INVALID_HANDLE_VALUE;
      throw OpenZGY::Errors::ZgyWindowsError(_name, ::GetLastError());
    }
    break;
  }
  _handle = INVALID_HANDLE_VALUE;
  // End specific to native vs. posix version.

  _name = std::string();
  _sync3timer.reset();
  _synctimer.reset();
  _rtimer.reset();
  _wtimer.reset();
  _mtimer.reset();
}

/**
 * \details: Thread safety: Yes, by locking.
 */
std::int64_t
LocalFileNativeWindows::xx_eof() const
{
  SimpleTimerEx mm(*_mtimer);
  std::lock_guard<std::mutex> lk(_mutex); // protect _eof
  mm.stop();
  return this->_eof;
}

/**
 * \details: Thread safety: Yes, called method is thread safe.
 */
std::vector<std::int64_t>
LocalFileNativeWindows::xx_segments(bool /*complete*/) const
{
  return std::vector<std::int64_t>{this->xx_eof()};
}

/**
 * \details: Thread safety: Yes.
 */
bool
LocalFileNativeWindows::xx_iscloud() const
{
  return false;
}

/**
 * \details: Thread safety: Yes, by setting a lock.
 */
void
LocalFileNativeWindows::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
#if XXX_DISABLE_ALL_IO
  memset(data, 0, size);
  *(float*)data = 1.0f;
  return;
#endif
  SimpleTimerEx tt(*_rtimer);
  _validate_read(data, offset, size, xx_eof(), _mode);
  // On windows, _read deals with int32 size so the compiler warning is
  // technically correct. In practice we never read that much in one go.
  if (size > std::numeric_limits<int>::max() || size > std::numeric_limits<DWORD>::max())
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Too large read");
  if (size <= 0)
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Bytes to write must be > 0");
  // Begin specific to native vs. posix version.
  OVERLAPPED ovlp;
  memset(&ovlp, 0, sizeof(ovlp));
  ovlp.Offset = (DWORD)(offset % (std::int64_t(1) << 32));
  ovlp.OffsetHigh = (DWORD)(offset >> 32);
  DWORD nRead(0);
  if (!::ReadFile(_handle, data, (DWORD)size, &nRead, &ovlp))
    throw OpenZGY::Errors::ZgyWindowsError(_name, ::GetLastError());
  std::int64_t nbytes = static_cast<std::int64_t>(nRead);
  // End specific to native vs. posix version.
  _check_short_read(offset, size, nbytes);
}

/**
 * \details: Thread safety: Yes, by setting locks.
 */
void
LocalFileNativeWindows::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
#if XXX_DISABLE_ALL_IO
  static std::shared_ptr<char> spare_buffer;
  static std::int64_t spare_size = 0;
  for (const ReadRequest& r : requests) {
    // Very much not thread safe, but this is a huge hack anyway.
    if (!spare_buffer || spare_size < r.size) {
      spare_buffer.reset(new char[r.size], std::default_delete<char[]>());
    }
    this->LocalFileNativeWindows::xx_read(spare_buffer.get(), r.offset, r.size, usagehint);
    _deliver(r.delivery, spare_buffer, 0, r.size, transient_ok);
  }
  // TODO, who is holding on to the buffers?
  return;
#endif
  // Note that I need the xx_read() in this class, not any overrides.
  // If xx_read() is overridden then whoever did that wouldn't expect
  // xx_readv() to change. The fact that I choose to implement one in
  // terms of the other is an implementation detail.
  // TODO-Performance: True async I/O. Won't help Petrel because it
  // currently reads brick or brick-column (in genlod) at a time.
  // And might not help that much, anyway.
  for (const ReadRequest& r : requests) {
    std::shared_ptr<char> data(new char[r.size], std::default_delete<char[]>());
    // Next line specific to native vs. posix version.
    this->LocalFileNativeWindows::xx_read(data.get(), r.offset, r.size, usagehint);
    _deliver(r.delivery, data, 0, r.size, transient_ok);
  }
}

/**
 * \details: Thread safety: Yes.
 * OpenZGY is in general not thread safe when writing, but in this
 * low level function it doesn't cost much to synchronize _eof
 * both here and the places it is read.
 */
void
LocalFileNativeWindows::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_wtimer);
#if XXX_DISABLE_ALL_IO
  _eof = std::max(_eof, offset + size);
  return;
#endif
  _validate_write(data, offset, size, _mode);
  // On windows, _write deals with int32 size so the compiler warning is
  // technically correct. In practice we never read that much in one go.
  if (size > std::numeric_limits<int>::max() || size > std::numeric_limits<DWORD>::max())
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Too large write");
  if (size <= 0)
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Bytes to write must be > 0");
  if (false)
    std::cout << "xx_write(*, " << std::hex
              << offset << ", " << size << ", hint=" << (int)usagehint
              << std::dec << ")\n";
  // Begin specific to native vs. posix version.
  OVERLAPPED ovlp;
  memset(&ovlp, 0, sizeof(ovlp));
  ovlp.Offset = (DWORD)(offset % (std::int64_t(1) << 32));
  ovlp.OffsetHigh = (DWORD)(offset >> 32);
  DWORD nWritten(0);
  if (!::WriteFile(_handle, data, static_cast<DWORD>(size), &nWritten, &ovlp))
    throw OpenZGY::Errors::ZgyWindowsError(_name, ::GetLastError());
  std::int64_t nbytes = static_cast<std::int64_t>(nWritten);
  {
    // Unlike _write, no need to lock the write call itself.
    // But _eof may still need protection.
    SimpleTimerEx mm(*_mtimer);
    std::lock_guard<std::mutex> lk(_mutex); // protect _eof
    mm.stop();
    _eof = std::max(_eof, offset + nbytes);
  }
  // End specific to native vs. posix version.
  if (nbytes != size)
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Short write");
  static int periodical{ 1 };
  if (enable_fsync() >= 3 && (periodical%1000) == 0) {
    tt.stop();
    SimpleTimerEx uu(*_sync3timer);
    // Next line specific to native vs. posix version.
    ::FlushFileBuffers(_handle);
    _rtimer->print();
    _wtimer->print();
    _synctimer->print();
    _sync3timer->print();
  }
  ++periodical;
}

/**
 * \details: Thread safety: Yes.
 */
std::int64_t
LocalFileNativeWindows::_real_eof() const
{
  return _get_filesize();
}

/**
 * \brief Non-virtual get file size, callable from constructor.
 * \details: Returns -1 instead of throwing on error.
 * Thread safety: Yes.
 */
std::int64_t 
LocalFileNativeWindows::_get_filesize() const
{
  // Function specific to native vs. posix version.
  if (_handle == INVALID_HANDLE_VALUE)
    return -1;
  DWORD dwHigh(0), dwLow(0);
  ::SetLastError(ERROR_SUCCESS);
  dwLow = ::GetFileSize(_handle, &dwHigh);
  DWORD err = ::GetLastError();
  if (err != ERROR_SUCCESS)
    return -1;
  return (std::int64_t(dwHigh) << 32) + dwLow;
}

namespace {
  /**
   * \details: Thread safety: Yes. add_factory() is synchronized.
   */
  class Register
  {
  public:
    Register()
    {
      FileFactory::instance().add_factory(LocalFileNativeWindows::xx_make_instance);
    }
  } dummy;
} // anonymous namespace for registration

} // namespace

#endif // Entire file is windows only
