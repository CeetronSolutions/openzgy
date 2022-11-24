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
class LocalFileWindows : public FileCommon
{
  LocalFileWindows(const LocalFileWindows&) = delete;
  LocalFileWindows& operator=(const LocalFileWindows&) = delete;
public:
  LocalFileWindows(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  virtual ~LocalFileWindows();
  static std::shared_ptr<IFileADT> xx_make_instance(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  // Methods from IFileBase
  virtual void deleteFile(const std::string& name, bool missing_ok) const override;
  virtual std::string altUrl(const std::string& name) const override;
  virtual std::string idToken() const override;
  // Methods from IFileADT:
  virtual void xx_close() override;
  virtual std::int64_t xx_eof() const override;
  virtual std::vector<std::int64_t> xx_segments(bool complete) const override;
  virtual bool xx_iscloud() const override;
  virtual void xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_readv(const ReadList& requests, bool parallel_ok=false, bool immutable_ok=false, bool transient_ok=false, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual std::int64_t _real_eof() const;
private:
  int _fd;
  mutable std::mutex _mutex;
};

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> FileCommon -> LocalFileWindows   ///////////////////////////
/////////////////////////////////////////////////////////////////////////////

LocalFileWindows::LocalFileWindows(const std::string& filename, OpenMode mode, const IOContext*)
  : FileCommon(filename, mode)
  , _fd(-1)
  , _mutex()
{
  switch (mode) {
  case OpenMode::ReadOnly:
    _fd = _open(filename.c_str(), O_RDONLY | _O_BINARY, 0666);
    break;
  case OpenMode::ReadWrite:
    _fd = _open(filename.c_str(), O_RDWR | _O_BINARY, 0666);
    break;
  case OpenMode::Truncate:
    _fd = _open(filename.c_str(), O_RDWR | _O_BINARY | O_CREAT | O_TRUNC, 0666);
    break;
  case OpenMode::Closed:
  default:
    _fd = -2;
    break;
  }
  if (_fd == -1)
    throw OpenZGY::Errors::ZgyIoError(filename, errno);

  if (_fd >= 0) {
    _eof = static_cast<std::int64_t>(_lseeki64(_fd, 0, SEEK_END));
    (void)_lseeki64(_fd, 0, SEEK_SET);
    if (false)
      std::cout << "Opened file \"" << filename
      << "\" size " << std::hex << _eof << std::dec << "\n";
  }
}

LocalFileWindows::~LocalFileWindows()
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
LocalFileWindows::xx_make_instance(const std::string& filename, OpenMode mode, const IOContext *iocontext)
{
  if (filename.find("://") == std::string::npos) {
    auto file = std::shared_ptr<IFileADT>(new LocalFileWindows(filename, mode, iocontext));

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
LocalFileWindows::deleteFile(const std::string& name, bool missing_ok) const
{
  if (_unlink(name.c_str()) < 0 && errno != ENOENT)
    throw std::runtime_error("Cannot delete \"" + name + "\"");
}

std::string
LocalFileWindows::altUrl(const std::string& name) const
{
  return name;
}

std::string
LocalFileWindows::idToken() const
{
  return std::string();
}

/**
 * \details: Thread safety: No. All other operations must be completed first.
 */
void
LocalFileWindows::xx_close()
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

  switch (mode) {

  default:
  case OpenMode::Closed:
    break;

  case OpenMode::ReadOnly:
  case OpenMode::ReadWrite:
  case OpenMode::Truncate:
    if (_close(_fd) < 0)
      throw OpenZGY::Errors::ZgyIoError(_name, errno);
    _fd = -2;
    break;
  }

  _fd = -2;
  _name = std::string();
  _rtimer.reset();
  _wtimer.reset();
  _mtimer.reset();
}

/**
 * \details: Thread safety: Yes, by locking.
 */
std::int64_t
LocalFileWindows::xx_eof() const
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
LocalFileWindows::xx_segments(bool /*complete*/) const
{
  return std::vector<std::int64_t>{this->xx_eof()};
}

/**
 * \details: Thread safety: Yes.
 */
bool
LocalFileWindows::xx_iscloud() const
{
  return false;
}

/**
 * \details: Thread safety: Yes, by setting a lock.
 */
void
LocalFileWindows::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_rtimer);
  _validate_read(data, offset, size, xx_eof(), _mode);
  // On windows, _read deals with int32 size so the cast warning
  // technically correct. In practice we never read that much in one go.
  if (size > std::numeric_limits<int>::max())
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Too large read");
  SimpleTimerEx mm(*_mtimer);
  std::lock_guard<std::mutex> lk(_mutex); // make seek + read atomic
  mm.stop();
#if XXX_DISABLE_ALL_IO
  memset(data, 0, size);
  *(float*)data = 1.0f;
  return;
#endif
  _lseeki64(_fd, offset, SEEK_SET);
  std::int64_t nbytes = static_cast<std::int64_t>(_read(_fd, data, static_cast<int>(size)));
  _check_short_read(offset, size, nbytes);
}

/**
 * \details: Thread safety: Yes, by setting locks.
 */
void
LocalFileWindows::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
#if XXX_DISABLE_ALL_IO
  static std::shared_ptr<char> spare_buffer;
  static std::int64_t spare_size = 0;
  for (const ReadRequest& r : requests) {
    // Very much not thread safe, but this is a huge hack anyway.
    if (!spare_buffer || spare_size < r.size) {
      spare_buffer.reset(new char[r.size], std::default_delete<char[]>());
    }
    this->LocalFileWindows::xx_read(spare_buffer.get(), r.offset, r.size, usagehint);
    _deliver(r.delivery, spare_buffer, 0, r.size, transient_ok);
  }
  // TODO, who is holding on to the buffers?
  return;
#endif
  // Note that I need the xx_read() in this class, not any overrides.
  // If xx_read() is overridden then whoever did that wouldn't expect
  // xx_readv() to change. The fact that I choose to implement one in
  // terms of the other is an implementation detail.
  for (const ReadRequest& r : requests) {
    std::shared_ptr<char> data(new char[r.size], std::default_delete<char[]>());
    this->LocalFileWindows::xx_read(data.get(), r.offset, r.size, usagehint);
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
LocalFileWindows::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_wtimer);
#if XXX_DISABLE_ALL_IO
  _eof = std::max(_eof, offset + size);
  return;
#endif
  _validate_write(data, offset, size, _mode);
  // On windows, _write deals with int32 size so the compiler warning is
  // technically correct. In practice we never read that much in one go.
  if (size > std::numeric_limits<int>::max())
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Too large write");
  if (false)
    std::cout << "xx_write(*, " << std::hex
              << offset << ", " << size << ", hint=" << (int)usagehint
              << std::dec << ")\n";
  SimpleTimerEx mm(*_mtimer);
  std::lock_guard<std::mutex> lk(_mutex); // make seek + write + _eof atomic
  mm.stop();
  _lseeki64(_fd, offset, SEEK_SET);
  std::int64_t nbytes = static_cast<std::int64_t>(_write(_fd, data, static_cast<int>(size)));
  if (nbytes < 0)
      throw OpenZGY::Errors::ZgyIoError(_name, errno);
  _eof = std::max(_eof, offset + nbytes);
  if (nbytes != size)
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Short write");
}

/**
 * \details: Thread safety: Yes.
 */
std::int64_t
LocalFileWindows::_real_eof() const
{
  struct _stat64 st;
  if (::_fstat64(_fd, &st) < 0)
    return -1;
  return static_cast<std::int64_t>(st.st_size);
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
      FileFactory::instance().add_factory(LocalFileWindows::xx_make_instance);
    }
  } dummy;
} // anonymous namespace for registration

} // namespace

#endif // Entire file is windows only
