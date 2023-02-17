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

#ifndef _WIN32 // Entire file is linux only

#include "file.h"
#include "../exception.h"
#include "timer.h"
#include "fancy_timers.h"
#include "environment.h"
#include "mtguard.h"
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
#include <unistd.h>
#include <mutex>
#include <atomic>
#include <omp.h>

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
class LocalFileLinux : public FileCommon
{
  LocalFileLinux(const LocalFileLinux&) = delete;
  LocalFileLinux& operator=(const LocalFileLinux&) = delete;
public:
  LocalFileLinux(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  virtual ~LocalFileLinux();
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
  int _fd;
  mutable std::mutex _mutex;
};

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> FileCommon -> LocalFileLinux   /////////////////////////////
/////////////////////////////////////////////////////////////////////////////

LocalFileLinux::LocalFileLinux(const std::string& filename, OpenMode mode, const IOContext*)
  : FileCommon(filename, mode)
  , _fd(-1)
  , _mutex()
{
  switch (mode) {
  case OpenMode::ReadOnly:
    _fd = ::open(filename.c_str(), O_RDONLY, 0666);
    break;
  case OpenMode::ReadWrite:
    _fd = ::open(filename.c_str(), O_RDWR, 0666);
    break;
  case OpenMode::Truncate:
    _fd = ::open(filename.c_str(), O_RDWR|O_CREAT|O_TRUNC, 0666);
    break;
  case OpenMode::Closed:
  default:
    _fd = -2;
    break;
  }
  if (_fd == -1)
    throw OpenZGY::Errors::ZgyIoError(filename, errno);

  if (_fd >= 0) {
    _eof = static_cast<std::int64_t>(::lseek(_fd, 0, SEEK_END));
    (void)::lseek(_fd, 0, SEEK_SET);
    if (false)
      std::cout << "Opened file \"" << filename
      << "\" size " << std::hex << _eof << std::dec << "\n";
  }
}

LocalFileLinux::~LocalFileLinux()
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
LocalFileLinux::xx_make_instance(const std::string& filename, OpenMode mode, const IOContext *iocontext)
{
  if (filename.find("://") == std::string::npos) {
    auto file = std::shared_ptr<IFileADT>(new LocalFileLinux(filename, mode, iocontext));
    // This is a no-op unless enabled by enviroment variables
    file = FileWithPerformanceLogger::inject(file, filename);

    // This is for ad-hoc testing ONLY. Enable the parallelizer as the
    // windows reader does. On Linux (i.e. in this file) it gives
    // better results to parallelize inside xx_readv.
    // Remember to add #include "file_parallelizer.h"
    //file = FileParallelizer::inject(file, 256);

    return file;
  }
  else
    return std::shared_ptr<IFileADT>();
}

void
LocalFileLinux::deleteFile(const std::string& name, bool missing_ok) const
{
  if (::unlink(name.c_str()) < 0 && errno != ENOENT)
    throw std::runtime_error("Cannot delete \"" + name + "\"");
}

std::string
LocalFileLinux::altUrl(const std::string& name) const
{
  return name;
}

std::string
LocalFileLinux::idToken() const
{
  return std::string();
}

/**
 * \details: Thread safety: No. All other operations must be completed first.
 */
void
LocalFileLinux::xx_close()
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
    if (mode != OpenMode::ReadOnly && enable_fsync() > 0) {
      SimpleTimerEx mm(*_synctimer);
      (void)fsync(_fd); // errors are not fatal.
    }
    if (::close(_fd) < 0)
      throw OpenZGY::Errors::ZgyIoError(_name, errno);
    _fd = -2;
    break;
  }

  _fd = -2;
  _name = std::string();
  _synctimer.reset();
  _rtimer.reset();
  _wtimer.reset();
  _mtimer.reset();
}

/**
 * \details: Thread safety: Yes, by locking.
 */
std::int64_t
LocalFileLinux::xx_eof() const
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
LocalFileLinux::xx_segments(bool /*complete*/) const
{
  return std::vector<std::int64_t>{this->xx_eof()};
}

/**
 * \details: Thread safety: Yes.
 */
bool
LocalFileLinux::xx_iscloud() const
{
  return false;
}

/**
 * \details: Thread safety: Yes, assuming that the linux ::pread is thread safe.
 */
void
LocalFileLinux::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_rtimer);
  _validate_read(data, offset, size, xx_eof(), _mode);
  ssize_t nbytes = ::pread(_fd, data, size, offset);
  _rtimer->addBytesRead(nbytes);
  _check_short_read(offset, size, nbytes);
}

/**
 * \details: Thread safety: Yes, assuming that the linux ::pread is thread safe.
 *
 * If the caller passes parallel_ok=true this means the caller allows and
 * even prefers that we deliver each request on a different thread. This
 * parallelization comes in addition to allowing multiple reads in parallel
 * at the OpenZGY API level.
 *
 * Caveat: Consider carefully whether you want both. If the
 * application uses OpenMP for multi threading then by default nested
 * parallel regions are disabled. You can change this. If the
 * application uses some other mechanism than OpenMP used here might
 * not realize that it is creating nested loops. Or maybe it does, if
 * it uses an application-wide thread pool?
 *
 * TODO-Test:
 * Caveat: Since finalize() is single threaded then it should probably
 * enable parallel here. One problem is that the application might
 * still be inside an OpenMP loop, using a lock to make sure that
 * finalize() runs unmolested. OpenMP will still see it is inside a
 * parallel region so it might refuse to make one here.
 */
void
LocalFileLinux::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
  // Note that I need the xx_read() in this class, not any overrides.
  // If xx_read() is overridden then whoever did that wouldn't expect
  // xx_readv() to change. The fact that I choose to implement one in
  // terms of the other is an implementation detail.
  // Note that the delivery function can retain a reference to the data.
  // This is allowed as long as the data is still short lived. If not then
  // this isn't a disaster due to the eviction code in _allocate().
  if (!parallel_ok || requests.size() < 2) {
    for (const ReadRequest& r : requests) {
      std::shared_ptr<void> data = _allocate(r.size);
      this->LocalFileLinux::xx_read(data.get(), r.offset, r.size, usagehint);
      _deliver(r.delivery, data, 0, r.size, transient_ok);
    }
  }
  else {
    // OpenMP needs signed loop variable on windows.
    const std::int64_t requestcount = requests.size();

    // It is pointless to use more threads than we have requests.
    const int threadcount = static_cast<int>
      (std::min(requestcount,static_cast<std::int64_t>(omp_get_max_threads())));

    MTGuard guard("local-read", threadcount);
#pragma omp parallel num_threads(threadcount)
    {
      std::int64_t datasize = 0;
      std::shared_ptr<char> data;
#pragma omp for
      for (std::int64_t ii=0; ii<requestcount; ++ii) {
        const ReadRequest& r = requests[ii];
        if (datasize < r.size || !data || !data.unique()) {
          datasize = 0;
          data.reset(new char[r.size], std::default_delete<char[]>());
          datasize = r.size;
        }
        guard.run([&](){
          this->LocalFileLinux::xx_read(data.get(), r.offset, r.size, usagehint);
          _deliver(r.delivery, data, 0, r.size, transient_ok);
        });
      }
    }
    guard.finished();
  }
}

/**
 * \details: Thread safety: Yes.
 * OpenZGY is in general not thread safe when writing, but in this
 * low level function it doesn't cost much to synchronize _eof
 * both here and the places it is read.
 */
void
LocalFileLinux::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_wtimer);
  _validate_write(data, offset, size, _mode);
  if (false)
    std::cout << "xx_write(*, " << std::hex
              << offset << ", " << size << ", hint=" << (int)usagehint
              << std::dec << ")\n";
  ssize_t nbytes = ::pwrite(_fd, data, size, offset);
  if (nbytes < 0)
      throw OpenZGY::Errors::ZgyIoError(_name, errno);
  _wtimer->addBytesWritten(nbytes);
  SimpleTimerEx mm(*_mtimer);
  std::lock_guard<std::mutex> lk(_mutex); // protect _eof
  mm.stop();
  _eof = std::max(_eof, offset + nbytes);
  if (nbytes != size)
    throw OpenZGY::Errors::ZgyInternalError(_name + ": Short write");
}

/**
 * \details: Thread safety: Yes. Uses xx_eof().
 * And should in any case be changed to call ::stat().
 */
std::int64_t
LocalFileLinux::_real_eof() const
{
  struct stat st;
  if (::fstat(_fd, &st) < 0)
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
      FileFactory::instance().add_factory(LocalFileLinux::xx_make_instance);
    }
  } dummy;
} // anonymous namespace for registration

} // namespace

#endif // Entire file is linux only
