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

#include "file.h"
#include "environment.h"
#include "fancy_timers.h"
#include "../exception.h"

#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <iostream>
#include <sstream>
#include <mutex>
#include <algorithm>
#include <random>

using OpenZGY::IOContext;
namespace InternalZGY {
#if 0
}
#endif

namespace {
  /**
   * Return non-zero if the private memory allocator should be used in certain
   * circumstances. The default is 1 which will enable those places which will
   * most likey benefit and which appear safe. Set to 0 to get an idea of how
   * large the saving is. Note that this might result in disappointment, aa the
   * system malloc() is fairly efficient already.
   */
  static int
  malloc_shortcut()
  {
    static int enable = Environment::getNumericEnv("OPENZGY_MALLOC_SHORTCUT", 1);
    return enable;
  }
}

/////////////////////////////////////////////////////////////////////////////
//    FileADT (base class)  /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

IFileBase::~IFileBase()
{
}

IFileADT::~IFileADT()
{
}

FileADT::~FileADT()
{
}

/**
 * \brief Human readable number.
 */
std::string
FileADT::_nice(std::int64_t n)
{
  return SummaryPrintingTimerEx::niceSize(n);
}

void
FileADT::_validate_read(void *data, std::int64_t offset, std::int64_t size, std::int64_t eof, OpenMode mode)
{
  switch (mode) {
  case OpenMode::ReadOnly:
  case OpenMode::ReadWrite:
  case OpenMode::Truncate:
    break;
  case OpenMode::Closed:
  default:
    throw OpenZGY::Errors::ZgyUserError("The file is not open for reading.");
  }

  // These test are more relevant in the Python version because there
  // it can also do a type check on the arguments.
  if (data == nullptr)
    throw OpenZGY::Errors::ZgyUserError("Trying to read into null buffer.");

  if (offset < 0)
    throw OpenZGY::Errors::ZgyUserError("Trying to read at negative offset.");

  if (size < 1)
    throw OpenZGY::Errors::ZgyUserError("Trying to read zero or negative bytes.");

  // The next one night be an internal error or a corrupted file,
  // but let's report only the immediate error and not try to guess.
  if (offset + size > eof)
    throw OpenZGY::Errors::ZgyEndOfFile("Offset " + _nice(offset) +
                                " size " + _nice(size) +
                                " is past EOF at " + _nice(eof));
}

void
FileADT::_validate_write(const void *data, std::int64_t offset, std::int64_t size, OpenMode mode)
{
  switch (mode) {
  case OpenMode::ReadWrite:
  case OpenMode::Truncate:
    break;
  case OpenMode::Closed:
  case OpenMode::ReadOnly:
  default:
    throw OpenZGY::Errors::ZgyUserError("The file is not open for write.");
  }

  if (offset < 0)
    throw OpenZGY::Errors::ZgyUserError("Trying to write to negative offset.");

  if (size < 1 || data == nullptr)
    throw OpenZGY::Errors::ZgyUserError("Trying to write zero or negative bytes.");
}

void
FileADT::_validate_readv(const ReadList& requests, std::int64_t eof, OpenMode mode)
{
  switch (mode) {
  case OpenMode::ReadOnly:
  case OpenMode::ReadWrite:
  case OpenMode::Truncate:
    break;
  case OpenMode::Closed:
  default:
    throw OpenZGY::Errors::ZgyUserError("The file is not open for reading.");
  }
  for (const ReadRequest& rr : requests) {
    if (rr.offset < 0)
      throw OpenZGY::Errors::ZgyUserError("Trying to read at negative offset.");
    if (rr.size < 1)
      throw OpenZGY::Errors::ZgyUserError("Trying to read zero or negative bytes.");
    if (rr.offset + rr.size > eof)
      throw OpenZGY::Errors::ZgyEndOfFile("Offset " + _nice(rr.offset) +
                                " size " + _nice(rr.size) +
                                " is past EOF at " + _nice(eof));
  }
}

/**
 * Convenience function to invoke a delivery functor with optional
 * pointer arithmetic, delivering just a part of the buffer.
 *
 * Optionally check that the called function did not retain a pointer
 * to the data if it promised not to do that.
 *
 * If the contract about not retaining references is broken then raise
 * an exception. Even if the code happens to have allocate a unique
 * buffer so it doesn't really care. Note: In some cases, e.g. if a
 * proper cache is involved, the functor might end up retaining a
 * pointer aliased to the entire cache. Which would be fatal and would
 * warrant an abort(). The code here will still just throw an
 * exception, though.
 *
 * If the functor states that it needs to retain a pointer then make
 * sure it gets a smart pointer that is aliased to the entire buffer.
 */

/**
 * For debugging, but might be kept in production code as an assert.
 *
 * Equivalent to calling fn(data, size) but abuse the smart_ptr
 * mechanism to check whether the functor kept references to the
 * data. And throw an exception if this is not safe.
 *
 * Note that this is meant to be called from inside _deliver().
 * It is not a replacement for _deliver becaise it is missing the
 * offset argument so it doesn't handle pointer arithmetic.
 * And this function doesn't allow fn to be empty.
 *
 * The function cannot simply check the reference count before
 * and after calling the functor because the data smart pointer
 * can be asjusted by another thread. So there needs to be some
 * kind of indirect pointer that is unique to this function but
 * still keeps a reference to the data and releases it when the
 * indirect pointer is released.
 */
void
FileADT::_deliver_now(
     const ReadRequest::delivery_t& fn,
     const ReadRequest::data_t& data,
     std::int64_t size,
     bool transient)
{
  //static std::atomic<int> pending(0);
  //pending.fetch_add(1);
  ReadRequest::data_t* keep = new ReadRequest::data_t(data);
  ReadRequest::data_t  fake(data.get(), [keep](const void* p) {
                            //std::cerr << "Delivery functor released pointer\n";
                            //if (pending.fetch_sub(1) == 1)
                            //  std::cerr << "All data references released\n";
                            delete keep;
                       });
  fn(fake, size);
  auto retained = fake.use_count() - 1;
  fake.reset();
  std::string msg = ("Delivery functor retained " +
                     std::to_string(retained) +
                     std::string(retained==1 ? " reference." : " references."));
  if (retained != 0 && transient) {
    std::cerr << (msg + " FATAL ERROR.\n");
    throw OpenZGY::Errors::ZgyInternalError(msg);
  }
  else {
    //std::cerr << (msg + "\n");
  }
}

void
FileADT::_deliver(
     const ReadRequest::delivery_t& fn,
     const ReadRequest::data_t& data,
     std::int64_t offset,
     std::int64_t size,
     bool transient)
{
  if (!data)
    throw OpenZGY::Errors::ZgyInternalError("Attempt to deliver null data");
  if (!fn)
    return; // Caller doesn't need the data. This is ok.
  if (offset == 0) {
    _deliver_now(fn, data, size, transient);
  }
  else {
    auto dumb_ptr  = static_cast<const char*>(data.get()) + offset;
    auto smart_ptr = std::shared_ptr<const void>(data, dumb_ptr); // aliased ptr
    _deliver_now(fn, smart_ptr, size, transient);
  }
}

/**
 * Rudimentary pool of scratch buffers to avoid alloc/dealloc overhead.
 * Allocated data is returned as a std::shared_ptr that takes care of release.
 * Once allocated the the memory might not be released to the CRT until the
 * application exists. The function should only be used for short lived data.
 * Also, do NOT create any std::weak_ptr instances referencing the result.
 */
std::shared_ptr<void>
FileADT::_allocate(std::int64_t size)
{
  const std::int64_t minsize{16*1024}, maxsize{1024*1024};
  const std::int64_t highwater{500}, lowwater{200};
  static std::vector<std::shared_ptr<void>> cache;
  static std::mutex mutex;
  static size_t hint{0};
  std::shared_ptr<void> result;
  if (size >= minsize && size <= maxsize && malloc_shortcut() >= 1) {
    std::lock_guard<std::mutex> lk(mutex);
    if (hint >= cache.size())
      hint = 0;
    // Start searching for an available entry at the point we left off last.
    // This might *slightly* improve performance at the cost of messier code.
    // And be careful about the corner case where the cache is empty.
    std::vector<std::shared_ptr<void>>::const_iterator it, start = cache.begin() + hint;
    for (it = start; it != cache.end(); ++it) {
        if (it->use_count() == 1) {
            hint = it - cache.begin() + 1;
        return *it;
      }
    }
    for (it = cache.begin(); it != start; ++it) {
        if (it->use_count() == 1) {
            hint = it - cache.begin() + 1;
        return *it;
      }
    }
    result = std::shared_ptr<void>(::malloc(maxsize), [](void* p){::free(p);});
    if (result) {
      // This eviction method is crude but is not expected to be needed at all.
      // Assuming that all allocated data really is short lived, but for some
      // reason there is so much of it that keeping it allocated after
      // everything settles down is not practical. Whether free or in-use
      // pointers are evicted doesn't really matter because in-use data will
      // be freed very soon anyway. The reason for the shuffle is that if
      // long lived data is allocated by accident it will eventually get
      // evicted from the cache and becomes regular allocated data.
      // To force the eviction code to be executed, try to read a large file
      // usinh 1,000 or so threads.
      if ((std::int64_t)cache.size() > highwater && highwater > lowwater) {
        //std::cerr << "FileADT::_allocate() is evicting entries." << std::endl;
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(cache.begin(), cache.end(), g);        cache.resize(lowwater);
        hint = 0;
      }
      cache.push_back(result);
    }
  }
  else {
    // Bypass everything and just do a normal malloc.
    result = std::shared_ptr<void>(::malloc(size), [](void* p){::free(p);});
  }
  if (!result) {
    // This should be a std::bad_alloc but it might not be legal to throw that
    // exception from user code. Using new instead of malloc would make the
    // issue moot, but malloc seems safer both with respect to type punning
    // and alignmemt issues.
    throw std::runtime_error("Failed to allocate " + std::to_string(size) + " bytes");
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
//    FileFactory   /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Try the registered factories in the order of registration until one
 * is found that can handle this file. Caveat: If registration is done
 * using static initializers then that order is undefined. So when a
 * factory decides whether to handle a file or not it should not make
 * assumptions about where it is in the list.
 *
 * Thread safety: The factory is synchronized using a lock.
 * The lock is dropped before the actual file open.
 */
std::shared_ptr<IFileADT>
FileFactory::create(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext)
{
  // Need to copy the registry because the factory's lock must not be
  // held while opening the file. Should not be a performance issue
  // since the open itself is usually expensive and infrequent,
  std::vector<factory_t> registry_copy;
  {
    std::lock_guard<std::mutex> lk(_mutex);
    registry_copy = _registry;
  }
  std::shared_ptr<IFileADT> result;
  for (const factory_t& f : registry_copy) {
    result = f(filename, mode, iocontext);
    if (result)
      break;
  }
  if (!result)
    throw OpenZGY::Errors::ZgyUserError("Don't know how to open \"" + filename + "\".");
  return result;
}

/**
 * Register a new factory.
 *
 * It is allowed to do this from a static constructor. When this
 * method is entered we know that the _registry data member has been
 * initialized because the static _instance is declared inside
 * FileFactory::instance() which means the compiler is responsible for
 * constructing it early enough.
 *
 * Similarly, _mutex has already been constructed when any member
 * function is called. So it can be used to synchronize access.
 * The mutex cannot be use to protect the instance itself.
 * Fortunately it doesn't need to, because the compiler does that.
 *
 * Thread safety: The factory is synchronized using a lock.
 */
void
FileFactory::add_factory(const factory_t& factory)
{
  std::lock_guard<std::mutex> lk(_mutex);
  _registry.push_back(factory);
}

FileFactory&
FileFactory::instance()
{
  static FileFactory _instance;
  return _instance;
}

FileFactory::FileFactory()
  : _registry()
{
}

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> FileCommon   ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

FileCommon::FileCommon(const std::string& filename, OpenMode mode)
  : FileADT()
  , _mode(mode)
  , _name(filename)
  , _eof(0)
{
  _synctimer.reset(new SummaryPrintingTimerEx("File.sync"));
  _rtimer.reset(new SummaryPrintingTimerEx(mode == OpenMode::ReadWrite || mode == OpenMode::Truncate ? "File.reread" : "File.read"));
  _wtimer.reset(new SummaryPrintingTimerEx("File.write"));
  _mtimer.reset(new SummaryPrintingTimerEx("File.mutex"));
}

/**
 * \brief Get the current file size for error reporting.
 *
 * The default implementation in the base class just assumes that the
 * xx_eof() that is (probably) maintained internally is correct.
 *
 * Normally used just for error reporting, so it doesn't need to be performant.
 */
std::int64_t
FileCommon::_real_eof() const
{
  return xx_eof();
}

/**
 * \brief Throw a descriptive error if there was something wrong with the read.
 * \details Currently works for local files only. TODO-Low fix?
 * If fixing this then make sure all implementations of xx_eof()
 * and _real_eof() are also thread safe.
 */
void
FileCommon::_check_short_read(std::int64_t offset, std::int64_t size, std::int64_t got) const
{
  using OpenZGY::Errors::ZgyEndOfFile;
  using OpenZGY::Errors::ZgyInternalError;

  if (got == size)
    return;

  std::string msg = ("Cannot read offset " + _nice(offset) +
                     " size " + _nice(size) + ": ");
  if (got > size) {
    // Likely some kind of OS error. Beware possible buffer overrun.
    throw ZgyInternalError(msg + "got too much data: " + _nice(got) + ".");
  }
  else if (offset + size > xx_eof()) {
    // This can only happen if I (bug!) forgot to call _validate_read.
    throw ZgyEndOfFile(msg + "past EOF at " + _nice(xx_eof()) + ".");
  }
  else if (_real_eof() < xx_eof()) {
    // This can happen if opening /dev/null for read/write,
    // or if a write failed due to a full disk (and was not checked),
    // or I somehow (bug!) failed to keep track of eof while writing.
    // Or maybe some other process truncated the file.
    throw ZgyEndOfFile(msg + "File is shorter than expected: " +
                       _nice(_real_eof()) + " / " +
                       _nice(xx_eof()) + ".");
  }
  else {
    // The os returned a short read for no apparent reason.
    // Maybe the file is a special device other than /dev/null.
    throw ZgyEndOfFile(msg + "short read for unknown reason.");
  }
}

} // namespace
