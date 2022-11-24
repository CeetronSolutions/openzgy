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

#pragma once

/**
\file: file.h
\brief Low level I/O, abstract layer.

This file contains the base class for low level I/O either to on-prem data
using the regular read and write methods of the OS or to a cloud back-end.
\code{.unparsed}
    // TODO-Low: To improve isolation, user visible context such as
    // OpenZGY::SeismicStoreIOContext could be be copied into an equivalent
    // InternalZGY::SDConfig. Base class InternalZGY::Config in this header.
    // Ditto for plain files, an InternalZGY::FileConfig defined here.
    // The downside is that it gets more tedious to maintain.
    InternalZGY::Config:
      InternalZGY::FileConfig(Config):
      InternalZGY::SDConfig(Config):

        * Details such as user credentials etc. established when the
          file is open. Specific to the backend type.
        * Note that there is currently no way to pass a configuration
          object along with every read and write request. This might
          have been useful for a server type application but would
          require the config parameter to ripple across at least 50
          existing methods. I doubt this would be worth the trouble.
        * Higher level code should only access the polymorphic IFileADT
          base class and the InternalZGY::FileFactory that will create
          an instance of the desired type.

                            IFileBase (Pure virtual interface)
                                |
                                v
                            IFileADT (Pure virtual interface)
                                |
                                v
                            FileADT (Protected static methods only)
                             /        \
                            /          \
                           /            \
                          /              \
                         /          FileCommon[WinAndLinux]
                        /              /          \
                       /              /            \
         SeismicStoreFile    LocalFileLinux    LocalFileWindows

\endcode
*/

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <mutex>

#include "../declspec.h"

namespace OpenZGY {
  class IOContext;
}

namespace InternalZGY {
#if 0
}
#endif

class SummaryPrintingTimerEx;

enum class OpenMode
{
  Closed = 0,
  ReadOnly,
  ReadWrite,
  Truncate,
};

enum class UsageHint
{
  Unknown    = 0x00,
  TextFile   = 0x01,
  Header     = 0x10,
  Data       = 0x20,
  Compressed = 0x40,
  Mixed      = 0x40,
};

/**
 * Single entry in a scatter/gather read request.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class ReadRequest
{
public:
  typedef std::shared_ptr<const void> data_t;
  typedef std::function<void(data_t, std::int64_t)> delivery_t;
  std::int64_t offset;
  std::int64_t size;
  delivery_t delivery;
  ReadRequest(std::int64_t offset_in, std::int64_t size_in, const delivery_t& delivery_in)
    : offset(offset_in)
    , size(size_in)
    , delivery(delivery_in)
  {
  }
};

typedef std::vector<ReadRequest> ReadList;
typedef std::vector<ReadList> ReadDoubleList;

/*
 * \brief Low level I/O utilities, not retated to an open file.
 *
 * This class declares an interface for utility functions specific to one
 * backend. The actual implementaton is in some derived IFileADT class,
 * allowing us to use the same factory mechanism to access the functions.
 *
 */
class OPENZGY_TEST_API IFileBase
{
public:
  virtual ~IFileBase();
  virtual void deleteFile(const std::string& name, bool missing_ok) const = 0;
  virtual std::string altUrl(const std::string& name) const = 0;
  virtual std::string idToken() const = 0;
};

/**
 * \brief Internal interface for I/O operations.
 *
 * Public methods are prefixed with xx_ for practical reasons.
 * It makes it more obvious that an invocation is being made
 * on the IFileADT interface. It also becomes simple to search
 * for usage.
 *
 * Thread safety: Interfaces and classes that only contain static
 * methods do not have race conditions.
 */
class OPENZGY_TEST_API IFileADT : public IFileBase
{
public:
  virtual ~IFileADT();

  /**
   * Read binary data from the file. Both size and offset are mandatory.
   * I.e. caller is not allowed to read "the entire file", and not
   * allowed to "read from where I left off the last time".
   * The actual reading will be done in a derived class.
   * The base class only validates the arguments.
   */
  virtual void xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) = 0;

  /**
   * Read binary data from multiple regions in the file. Each part
   * of the request specifies offset, size, and a delivery functor
   * which will be invoked to pass back the returned bulk.
   *
   * Arguments:
   *     parallel_ok:  If true then the delivery functor might be called
   *                   simultaneously from multiple worker threads.
   *                   The function itself will block until all the data
   *                   has been read or an error occurs.
   *     immutable_ok: If true the caller promises that the delivery
   *                   functor will not try to modify the data buffer.
   *                   Pass False e.g. if the functor may need to byteswap
   *                   the data it has read from file.
   *                   With the current implementation the bulk layer
   *                   will uncondiionally pass false because it doesn't
   *                   know yet whether byeswap and/or subtiling is needed.
   *                   With the current implementation this doesn't add
   *                   much cost to the cloud reader so this is probably ok.
   *     transient_ok: If true the caller promises that the delivery
   *                   functor will not keep a reference to the data buffer
   *                   after the functor returns.
   *                   With smart pointers it is possible to check whether
   *                   the delivery functor kept its promise and signal
   *                   a fatal error if it didn't. The reason that the code
   *                   doesn't just allow keeping a pointer and look at the
   *                   refcount on return is that this might make a future
   *                   cache module less efficient.
   *
   * The delivery functor is called as
   *     fn(const std::shared_ptr<const void>& data, std::int64_t size)
   *
   * size can in some cases be more than originally requested due to
   * caching and possibly less if end of file was encountered.
   *
   * FUTURE: a new argument partial_ok may be set to True if it is ok to
   * call the delivery functor with less data than requested, and to keep
   * calling it until all data has been delivered. The signature of the
   * delivery functor gets changed to fn(data, offset, size). Offset is the
   * absolute file offset. I.e. not relative to the requested offset.
   * Passing partial_ok=True might elide some buffer copies if the
   * caller is doing something simple (such as reading an uncompressed
   * brick) where partial copies are possible, and the backend is in the
   * cloud, and a longer lived cache is being maintained, and the cache
   * block size is smaller than the requested size. That is a lot of ifs.
   * There was some code to handle partial_ok but it has been removed.
   * Get it from the git history if you really want it.
   */
  virtual void xx_readv(const ReadList& requests, bool parallel_ok=false, bool immutable_ok=false, bool transient_ok=false, UsageHint usagehint=UsageHint::Unknown) = 0;

  /**
   * Write binary data to the file. Offset is mandatory. I.e. caller
   * is not allowed to "write to where I left off the last time".
   * The actual writing will be done in a derived class.
   * The base class only validates the arguments.
   */
  virtual void xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) = 0;

  /**
   * Close the file. This should be done exactly once if the file was
   * opened normally, and not at all if the file was "opened" with mode
   * OpenMode::Closed. The latter is used when a backend handle is needed
   * e.g. for deleting files and the handle doesn't represent an open file.
   * Explicitly calling xx_close() more than once will raise an error.
   *
   * If the application forgets to close then the destructor will do so.
   * But in that case any exceptions will be caught and swallowed.
   *
   * Thread safety: Must not be called concurrently with any other method.
   */
  virtual void xx_close() = 0;

  /**
   * Return the current end-of-file, i.e. the file size.
   */
  virtual std::int64_t xx_eof() const = 0;

  /**
   * Return the size of each segment of the file, if the backend
   * has a concept of multiple segments in one file. Otherwise
   * just return xx_eof() in the first and only slot.
   *
   * If complete=false return at most 3 numbers: The first, second,
   * and last segment size. Currently all segments except the
   * first and last are required to have the same size, so by
   * combining the results of xx_segments() and xx_eof() it is
   * possible to compute the rest of the information.
   *
   * If the file is open for write then the last number will be
   * the in-memory buffer. That can be zero and it can also be
   * larger than the preferred segment size.
   *
   * Currently the only reason this is needed (apart from debug)
   * is a couple of consistency checks when opening a cloud file
   * for update. This is unfortunate because I really wanted to
   * keep this api as small as possible.
   */
  virtual std::vector<std::int64_t> xx_segments(bool complete) const = 0;

  /**
   * Return true if the file is on the cloud.
   * This might trigger some optimizations.
   */
  virtual bool xx_iscloud() const = 0;
};

/**
 * \brief Internal interface for I/O operations.
 *
 * The class adds contains some protected static convenience
 * methods that specializations might need. All instance methods
 * apart from the destructor remain pure virtual.
 *
 * Thread safety: Interfaces and classes that only contain static
 * methods do not have race conditions.
 */
class OPENZGY_TEST_API FileADT : public IFileADT
{
public:
  virtual ~FileADT();

protected:
  static std::string _nice(std::int64_t n);
  static void _validate_read(void *data, std::int64_t offset, std::int64_t size, std::int64_t eof, OpenMode mode);
  static void _validate_write(const void *data, std::int64_t offset, std::int64_t size, OpenMode mode);
  static void _validate_readv(const ReadList& requests, std::int64_t eof, OpenMode mode);
public: // Actually internal. Used by ConsolidateRequests.
  static void _deliver(const ReadRequest::delivery_t& fn, const ReadRequest::data_t&  data, std::int64_t offset, std::int64_t size, bool transient);
  static void _deliver_now(const ReadRequest::delivery_t& fn, const ReadRequest::data_t&  data, std::int64_t size, bool transient);
  static std::shared_ptr<void> _allocate(std::int64_t size);
};

/**
 * Factory for instanciating an appropriate concrete IFileADT instance.
 *
 * Thread safety: Synchronized using a lock.
 */
class OPENZGY_TEST_API FileFactory
{
public:
  typedef std::function<std::shared_ptr<IFileADT>(const std::string&, OpenMode, const OpenZGY::IOContext*)> factory_t;

public:
  std::shared_ptr<IFileADT> create(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  void add_factory(const factory_t& factory);
  static FileFactory& instance();

private:
  FileFactory();
  FileFactory(const FileFactory&) = delete;
  FileFactory& operator=(const FileFactory&) = delete;

private:
  std::vector<factory_t> _registry;
  std::mutex _mutex;
};

/**
 * \brief Implementation of some methods that might be shared.
 *
 * Using this class is optional. Concrete classes can inherit
 * directly from FileADT or even IFileADT if they want to.
 *
 * Thread safety: By design, all IFileADT specializations are expected to
 * allow concurrent reads but no guaranteed about anything else.
 * TODO-Worry: Need to analyze individual methods for thread safety issues.
 * This has been done. But maybe something slipped thru the cracks.
 */
class FileCommon : public FileADT
{
protected:
  OpenMode _mode;         // Only xx_close() allowed to change these,
  std::string _name;      // and that one is non threadsafe anyway.
  /**
   * Keep track of EOF: Initially set on file open; later kept
   * up to date if we write to the file.
   * Thread safety: Synchronized by the per-file mutex.
   */
  std::int64_t _eof;
  std::shared_ptr<SummaryPrintingTimerEx> _rtimer; // Access is thread safe
  std::shared_ptr<SummaryPrintingTimerEx> _wtimer; // Access is thread safe
  std::shared_ptr<SummaryPrintingTimerEx> _mtimer; // Access is thread safe

public:
  FileCommon(const std::string& filename, OpenMode mode);
  // NEW functions
  // TODO-Worry: All implementations need to be thread safe.
  virtual std::int64_t _real_eof() const;
  // TODO-Worry: calls xx_eof() and _real_eof() which need to be threadsafe.
  // Currently not true for SeismicStoreFileDelayedWrite, but this
  // check method won't be used there.
  void _check_short_read(std::int64_t offset, std::int64_t size, std::int64_t got) const;
};

}
