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

#include "api.h"
#include "iocontext.h"
#include "writerargs.h"
#include "exception.h"
#include "safewriter.h"
#include "impl/enum_mapper.h"
#include "impl/enum.h"
#include "impl/file.h"
#include "impl/meta.h"
#include "impl/bulk.h"
#include "impl/databuffer.h"
#include "impl/transform.h"
#include "impl/lodalgo.h"
#include "impl/statisticdata.h"
#include "impl/histogramdata.h"
#include "impl/histogrambuilder.h"
#include "impl/expandablebuilder.h"
#include "impl/genlod.h"
#include "impl/compression.h"
#include "impl/guid.h"
#include "impl/logger.h"
#include "impl/environment.h"
#include "impl/apirecorder.h"
#include "impl/writerargsimpl.h"

#include <tuple>
#include <list>
#include <sstream>
#include <fstream>
#include <atomic>
#include <mutex>
#include <assert.h>
#include <chrono>

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

/**
 * \file api.cpp
 * \brief Implements the pure interfaces of the API.
 *
 * \internal TODO-Doc consider moving all doxygen comments to the interface
 * and exclude the concrete types. Unfortunately than might make the comments
 * harder to keep up to date.
 */

namespace {
  /**
   * \brief Make a fake shared_ptr from a plain unsafe pointer.
   *
   * Used when a method or a data member expects a smart pointer while
   * the application has an unsafe pointer is guaranteed to remain
   * valid long enough. The returned smart pointer does nothing on
   * cleanup.
   *
   * CAVEAT: Using this method indicates a code smell.
   * Note: The method is duplicated in bulk.cpp, databuffer.cpp,
   * and api.cpp. There isn't much point in consolidatimg them.
   */
  template<typename T>
  static inline std::shared_ptr<T>
  fake_shared(T* p)
  {
    return std::shared_ptr<T>(p, [](T* p){});
  }
}

namespace {
  /**
   * Attach a logger. Used by the ZgyReader and ZgyWriter constructors.
   *
   * If the iocontext belongs to seismic store and if it has a logger
   * set then use that logger also for the higher level ZgyInteralMeta
   * and ZgyInternalBulk instances. And the other was around: If a
   * logger is created here, and a seismic store iocontext exists
   * without a logger, then use the iocontext to push the logger down
   * to the low level code.
   *
   * TODO-Medium: Bad code smell. The iocontext is supposed to be
   * populated by the application code only, and used by the file
   * layer only.
   */
  std::tuple<std::function<bool(int, const std::string&)>, std::shared_ptr<OpenZGY::IOContext>> setupLogging(const OpenZGY::IOContext *iocontext)
  {
    std::function<bool(int, const std::string&)> logger;
    std::shared_ptr<OpenZGY::IOContext> ctxt = iocontext ? iocontext->clone() : nullptr;
    auto iocontext_sd = dynamic_cast<OpenZGY::SeismicStoreIOContext*>(ctxt.get());
    if (iocontext_sd && iocontext_sd->getLogger()) {
      logger = iocontext_sd->getLogger();
    }
    else {
      logger = InternalZGY::LoggerBase::standardCallback
        (InternalZGY::LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"),
         "openzgy-api: ", "");
      // logger = std::ref(logger);
      if (iocontext_sd)
        iocontext_sd->logger(logger);
    }
    return std::make_tuple(logger, ctxt);
  }
}

namespace {
  class ZombieCheck {
  public:
    /**
     * The default is to silently sweep errors under the carpet,
     * after trying to prevent a crash due to code in the caller
     * that is not thread safe. E.g. closing a file in one thread
     * and reading from it in another.
     */
    static int complain_mode()
    {
      static int result = InternalZGY::Environment::getNumericEnv
        ("OPENZGY_COMPLAIN_IF_INVALID_MT", 0);
      return result;
    }

    /**
     * Report a serious problem with code in the caller that is not
     * thread safe. Or sweep it under the carpet. Throwing an exception
     * is also an option, but is discouraged. Because the error isn't
     * in this particular call. It just got caught here.
     */
    static void failure(const std::string& msg, const std::string& where)
    {
      if (complain_mode() & 0x01)
        std::cerr << (msg + " in " + where + "\n") << std::flush;
      if (complain_mode() & 0x02)
        assert(false && "OpenZGY ZombieCheck failed");
      if (complain_mode() & 0x04)
        throw OpenZGY::Errors::ZgyUserError(msg + " in " + where + "\n");
    }

    /**
     * Called when destructing a file, at which time there should be
     * no extra references to the internal accessor and file pointers.
     * If there are, this means there is probably an ongoing read
     * in a different thread. Possibly *this has beed deleted as well.
     * The latter triggers undefined behaviour. But if we are lucky
     * it won't actually crash.
     *
     * \internal The expected ptr.use_count() is 2, because the caller
     * has a reference and the "ptr" argument has another. That would
     * be the case even if ptr was declared as a const reference.
     * Presumably because the type won't match what the caller has.
     *
     * Keep in mind that _accessor also has a reference to _fd, so _fd
     * won't be unique until after _accessor has been deleted.
     */
    static void checkUniquePtr(std::shared_ptr<const void> ptr, const char *where)
    {
      if (ptr && ptr.use_count() > 2)
        failure("OpenZGY detected " +
                std::to_string(ptr.use_count() - 2) +
                " extra references", where);
    }

    /**
     * Called when reading or writing a file, at which time there
     * should be multiple references to the internal accessor
     * instance. One in the _accessor instance member and one in the
     * explicit local copy held by read() or write(). If there isn't
     * then *this has probably been destructed causing the _accessor
     * pointer to be released. A crash was probebly prevented by the
     * local copy.
     *
     * The check might have been done on _fd as well, but as long as
     * the accessor is alive it has its own reference to the file
     * instance.
     *
     * \internal when testing use_count, take into account that "ptr"
     * also holds a reference.
     */
    static void checkNonUniquePtr(std::shared_ptr<const void> ptr, const char *where)
    {
      // Expected at least one reference in the calling application,
      // one local reference in the function, and one for the "ptr"
      // argument to this function.
      if (ptr && ptr.use_count() < 3)
        failure("OpenZGY detected unexpected free", where);
    }
  };
} // namespace

/**
 * Parse a list of integers separated by one of comma,*,x,space,tab.
 * Returns an empty list if the environment variable is not found.
 * Throws if the string is not well formed.
 * Might want to move this inside class Environment.
 */
std::vector<int>
getNumericEnvList(const char* name, int expect_size)
{
  std::vector<int> result;
  std::string stl_str = InternalZGY::Environment::getStringEnv(name);
  if (stl_str.empty()) {
    return result;
  }
  else {
    const char* str = stl_str.c_str();
    for (;;) {
      char* end = nullptr;
      result.push_back(static_cast<int>(strtol(str, &end, 10)));
      if (end == str)
        throw std::runtime_error
          (std::string(name) +
           ": expected a number, found " +
           std::string("\"") + std::string(str) + "\".");
      else if (*end == '\0')
        break;
      else if (strchr(" \t,x*", *end) == nullptr)
        throw std::runtime_error
          (std::string(name) +
           ": numbers should be separated by one of space,tab,comma,x,*. Not " +
           std::string("'") + std::string(end).substr(0, 1) + "'.");
      str = end + 1;
    }
  }
  if (!result.empty() && expect_size >= 1 && (int)result.size() != expect_size)
    throw std::runtime_error
    (std::string(name) +
      ": expected " + std::to_string(expect_size) + " numbers. Not " +
      std::to_string(result.size()) + ".");
  return result;
}

namespace OpenZGY {
  class IOContext;
}

namespace OpenZGY { namespace Formatters {
    std::string enumToString(SampleDataType);
    std::string enumToString(UnitDimension);
    std::string enumToString(DecimationType);
    std::string enumToString(FinalizeAction);
    std::ostream& operator<<(std::ostream&, SampleDataType);
    std::ostream& operator<<(std::ostream&, UnitDimension);
    std::ostream& operator<<(std::ostream&, DecimationType);
    std::ostream& operator<<(std::ostream&, FinalizeAction);
  }
}

namespace OpenZGY { namespace Impl {
#if 0
}}
#endif

/**
 * Prevent the application from shooting itself in the foot by
 * opening a file more than once unless all opens are read-only.
 *
 * This will ONLY catch opens within a single process.
 * Also it is fairly simple to fool the check because e.g.
 * ./foo.zgy and ././foo.zgy will be considered different.
 *
 * The class can operate in complain-only mode by setting
 * OPENZGY_THROW_IF_FILE_LOCKED=0.
 *
 * A possible future option when in complain-only mode is to
 * assume there might be leaked handles open for read that
 * cannot be accessed anymore and thus are harmless. Output
 * an extra message if somebody tries to read from an assumed
 * leaked handle.
 */
class LockedInProcess
{
  /**
   * Describes the current lock. Note that the list of locks contains
   * a copy of LockedInProcess::entry_, not a reference to it.
   */
  struct Entry
  {
    std::string name_;
    bool writable_;
    std::int64_t seq_;
  };
  typedef std::list<Entry> list_t;
  typedef std::function<bool(int, const std::string&)> logger_t;

  Entry entry_;
  logger_t logger_;
  std::string myname_;

  static list_t locks_;
  static std::int64_t last_seq_;
  static std::mutex mutex_;

public:
  /**
   * Check for conflicting read and write locks.
   * Will normally throw if the application is about to shoot itself
   * in the foot, but can be changed in reportError to get the old
   * behavion of "the customer is always right".
   */
  LockedInProcess(const std::string& name, bool writable, const logger_t& logger)
    : entry_{name, writable, -1}
    , logger_(logger)
    , myname_()
  {
    std::lock_guard<std::mutex> lock(mutex_);
    entry_.seq_ = ++last_seq_;
    myname_ = "\"" + entry_.name_ + "\" (id " + std::to_string(entry_.seq_) + ")";
    checkLocks();
    locks_.push_back(entry_);
  }

  /**
   * Remove our lock.
   */
  ~LockedInProcess()
  {
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto it = locks_.begin(); it != locks_.end(); ++it) {
      if (it->seq_ == this->entry_.seq_) {
        reportOpenClose("Ready to close: " + myname_);
        locks_.erase(it);
        break;
      }
    }
  }

private:
  static bool throwOnError()
  {
    static bool result = InternalZGY::Environment::getNumericEnv
      ("OPENZGY_THROW_IF_FILE_LOCKED", 0) > 0;
    return result;
  }

  /**
   * List all file descriptors currently locking this file.
   * Caller must hold the mutex.
   */
  std::string lockedBy(const std::string& name) const
  {
    std::stringstream ss;
    ss << "\"" << name<< "\" ids";
    int count = 0;
    for (const auto& it : locks_) {
      if (it.name_ == name) {
        ss << " (" << it.seq_ << ":" << (it.writable_ ? "W" : "R") << ")";
        ++count;
      }
    }
    if (count == 0) {
      ss << " (none)";
    }
    return ss.str();
  }

  /**
   * Check for conflicting read and write locks.
   * Caller must hold the mutex.
   */
  void checkLocks() const
  {
    std::int64_t rcount{0};
    std::int64_t wcount{0};
    for (const auto& it : locks_) {
      if (it.name_ == entry_.name_ && it.seq_ != entry_.seq_) {
        if (it.writable_)
          ++wcount;
        else
          ++rcount;
      }
    }

    if (wcount > 0)
      reportError("Already opened for write: " + lockedBy(entry_.name_));
    else if (entry_.writable_ && rcount > 0)
      reportError("Already opened for read: " + lockedBy(entry_.name_));
    else if (entry_.writable_)
      reportOpenClose("Ready to open for write: " + myname_);
    else
      reportOpenClose("Ready to open for read: " + myname_);
  }

  /**
   * This is a good place to set a breakpoint.
   */
  void reportError(const std::string& msg) const
  {
    logger_(0, "WARNING: " + msg);
    if (throwOnError())
      throw Errors::ZgyUserError(msg);
  }

  /**
   * This is a good place to set a breakpoint.
   */
  void reportOpenClose(const std::string& msg) const
  {
    logger_(1, msg);
  }
};

LockedInProcess::list_t LockedInProcess::locks_;
std::int64_t LockedInProcess::last_seq_;
std::mutex LockedInProcess::mutex_;

/**
 * \brief High level API for reading and writing ZGY files.
 *
 * The base class contains properties common to both reader and writer.
 *
 * The constructor should logically have had a ZgyInternalMeta parameter
 * for accessing the implementation layer. But due to the way the code
 * is structured the _meta pointer needs to be set in the constructor
 * for the leaf types. The _meta pointer is guaranteed to not be empty
 * except while constructig the instance.
 *
 * Both ZgyReader and ZgyWriter need to retain a pointer to the file
 * descriptor. Primarily in order to explcitly close it. Relying on the
 * shared_ptr to do this when going out of scope is dangerous because
 * of exception handling. ZgyWriter additionally needs it when flushing
 * metadata to disk. Since there is no file descriptor in ZgyInternalMeta.
 *
 * Both ZgyReader and ZgyWriter need a ZgyInternalBulk pointer to do
 * bulk I/O. That instance in turn keeps private references to the
 * file descriptor and the metadata but is not supposed to expose them.
 *
 * The IFileADT and ZgyInternalBulk pointers can be declared here since
 * both the reader and the writer needs them. Or removed from here and
 * duplicated in ZgyReader and ZgyWriter.
 *
 * Thread safety:
 * The following applies both to this class and derived types.
 * Modification may lead to a data race. The intent for this class is
 * that any calls a ZgyReader instance might make will not modify data.
 * Except possibly changes done while the file is being opened.
 * To help enforce this, both const and mutable pointers are used for
 * ZgyInternalMeta and ZgyInternalBulk with the latter being empty
 * when instanciated from the ZgyReader constructor instead of ZgyWriter.
 * FUTURE: The instance holds one or more pointers to the data it
 */
class ZgyMeta: virtual public IZgyMeta
{
protected:
  /** \brief Handle to the internal metadata layer which this class wraps. */
  std::shared_ptr<const InternalZGY::ZgyInternalMeta> _meta;

  /** \brief For logging and debugging only. Identify this instance. */
  std::string _myname;

public:

  std::array<int64_t,3>
  size() const override
  {
    return _meta->ih().size();
  }

  SampleDataType
  datatype() const override
  {
    return EnumMapper::mapRawDataTypeToSampleDataType(_meta->ih().datatype());
  }

  std::array<float32_t,2>
  datarange() const override
  {
    if (_meta->ih().datatype() == InternalZGY::RawDataType::Float32)
      return std::array<float32_t,2>{_meta->ih().smin(), _meta->ih().smax()};
    else
      return _meta->ih().safe_codingrange();
  }

  std::array<float32_t,2>
  raw_datarange() const override
  {
    if (_meta->ih().datatype() == InternalZGY::RawDataType::Float32)
      return std::array<float32_t,2>{_meta->ih().smin(), _meta->ih().smax()};
    else
      return _meta->ih().raw_codingrange();
  }

  UnitDimension
  zunitdim() const override
  {
    return EnumMapper::mapRawVerticalDimensionToUnitDimension(_meta->ih().vdim());
  }

  UnitDimension
  hunitdim() const override
  {
    return EnumMapper::mapRawHorizontalDimensionToUnitDimension(_meta->ih().hdim());
  }

  std::string
  zunitname() const override
  {
    return _meta->ih().vunitname();
  }

  std::string
  hunitname() const override
  {
    return _meta->ih().hunitname();
  }

  float64_t
  zunitfactor() const override
  {
    return _meta->ih().vunitfactor();
  }

  float64_t
  hunitfactor() const override
  {
    return _meta->ih().hunitfactor();
  }

  float32_t
  zstart() const override
  {
    return _meta->ih().orig()[2];
  }

  float32_t
  zinc() const override
  {
    return _meta->ih().inc()[2];
  }

  std::array<float32_t,2>
  annotstart() const override
  {
    return std::array<float32_t,2>{
      _meta->ih().orig()[0],
      _meta->ih().orig()[1]};
  }

  std::array<float32_t,2>
  annotinc() const override
  {
    return std::array<float32_t,2>{
      _meta->ih().inc()[0],
      _meta->ih().inc()[1]};
  }

  const corners_t
  corners() const override
  {
    return _meta->ih().ocp_world();
  }

  const corners_t
  indexcorners() const override
  {
    return _meta->ih().ocp_index();
  }

  const corners_t
  annotcorners() const override
  {
    return _meta->ih().ocp_annot();
  }

  std::array<int64_t,3>
  bricksize() const override
  {
    // THIS MIGHT CONFUSE THE CALLING APPLICATON. USE ONLY FOR TESTING.
    // Mega-kludge. Lie about the physical brick size.
    // This (sort of) works because the API layer in OpenZGY is supposed to
    // completely hide the fact that the file is bricked. So, bricksize()
    // might be considered just as a hint about the optimal size to read.
    // Petrel and possibly others uses this. This trick is the only way to
    // get Petrel to read more than one physical brick at a time.
    // This could allow more multi-threading compared to just changing the
    // physical brick size.
    // Caveat: This also changes the brick size in BASE. Probably.
    // Caveat: The override also affects ZgtWriterArgs::metafrom().
    // Caveat: If bricksize() in the API layer gets called from the
    // implementation layer, this will be very bad.
    // Caveat: Utilities such as ZgyDump will also pick up this lie,
    // reporting the wrong brick size.
    // Changing OpenZgyBulkAccessor constructor in Salmon might have been better.
    static std::vector<int> force =
      getNumericEnvList("OPENZGY_HACK_LOGICAL_BRICKSIZE", 3);
    if (force.size() == 3) {
      static std::atomic_flag warned;
      if (!warned.test_and_set()) {
        std::stringstream ss;
        ss << "ZGY logical size forced to "
           << force[0] << "," << force[1] << "," << force[2];
        //this->_logger(0, ss.str());
        // _logger only in leaf classes. But this is a major hack anyway.
        std::cerr << (ss.str() + "\n");
      }
      return std::array<int64_t, 3>{force[0], force[1], force[2]};
    }
    return _meta->ih().bricksize();
  }

  std::vector<std::array<int64_t,3>>
  brickcount() const override
  {
    return _meta->ih().lodsizes();
  }

  int32_t
  nlods() const override
  {
    // nlods() is meant to be used from application code. Since ZgyWriter
    // has no methods with a "lod" parameter, there is no lowres accessible
    // to the application during write even when the internals have it.
    // See ZgyInternalBulk::validateUserPosition.
    return this->_meta->cached_nlods();
  }

  // Currently not needed by any client.
  //std::string dataid() const override
  //{
  //  return InternalZGY::GUID(_meta->ih().dataid()).toString();
  //}

  std::string verid() const override
  {
    return InternalZGY::GUID(_meta->ih().verid()).toString();
  }

  // Currently not needed by any client.
  //std::string previd() const override
  //{
  //  return InternalZGY::GUID(_meta->ih().previd()).toString();
  //}

  /**
   * Output general information about the ZgyReader or ZgyWriter.
   * This is primarily meant for debugging and testing, because
   * the application has no control over the output format.
   * Currently this method is not part of the api but it is
   * used internally by dump(). See also very similar methods
   * in tools/zgydumpc.cpp and test/test_api.cpp
   *
   * When used for test coverage, this should be accessing all
   * methods implemented in ZgyMeta, plus the statistics(),
   * histogram(), and filestats() methods, plus members of the
   * first two but mostly not members of FileStatistics.
   */
  static std::string toString(
       const IZgyMeta& meta,
       const SampleStatistics& stat,
       const SampleHistogram&  hist,
       const FileStatistics& info)
  {
    std::stringstream ss;
    using namespace OpenZGY;
    using namespace OpenZGY::Formatters;

    ss << "File format and version        = "
       << meta.datatype() << " ZGY version " << info.fileVersion() << "\n";
    //<< "Data identifier                = " << meta.dataid() << "\n";
    ss << "Current data Version           = "
       << meta.verid() << "\n";
    //<< "Previous data version          = " << meta.previd() << "\n";
    ss << "Size I,J,K                     = "
       << "(" << meta.size()[0] << ", "
       << meta.size()[1] << ", " <<
      meta.size()[2] << ")\n";
    ss << "Brick size I,J,K               = "
       << "(" << meta.bricksize()[0] << ", "
       << meta.bricksize()[1] << ", " <<
      meta.bricksize()[2] << ")\n";
    ss << "Number of bricks I,J,K         = "
       << "(" << meta.brickcount()[0][0]
       << ", " << meta.brickcount()[0][1]
       << ", " << meta.brickcount()[0][2]
       << ")\n";
    ss << "Number of LODs                 = "
       << meta.nlods() << "\n";
    ss << "Coding range min/max           = "
       << std::setprecision(6)
       << meta.datarange()[0] << " " << meta.datarange()[1]
       << " (raw: "
       << meta.raw_datarange()[0] << " " << meta.raw_datarange()[1]
       << ") "
       << meta.size()[0] * meta.size()[1] * meta.size()[2] << "\n";
    ss << "Statistical min/max/count/avg  = "
       << std::setprecision(6)
       << stat.min << " " << stat.max << " " << stat.cnt
       << " " << stat.sum / stat.cnt
       << " " << std::sqrt(stat.ssq / stat.cnt)
       << "\n";
    ss << "Histogram range min/max/count  = "
       << std::setprecision(6)
       << hist.minvalue << " " << hist.maxvalue << " " << hist.samplecount
       << " bincount " << hist.bins.size() << "\n";
    ss << "Inline start/increment/count   = "
       << meta.annotstart()[0]
       << " " << meta.annotinc()[0]
       << " " << meta.size()[0] << "\n";
    ss << "Xline  start/increment/count   = "
       << meta.annotstart()[1]
       << " " << meta.annotinc()[1]
       << " " << meta.size()[1] << "\n";
    ss << "Sample start/increment/count   = "
       << meta.zstart()
       << " " << meta.zinc()
       << " " << meta.size()[2] << "\n";
    //ss << "Horizontal projection system   = "
    //   << "?\n" // {r._accessor._metadata._ih._hprjsys};
    ss << "Horizontal dim/factor/name     = "
       << meta.hunitdim()
       << " " << meta.hunitfactor()
       << " '" << meta.hunitname() << "'\n";
    ss << "Vertical dim/factor/name       = "
       << meta.zunitdim()
       << " " << meta.zunitfactor()
       << " '" << meta.zunitname() << "'\n";
    ss << "Ordered Corner Points Legend   = "
       << "[  <i>,   <j>] { <inline>,   <xline>} (  <easting>,  <northing>)\n";
    for (int ix = 0; ix < 4; ++ix) {
      std::stringstream tt;
      tt << "Ordered Corner Point " << ix << "         = "
         << "["
         << std::setw(5) << meta.indexcorners()[ix][0] << ", "
         << std::setw(5) << meta.indexcorners()[ix][1] << "] {"
         << std::setw(9) << meta.annotcorners()[ix][0] << ", "
         << std::setw(9) << meta.annotcorners()[ix][1] << "} ("
         << std::fixed << std::setprecision(2)
         << std::setw(11) << meta.corners()[ix][0] << ", "
         << std::setw(11) << meta.corners()[ix][1] << ")"
         << "\n";
      ss << tt.str();
    }
    return ss.str();
  }

  void
  dump(std::ostream& os) const override
  {
    const SampleStatistics stat = statistics();
    const SampleHistogram  hist = histogram();
    std::shared_ptr<const FileStatistics> info = filestats();
    os << toString(*this, stat, hist, *info) << std::flush;
  }

  SampleStatistics
  statistics() const override
  {
    const InternalZGY::IInfoHeaderAccess& ih = _meta->ih();
    return SampleStatistics(ih.scnt(), ih.ssum(), ih.sssq(), ih.smin(), ih.smax());
  }

  SampleHistogram
  histogram() const override
  {
    const InternalZGY::IHistHeaderAccess& hh = _meta->hh();
    if (hh.bincount() > 0 && hh.bins() != nullptr)
      return SampleHistogram(hh.samplecount(), hh.minvalue(), hh.maxvalue(),
                             std::vector<std::int64_t>(hh.bins(),
                                                       hh.bins() + hh.bincount()));
    else
      return SampleHistogram();
  }

  /**
   * For debugging the code and for forensic analysis of corrupt files.
   * In OpenZGY this information is deliberately not avaiable in the
   * public API to avoid writing code that knows too much about internals.
   * TODO-Low: consider a way for "zgydumpc" to access this function.
   */
  void dumplut(
       std::ostream& os,
       std::int64_t ii,
       std::int64_t jj,
       std::int64_t kk,
       int lod) const
  {
    using InternalZGY::BrickStatus;
    InternalZGY::LookupTable::LutInfo info(BrickStatus::Missing,0,0,0);
    if (kk >= 0) {
      info = InternalZGY::LookupTable::getBrickFilePosition
        (ii, jj, kk, lod,
         _meta->ih().lodsizes(),
         _meta->ih().brickoffsets(),
         _meta->blup().lup(),
         _meta->blup().lupend(),
         _meta->ih().bytesperbrick());
      os << "[" << lod << "]"
         << "[" << ii << "]"
         << "[" << jj << "]"
         << "[" << kk << "] ";
    }
    else {
      info = InternalZGY::LookupTable::getAlphaFilePosition
        (ii, jj, lod,
         _meta->ih().lodsizes(),
         _meta->ih().alphaoffsets(),
         _meta->alup().lup(),
         _meta->ih().bytesperalpha());
      os << "[" << lod << "]"
         << "[" << ii << "]"
         << "[" << jj << "]"
         << "[" << "*" << "] ";
    }
    switch (info.status) {
    case BrickStatus::Missing:
      os << "missing";
      break;
    case BrickStatus::Constant:
      os << "constant";
      // Requires private function ZgyInternalBulk::decodeConstant,
      // arguably that functionality should have been in LookupTable.
      //os << " " << InternalZGY::ZgyInternalBulk::decodeConstant
      //  (info.raw_constant, _meta->ih().datatype());
      break;
    case BrickStatus::Normal:
      os << info.offset_in_file << " " << info.size_in_file;
      break;
    case BrickStatus::Compressed:
      os << info.offset_in_file << " " << info.size_in_file
         << " (compressed)";
      break;
    default:
      os << "ERROR";
      break;
    }
    if (kk < 0)
      os << " (alpha)";
    os << "\n";
  }

  void dumpluts(std::ostream& os) const
  {
    for (int lod = 0; lod < (int)_meta->ih().lodsizes().size(); ++lod) {
      std::array<std::int64_t,3> size = _meta->ih().lodsizes()[lod];
      for (std::int64_t ii = 0; ii < size[0]; ++ii) {
        for (std::int64_t jj = 0; jj < size[1]; ++jj) {
          for (std::int64_t kk = 0; kk < size[2]; ++kk) {
            dumplut(os, ii, jj, kk, lod);
          }
        }
      }
    }
    for (int lod = 0; lod < (int)_meta->ih().lodsizes().size(); ++lod) {
      std::array<std::int64_t,3> size = _meta->ih().lodsizes()[lod];
      for (std::int64_t ii = 0; ii < size[0]; ++ii) {
        for (std::int64_t jj = 0; jj < size[1]; ++jj) {
          dumplut(os, ii, jj, -1, lod);
        }
      }
    }
  }

  std::shared_ptr<const FileStatistics> filestats_nocache() const
  {
    using InternalZGY::LookupTable;
    using InternalZGY::BrickStatus;

    const std::int64_t bytesperalpha = _meta->ih().bytesperalpha();
    const std::int64_t bytesperbrick = _meta->ih().bytesperbrick();
    const std::vector<std::uint64_t>& alup = _meta->alup().lup();
    const std::vector<std::uint64_t>& blup = _meta->blup().lup();
    const std::vector<std::uint64_t>& bend = _meta->blup().lupend();

    FileStatistics result;
    result._data_start = std::numeric_limits<std::int64_t>::max();
    result._file_version = _meta->fh().version();
    result._alpha_normal_size_per_entry = bytesperalpha;
    result._brick_normal_size_per_entry = bytesperbrick;
    // result._file_size = _fd->xx_size(); Available in ZgyReader and ZgyWriter only.

    // TODO-Low: Fix this kludge.
    // I happen to know that in V3 and V4 the headers are all stored
    // consecutively and the brick lookup table comes last.
    result._header_size = _meta->oh().bricklupoff() + _meta->oh().bricklupsize();

    for (std::size_t ix = 0; ix < alup.size(); ++ix) {
      LookupTable::LutInfo info =
        LookupTable::getAlphaFilePositionFromIndex(ix, alup, bytesperalpha);
      switch (info.status) {
      case BrickStatus::Missing:  result._alpha_missing_count  += 1; break;
      case BrickStatus::Constant: result._alpha_constant_count += 1; break;
      case BrickStatus::Normal:
        result._alpha_normal_count   += 1;
        result._data_start = std::min(result._data_start, info.offset_in_file);
        break;
      case BrickStatus::Compressed:
        result._alpha_compressed_count += 1;
        result._alpha_compressed_size += info.size_in_file;
        result._data_start = std::min(result._data_start, info.offset_in_file);
        break;
      }
    }
    for (std::size_t ix = 0; ix < blup.size(); ++ix) {
      LookupTable::LutInfo info =
        LookupTable::getBrickFilePositionFromIndex(ix, blup, bend, bytesperbrick);
      switch (info.status) {
      case BrickStatus::Missing:  result._brick_missing_count  += 1; break;
      case BrickStatus::Constant: result._brick_constant_count += 1; break;
      case BrickStatus::Normal:
        result._brick_normal_count   += 1;
        result._data_start = std::min(result._data_start, info.offset_in_file);
        break;
      case BrickStatus::Compressed:
        result._brick_compressed_count += 1;
        result._brick_compressed_size += info.size_in_file;
        result._data_start = std::min(result._data_start, info.offset_in_file);
        break;
      }
    }

    // TODO-Low: Keep track of wasted_size and padding_size.
    // Padding gets added in ZgyInternalMeta::flushMeta(). I need to
    // replicate the logic here. The alternative is to scan for the
    // lowest brick offset. But even that isn't completely reliable
    // because there might be wasted blocks between end of padding and
    // start of first block. And, do I really care at all?
    //result._padding_size = roundup(result._header_size,
    //                              result._brick_size_per_entry);
    //result._wasted_size = result._file_size - result._usedSize();

    // DERIVED INFORMATION:
    // The following could also have been generated on the fly in some
    // member function. I pre-calculate it here instead, to limit the
    // amount of code visible in the public api.h header file.

    // File size not including padding and holes.
    result._used_size =
      ((result._alpha_normal_count * result._alpha_normal_size_per_entry) + result._alpha_compressed_size +
       (result._brick_normal_count * result._brick_normal_size_per_entry) + result._brick_compressed_size +
       result._header_size);

    // As used_size if the file is/was uncompressed.
    result._used_if_uncompressed =
      (((result._alpha_normal_count + result._alpha_compressed_count) * result._alpha_normal_size_per_entry) +
       ((result._brick_normal_count + result._brick_compressed_count) * result._brick_normal_size_per_entry) +
       result._header_size);

    // Is there at least one brick flagged as compressed?
    result._is_compressed =
      (result._alpha_compressed_count + result._brick_compressed_count > 0);

    // Relative size of this possibly compressed file compared to uncompressed.
    result._compression_factor =
      (result._used_if_uncompressed > 0 ?
       result._used_size / (double)result._used_if_uncompressed :
       1);
    // Slightly different definition of compression factor.
    // Doesn't work because file_size not set yet.
    // Besides, I like the other one better.
    //result._compression_factor =
    //  (result._file_size > 0 ?
    //   result._used_size + (result._file_size-(result._used_if_uncompressed)) / (double)result._file_size:
    //   1);

    return std::shared_ptr<const FileStatistics>(new FileStatistics(result));
  }

  /**
   * This method should be overwridden in ZgyReader and ZgyWriter.
   * If it is ok to have ZgyMeta be an abstract type then this
   * implementation can simply be removed.
   *
   * Calling this generic version of filestats() method will not
   * populate the file size and will not do any caching.
   */
  std::shared_ptr<const FileStatistics> filestats() const override
  {
    return filestats_nocache();
  }

  std::string myname() const { return _myname; }

protected:
  /** Only call this once, from constructors, to init _myname */
  static std::string makemyname(const std::string& filename)
  {
    static std::atomic<int> _debug_next_sequence(101);
    int seq = _debug_next_sequence.fetch_add(1);
    std::stringstream namess;
    namess << "\"" << filename << "\" (seq " << seq << ")";
    return namess.str();
  }
};

/**
 * \brief Add coordinate conversion to the concrete ZgyMeta class.
 *
 * \details Thread safety: See the base class. This specialization
 * does not add any additional concerns because all the new members
 * are const. So they are safe to call concurrently with other read
 * operations.
 */
class ZgyMetaAndTools: public ZgyMeta, virtual public IZgyTools
{
public:
  void
  transform(const corners_t& A, const corners_t& B, std::vector<std::array<float64_t,2>>& data) const override
  {
    // If this method is needed then it is fairly simple to implement.
    // Change InternalZGY::generalTransform() to accept a vector of array.
    // transform1() can then be re-implemented in terms of transform.
    throw std::runtime_error("Not implemented: ZgyMetaAndToole::transform()");
  }

  std::array<float64_t,2>
  transform1(const corners_t& A, const corners_t& B, const std::array<float64_t,2>& point) const override
  {
    float64_t x{point[0]}, y{point[1]};
    if (!InternalZGY::generalTransform
        (A[0][0], A[0][1], A[1][0], A[1][1], A[2][0], A[2][1],
         B[0][0], B[0][1], B[1][0], B[1][1], B[2][0], B[2][1],
         &x, &y, 1))
      throw Errors::ZgyInternalError("Transform is not well defined due to colinear or coincident control points");
    return std::array<float64_t,2>{x, y};
  }

  std::array<float64_t,2>
  annotToIndex(const std::array<float64_t,2>& point) const override
  {
    return transform1(annotcorners(), indexcorners(), point);
  }

  std::array<float64_t,2>
  annotToWorld(const std::array<float64_t,2>& point) const override
  {
    return transform1(annotcorners(), corners(), point);
  }

  std::array<float64_t,2>
  indexToAnnot(const std::array<float64_t,2>& point) const override
  {
    return transform1(indexcorners(), annotcorners(), point);
  }

  std::array<float64_t,2>
  indexToWorld(const std::array<float64_t,2>& point) const override
  {
    return transform1(indexcorners(), corners(), point);
  }

  std::array<float64_t,2>
  worldToAnnot(const std::array<float64_t,2>& point) const override
  {
    return transform1(corners(), annotcorners(), point);
  }

  std::array<float64_t,2>
  worldToIndex(const std::array<float64_t,2>& point) const override
  {
    return transform1(corners(), indexcorners(), point);
  }
};

/**
 * \brief Concrete implementation of IZgyReader.
 * \details Thread safety: See the base class.
 */
class ZgyReader : public ZgyMetaAndTools, virtual public IZgyReader
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

private:
  std::shared_ptr<InternalZGY::IFileADT> _fd;
  std::shared_ptr<const InternalZGY::ZgyInternalBulk> _accessor;
  std::shared_ptr<const FileStatistics> _filestats;
  mutable std::mutex _filestats_mutex;
  std::shared_ptr<LockedInProcess> _locked_in_process;
  LoggerFn _logger;

public:
  /**
   * \copydoc IZgyReader::open()
   */
  ZgyReader(const std::string& filename, const IOContext* iocontext)
  {
    _myname = makemyname(filename);
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(iocontext);
    _logger(1, "Open R " + myname());
    _locked_in_process.reset(new LockedInProcess(filename, false, _logger));

    // Note: Currently the Python version of this constructor has an
    // additional "update" parameter which, if set, causes the underlying
    // file to be opened in read/write mode. This is a kludge to facilitate
    // a stand alone lowres generator. Which was only ever used for testing.
    _fd = InternalZGY::FileFactory::instance().create(
        filename, InternalZGY::OpenMode::ReadOnly, ctxt.get());

    // Set the protected metadata information in the ZgyMeta base class.
    // ZgyInternalMeta does not retain the file descriptor. The file is
    // only used inside the constructor to populate the headers.
    _meta.reset(new InternalZGY::ZgyInternalMeta(this->_fd, _logger));

    // At the implementation level the bulk- and meta access are separate,
    // and the bulk accessor needs some of the meta information to work.
    // The accessor will have its own metadata pointer which it holds on to.
    // It also holds on to the file descriptor. Unlike the metadata reader
    // it is obviously not possible to read everything up front.
    // Note: the ZgyInternalBulk constructor will mostly just copy the
    // arguments to the instance. The exception is setting up tracking
    // of dirty bricks, statistics, and histogram which is N/A for reading.
    // Note: compressed_write is obviously false, as we aren't writing.
    _accessor.reset(new InternalZGY::ZgyInternalBulk
                    (_fd, _meta, nullptr, /*compressed_write=*/false,
                     std::array<float,2>{0,0}, _logger));

    if (_logger(4, "")) {
      std::stringstream ss;
      dumpluts(ss);
      _logger(4, ss.str());
    }
  }

  ~ZgyReader()
  {
    try {
      _logger(1, "Destroy R " + myname() + std::string((_accessor || _fd) ? "" : " (no-op)"));
      //close(); // See close() for why this may be a bad idea.
      ZombieCheck::checkUniquePtr(_accessor, "~ZgyReader bulk");
      _accessor.reset();
      ZombieCheck::checkUniquePtr(_fd, "~ZgyReader file");
      _fd.reset();
      // Debatable, but since nobody will be able to read from this
      // instance it shouldn't be a problem to remove the read lock.
      // If the checks for unique pointers failed, we might be in trouble.
      _locked_in_process.reset();
    }
    catch (const std::exception& ex) {
      // We should never get here!
      // Caller should have done an explicit close() so it can handle
      // exceptions itself. Exceptions thrown from a destructors are evil.
      // Trying to catch exceptions from the two lines above might already
      // be too late. The destructor in _fd does a similar operation
      // (blind catch with logging) which makes it even less likely
      // thet we get here.
      _logger(1, "ERROR closing a file opened for read: " +
              std::string(ex.what() ? ex.what() : "(null)"));
    }
  }

  /**
   * \brief Read an arbitrary region.
   *
   * \details \callgraph \callergraph
   * The data is read into a buffer provided by the caller.
   * The method will apply conversion storage -> float if needed.
   *
   * Data is ordered inline(slowest), crossline, vertical(fastest).
   *
   * The start position refers to the specified lod level.
   * At lod 0 start + data.size can be up to the survey size.
   * At lod 1 the maximum is just half that, rounded up.

   * It is valid to pass a size that includes the padding area
   * between the survey and the end of the current brick. But not
   * more. In other words, the limit for lod 0 is actually
   * reader()->size() rounded up to a multiple of reader->bricksize().
   */
  void read(const size3i_t& start, const size3i_t& size, float* data, int lod) const override
  {
    // This comment applies to all the read and write overloads.
    // Keep in mind that a file must not be closed while being read
    // in another thread. Keeping a local reference to the accessor's
    // implementation does NOT remove this rule. It can, however,
    // make it less likely that a buggy application will crash
    // accessing freed memory. There is STILL A WINDOW where this
    // might happen. But it is hopefully smaller.
    auto accessor = _accessor;
    throw_if_not_readable();
    if (!accessor->expeditedRead(start, size, data, lod, InternalZGY::RawDataType::Float32)) {
      std::shared_ptr<float> fakeshared = fake_shared(data);
      auto databuffer = std::make_shared<InternalZGY::DataBufferNd<float,3>>(fakeshared, size);
      accessor->readToExistingBuffer(databuffer, start, lod, true);
      databuffer.reset();
      if (fakeshared.use_count() != 1)
        throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    }
    ZombieCheck::checkNonUniquePtr(accessor, "read accessor");
  }

  /**
   * \brief Read an arbitrary region with no conversion.
   *
   * As the read overload with a float buffer but only works for files with
   * SampleDataType::int16 and does not scale the samples.
   */
  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod) const override
  {
    auto accessor = _accessor;
    throw_if_not_readable();
    if (!accessor->expeditedRead(start, size, data, lod, InternalZGY::RawDataType::SignedInt16)) {
      std::shared_ptr<std::int16_t> fakeshared = fake_shared(data);
      auto databuffer = std::make_shared<InternalZGY::DataBufferNd<std::int16_t,3>>(fakeshared, size);
      accessor->readToExistingBuffer(databuffer, start, lod, false);
      databuffer.reset();
      if (fakeshared.use_count() != 1)
        throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    }
    ZombieCheck::checkNonUniquePtr(accessor, "read accessor");
  }

  /**
   * \brief Read an arbitrary region with no conversion.
   *
   * As the read overload with a float buffer but only works for files with
   * SampleDataType::int8 and does not scale the samples.
   */
  void read(const size3i_t& start, const size3i_t& size, std::int8_t* data, int lod) const override
  {
    auto accessor = _accessor;
    throw_if_not_readable();
    if (!accessor->expeditedRead(start, size, data, lod, InternalZGY::RawDataType::SignedInt8)) {
      std::shared_ptr<std::int8_t> fakeshared = fake_shared(data);
      auto databuffer = std::make_shared<InternalZGY::DataBufferNd<std::int8_t,3>>(fakeshared, size);
      accessor->readToExistingBuffer(databuffer, start, lod, false);
      databuffer.reset();
      if (fakeshared.use_count() != 1)
        throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    }
    ZombieCheck::checkNonUniquePtr(accessor, "read accessor");
  }

  /**
   * \brief Get hint about all constant region.
   *
   * \details \callgraph \callergraph
   * Check to see if the specified region is known to have all samples
   * set to the same value. Returns a pair of (is_const, const_value).
   * If it returns is_const=true you know what the region contains;
   * otherwise you need to use the regular read() method.
   *
   * For int8 and int16 files the caller may specify whether to scale
   * the values or not. Even if unscaled the function returns the value
   * as a double.
   *
   * The function only makes inexpensive checks so it might return
   * is_const=false even if the region was in fact constant. It will not
   * make the opposite mistake. This method is only intended as a hint
   * to improve performance. The following figure might clarify how it works.
   * \image html readconst-fig1.png
   * \image latex readconst-fig1.png
   *
   * If the region requested in readconst() is exactly one ZGY brick
   * then this just checks whether that brick is known to be constant.
   * If called with a larger, arbitrary region it fails unless all
   * bricks in that region are known to be constant and all have the
   * same value. It will not tell you which brick(s) were not const.
   * If you need to know then you must loop over each ZGY brick and
   * make multiple calls to readconst().
   */
  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, int lod, bool as_float) const override
  {
    auto accessor = _accessor;
    throw_if_not_readable();
    auto result = accessor->readConstantValue(start, size, lod, as_float);
    ZombieCheck::checkNonUniquePtr(accessor, "read accessor");
    return result;
  }

  /**
   * \brief Close the file and release resources.
   *
   * Unlike ZgyWriter::close(), forgetting to close a file that was
   * only open for read is not a major faux pas. It can still lead to
   * problems though.
   *
   *   \li The destructor of _fd will catch and ignore any exceptions
   *       because if a destructor throws then this will normally
   *       cause the application to crash.
   *
   *   \li If a callback is used to refresh the token this will not
   *       happen in our destructor (it won't call xx_close) because
   *       it is too risky to invoke the callback this late. It might
   *       not be valid any longer. This means that if the token has
   *       expired since the last read then the close will fail.
   *       Exactly why SDAPI requires a token just to close a file is
   *       a different question. Possibly this is for removing any
   *       read locks.
   */
  void close() override
  {
    _logger(1, "Close R " + myname() + std::string((_accessor || _fd) ? "" : " (no-op)"));
    auto accessor = _accessor;
    if (accessor) {
      _accessor.reset();
      ZombieCheck::checkUniquePtr(accessor, "ZgyReader close");
      accessor.reset();
    }
    auto victim = _fd;
    if (victim) {
      _fd.reset();
      victim->xx_close();
      ZombieCheck::checkUniquePtr(victim, "ZgyReader close file");
      victim.reset();
    }
    // Maybe do this also if the ZombieCheck threw an exception?
    // Not a big deal, because it is discouraged to configure
    // ZombieCheck that way.
    _locked_in_process.reset();
    // Metadata remains accessible. Not sure whether this is a good idea.
  }

  /**
   * Get the file statistics of a file currently opened for read.
   * The result does not change while the file is open, so it can
   * be cached here.
   */
  std::shared_ptr<const FileStatistics> filestats() const override
  {
    if (!_filestats) {
      std::shared_ptr<FileStatistics> result
        (new FileStatistics(*filestats_nocache()));
      // The base class has no _fd member so I need to set the size here.
      result->_file_size = _fd->xx_eof();
      result->_segment_sizes = _fd->xx_segments(false);
      result->_data_start = std::min(result->_data_start, result->_file_size);
      {
        // Too bad there is no proper atomic_shared_ptr yet.
        std::lock_guard<std::mutex> lk(_filestats_mutex);
        // The cache is semantically const.
        const_cast<ZgyReader*>(this)->_filestats = result;
      }
    }
    std::lock_guard<std::mutex> lk(_filestats_mutex);
    return _filestats;
  }

private:
  /**
   * Test for common user errors as early as possible.
   */
  void throw_if_not_readable() const
  {
    if (!_fd || !_accessor || !_meta)
      throw Errors::ZgyUserError("ZGY file not open for read");
  }
};

/**
 * \brief Concrete implementation of IZgyWriter.
 *
 * \details Thread safety: This class is single threaded. IZgyWriter::open()
 * may choose to return a wrapper that serializes all access, just in case
 * the application didn't read the memo about no concurrent access.
 * See class ZgySafeWriter for details.
 *
 * Const-correctness: Most methods are non-const because writes will update
 * metadata. If not in ZgyWriter itself then in the classes it refers to.
 * The class refers to a mutable ZgyInternalBulk instance. And a mutable
 * ZgyInternalMeta instance that is actually the same as the one in the
 * base class except this one is not declared const. Now for the minor
 * code smell: Any const method (currently just errorflag) will still be
 * able to access those _rw pointers so the "const" declaration hardly
 * protects against anything. There are ways of handling this better.
 * But with just a single method affected I'll let it slide.
 */
class ZgyWriter : public ZgyMetaAndTools, virtual public IZgyWriter
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

private:
  std::shared_ptr<InternalZGY::IFileADT> _fd;
  std::shared_ptr<InternalZGY::ZgyInternalBulk> _accessor_rw;
  std::shared_ptr<InternalZGY::ZgyInternalMeta> _meta_rw;
  compressor_t _compressor;
  compressor_t _lodcompressor;
  std::shared_ptr<LockedInProcess> _locked_in_process;
  LoggerFn _logger;

public:
  struct OpenForUpdate{};

  /**
   * \copydoc IZgyWriter::open()
   */
  explicit ZgyWriter(const ZgyWriterArgsV3& args)
    : _fd()
    , _accessor_rw()
    , _compressor(args.impl().compressor_)
    , _lodcompressor(args.impl().lodcompressor_ ?
                     args.impl().lodcompressor_ :
                     args.impl().compressor_)
  {
    // Mutable copy needed only for the bs_override kludge, below.
    InternalZGY::ZgyWriterArgsV3Impl iargs(args.impl());
    _myname = makemyname(iargs.filename_);
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(iargs.iocontext_.get());
    _logger(1, "Open W " + myname());
    _locked_in_process.reset(new LockedInProcess(iargs.filename_, true, _logger));

    for (int ii=0; ii<3; ++ii) {
      if (iargs.size_[ii] < 1)
        throw Errors::ZgyUserError("Survey size cannot be empty or negative.");
      else if (iargs.size_[ii] >= std::numeric_limits<std::int32_t>::max())
        throw Errors::ZgyUserError("Survey size is too large.");
      }
    const bool compress = iargs.compressor_ || iargs.lodcompressor_;
    // This is both pointless and expensive to support.
    if (compress && iargs.datatype_ != InternalZGY::RawDataType::Float32)
      throw Errors::ZgyUserError("Compressed files need to be stored as float.");

    // Set the protected metadata information in the ZgyMeta base class.
    // Also store a mutable copy in the current class.
    // ZgyInternalMeta does not contain any file descriptor. This means
    // we need to hold on to _fd ourselves and provide it when it is time
    // to flush metadata to disk.

    // Mega-kludge. Allow overriding the brick size requested by the caller.
    // THIS MIGHT CONFUSE THE CALLING APPLICATON. USE ONLY FOR TESTING.
    std::vector<int> bs_override = getNumericEnvList("OPENZGY_HACK_OVERRIDE_BRICKSIZE", 3);
    if (bs_override.size() == 3) {
      iargs.bricksize_ = std::array<std::int64_t, 3>
        {bs_override[0], bs_override[1], bs_override[2]};
      iargs.have_bricksize_ = true; // Don't override later.
      std::stringstream ss;
      ss << "ZGY physical brick size forced to "
         << bs_override[0] << "," << bs_override[1] << "," << bs_override[2];
      _logger(0, ss.str());
    }
    _meta_rw.reset(new InternalZGY::ZgyInternalMeta(iargs, compress, _logger));
    _meta = _meta_rw; // in base class

    // The file creation was deferred until after the consistency checks.
    _fd = InternalZGY::FileFactory::instance().create(
        iargs.filename_, InternalZGY::OpenMode::Truncate, ctxt.get());
    _meta_rw->flushMeta(_fd);
    // Compress or not at this level only controls alignment etc.
    // The actual compression functor is passed to each write.
    // Note: the ZgyInternalBulk constructor will mostly just copy the
    // arguments to the instance. The exception is setting up tracking
    // of dirty bricks, statistics, and histogram which is not trivial.
    _accessor_rw.reset(new InternalZGY::ZgyInternalBulk
                       (_fd, _meta_rw, _meta_rw, compress,
                        iargs.historange_, _logger));

    enableIncrementalLOD(iargs.decimation_, iargs.internal_lod_mode_);
  }

  /**
   * \brief Open an existing file for update.
   * \details
   * See the ZgyReader constructor and the "truncate" ZgyWriter constructor.
   * Understandably this nethod contains elements from both.
   */
  ZgyWriter(const ZgyWriterArgsV3& args, OpenForUpdate)
    : _fd()
    , _accessor_rw()
    , _compressor(args.impl().compressor_)
    , _lodcompressor(args.impl().lodcompressor_ ?
                     args.impl().lodcompressor_ :
                     args.impl().compressor_)
  {
    const InternalZGY::ZgyWriterArgsV3Impl& iargs(args.impl());
    _myname = makemyname(iargs.filename_);
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(iargs.iocontext_.get());
    _logger(1, "Open U " + myname());
    _locked_in_process.reset(new LockedInProcess(iargs.filename_, true, _logger));

    // TODO-WIP-BrickedAPI: BUG: compress should also have been true
    // if the file was written compressed at some time but isn't now.
    // "compress" controls whether overwriting real data is allowed.
    // update mode still needs to be UpdateMode::Constant, not
    // UpdateMode::Always, because overwriting compressed data with
    // uncompressed will also leak the old brick. See Unit test
    // api.lodmode9.
    const bool compress = _compressor || _lodcompressor;

    _fd = InternalZGY::FileFactory::instance().create(
        iargs.filename_, InternalZGY::OpenMode::ReadWrite, ctxt.get());

    // Set the protected metadata information in the ZgyMeta base class.
    // The metadata ctor is very different on update. bulkdata is not.
    // ZgyInternalMeta does not retain the file descriptor. The file is
    // only used inside the constructor to populate the headers.

    _meta_rw.reset(new InternalZGY::ZgyInternalMeta(_fd, iargs, compress, _logger));
    _meta = _meta_rw;

    checkValidLodMode(iargs.internal_lod_mode_, _meta_rw.get());

    if (compress && _meta->ih().datatype() != InternalZGY::RawDataType::Float32)
      throw Errors::ZgyUserError("Compressed files need to be stored as float.");

    // At the implementation level the bulk- and meta access are separate,
    // and the bulk accessor needs some of the meta information to work.
    // The accessor will have its own metadata pointer which it holds on to.
    // It also holds on to the file descriptor. Unlike the metadata reader
    // it is obviously not possible to read everything up front.
    _accessor_rw.reset(new InternalZGY::ZgyInternalBulk
                       (_fd, _meta_rw, _meta_rw, compress,
                        iargs.historange_, _logger));

    // Consistency checks: Only files uploaded by OpenZGY can be updated.
    // See also the consistency checks in the SeismicStoreFileDelayedWrite
    // constructor regarding the segment size.
    if (_fd->xx_iscloud()) {
      const std::int64_t headersize = this->_meta_rw->flushMeta(nullptr);
      const std::shared_ptr<const FileStatistics> fs = filestats();
      const std::vector<std::int64_t> segsizes = fs->segmentSizes();
      if (segsizes.size() != 0) {
        if (fs->dataStart() >= 0 && fs->dataStart() < segsizes[0]) {
          // One or more bricks or tiles were found in the first segment.
          // Most likely the file was uploaded by sdutil in a single chunk,
          // or it may have been written by the old ZGY-Cloud.
          // Distinguishing those two is not always possible because
          // ZGY-Cloud can also put everything in the same segment.
          throw Errors::ZgyUpdateRules
            ("Only files uploaded by OpenZGY can be updated.");
        }
        if (headersize != segsizes[0]) {
          // Even when there is no data in the header area, the header
          // segment must be exactly the expected size. Most likely
          // this is a file containing no data and uploaded by sdutil
          // or the old ZGY-Cloud. If there is more than one segment
          // or eof is > headersize then there is something weird going on.
          // Probably not useful to report on that case though.
          throw Errors::ZgyUpdateRules
            ("Only files uploaded by OpenZGY can be updated. Bad Header size.");
        }
      }
    }

    enableIncrementalLOD(iargs.decimation_, iargs.internal_lod_mode_);
  }

  /**
   * \brief Automatically close the file when it goes out of scope.
   *
   * Application code is encouraged to close the file explicitly.
   * The destructor is just there as a fallback. Errors caught
   * in the fallback will be logged to stderr and otherwise ignored.
   */
  virtual ~ZgyWriter()
  {
    _logger(1, "Destroy W " + myname() + std::string((_accessor_rw || _fd) ? "" : " (no-op)"));
    try {
      close();
    }
    catch (const Errors::ZgyError& ex) {
      _logger(0, "ERROR: Uncaught ZGY exception closing file: " +
              std::string(ex.what() ? ex.what() : "(null)"));
    }
    catch (const std::exception& ex) {
      _logger(0, "ERROR: Uncaught general exception closing file: " +
              std::string(ex.what() ? ex.what() : "(null)"));
    }
  }

private:

  /**
   * Set up instances for writing low resolution data in parallel with lod0.
   * To be called from both constructors. Must NOT be virtual.
   */
  void enableIncrementalLOD(
       const std::vector<DecimationType>& decimation_in,
       const InternalZGY::InternalLodMode& internal_lod_mode)
  {
    // ZgyInternalBulk::writeAlignedBrickList() needs to create lowres
    // bricks as soon as possible. It is easier to hold on to a single
    // genlod instance to use as long as the file remains open. Create
    // it here, not in ZgyInternalBulk, because there is information
    // that genlod holds on to that ZgyInternalBulk should not know
    // about.
    // TODO-WIP-BrickedAPI: It is all to easy to reverse force and level,
    // consider wrapping them in a struct.
    // If LodMode::Never, all force & level are zero and tryToCall()
    // is a no-op. So, there is no need for the instance at all.
    if ((internal_lod_mode.incr.force != 0 && internal_lod_mode.incr.level != 0) ||
        (internal_lod_mode.last.force != 0 && internal_lod_mode.last.level != 0)) {
      // Share the stats/histo that ZgyInternalBulk created and owns.
      std::shared_ptr<InternalZGY::StatisticData> stats;
      std::shared_ptr<InternalZGY::HistogramData> histo;
      std::tie(stats, histo) = _accessor_rw->newTrackedChanges();
      if (stats && histo) {
        _logger(2, "Set new GenLod instance.");
        // Environment variables for LOD algorithms have precedence over
        // both the fallback and any explicitly specified overrides.
        const std::vector<DecimationType> decimation =
          _lodAlgoFromEnvironment(decimation_in);

        if (_logger(2, "")) {
          std::stringstream ss;
          ss << "enableIncrementalLOD decimation in [";
          for (const auto& it : decimation_in)
            ss << " " << (int)it;
          ss << "] actual [";
          for (const auto& it : decimation)
            ss << " " << (int)it;
          ss << "]\n";
          ss << "LOD incr force " << internal_lod_mode.incr.force
             << " level " << internal_lod_mode.incr.level
             << ", LOD last force " << internal_lod_mode.last.force
             << " level " << internal_lod_mode.last.level;
          _logger(2, ss.str());
        }

        // Save maxlevel & force for incremental & final in the bulk
        // instance, in the same way the the genlod instance is handled.
        // Pass the information to the tryToCall() method.
        _accessor_rw->newSetGenLodInstance
          (std::make_shared<InternalZGY::GenLodC>
           (_accessor_rw, _meta_rw, _lodcompressor,
            EnumMapper::mapDecimationTypeToLodAlgorithm(decimation),
            histo, _logger),
           internal_lod_mode);
      }
      else {
        _logger(1, "No new GenLod instance. New tracker is disabled.");
      }
    }
    else {
      _logger(1, "No new GenLod instance. Disabled by user.");
    }
  }

  static LodMode
  getUserLodMode(const InternalZGY::InternalLodMode& lodmode)
  {
    if ((lodmode.incr.level == 0 || lodmode.incr.force == 0) &&
        (lodmode.last.level == 0 || lodmode.last.force == 0))
      return LodMode::Never; // Guaranteed nothing can be produced.

    if (lodmode.last.level < 0 && lodmode.last.force >= 3)
      return LodMode::Rebuild; // In finalize. Clobbers data from incr level.

    if (lodmode.last.level < 0 && lodmode.last.force == 2)
      return LodMode::Early; // Or Early1, but we don't care.

    // The user must have set the 4 individual variables.
    // The settings might not make sense, but we know that *some*
    // lowres will be produced, and it isn't a rebuild.
    // So: Early, Early1, Build (only dirty) only in finalize,
    // or a mode that probably writes an inconsistent file.
    // Either way, treating ity as Early is reasonable.

    return LodMode::Early;
  }

  /**
   * \brief See if the requested LodMode makes sense for a file to be updated.
   *
   * \details
   * This is surprisingly difficult.
   *
   * See LookupTable::hasBrickLOD() that does something similar.
   *
   * Using the following rules:
   *   - Only allow the predefined LodMode tags, not the lower level settings.
   *   - Do not allow application to leave lowres bricks present but stale.
   *   - Do not allow application to delete existing lowres with LodMode::Never.
   *     This would hardly be useful, and tricky to implement. Potentially a
   *     serious performance hog and/or leak when clearing all lowres bricks.
   *   - On open for update, trust that none of the existing lowres bricks
   *     are stale. If the application suspects stale bricks then
   *     LodMode::Rebuild must be requested. Beware that files written by
   *     an older OpenZGY version might have had its lowres logically
   *     deleted. That is an obscure case where LodMode::Rebuild is needed.
   *   - If an invalid LodMode is spotted, an exception will be thrown instead
   *     of silently switching to LodMode::rebuild. Caveat, the code won't
   *     always be able to catch invalid modes. E.g. in the above case.
   *
   * Implementation:
   *
   * <ol>
   * <li>
   * If the file is opened for create, not update, definitely all LodMode
   * are good. An explicit test might not be needed except for a tiny
   * performance improvement. This case is caught by the all-constant test.
   *
   * <li>
   * if the requested mode is LodMode::Rebuild, definitely ok.
   *
   * <li>
   * if nlods==1, no lowres is possible, definitely all LodMode is good.
   * The survey might or might not be empty.
   *
   * <li>
   * if lod[0..nlods-1] is surely all constant, surely ok for any LodMode::
   * Lowres may or may not have been generated;reading lowres will work
   * either way. Note that this doesn't check for inflated bricks. Might
   * end up stricter than needed due to the tests below.
   *
   * <li>
   * The single lod[nlods-1] non-constant: The file has lowres, and the
   * code will trust that it is not stale. All except LodMode::Never is
   * allowed. The latter would leave stale data behind.
   *
   *   - Obscure case: If lod[nlods-1] only appears to be non-constant due
   *     to the brick being inflated, the code will still trust that the
   *     rest isn't stale. This is hopefully a valid assumption.
   *
   *   - Obscure case: If lod[nlods-1] only appears to be non-constant and
   *     the entire survey is actually constant, LodMode::Never could have
   *     been accepred.
   *
   * <li>
   * if lod[1..nlods-1] is for sure all constant and level 0 is
   * apparently not, lowres is missing. This means that previous writes did
   * write actual data and elected to use LodMode::Never. Allow Never (keep
   * it missing) or Rebuild, but Early/Early1/Default won't work
   *
   *   - Obscure case: level0 might be all-const but have inflated bricks,
   *     The survey is logically empty. Refusing Early when we technically
   *     don't need to.
   *
   * <li>
   * lod[nlods-1] surely constant and lod[1..nlods-2] apparently
   * non-constant: This is ambiguous. It might be that:
   *
   *   - The decimation algorithm might end up with a constant value top
   *     brick, in spite of the lowres being generated normally. This means
   *     the file has lowres, and LodMode::Never is not allowed.
   *     Obscure issue: The old code might incorrectly believe the file
   *     had no lowres at all.
   *
   *   - Some constant but inflated lowres bricks exist. It is unclear what
   *     to do.
   *
   *   - Lowres was logically deleted by an older OpenZGY. Must use
   *     LodMode::Rebuild because the new code sees this as stale data.
   *
   * <li>
   * To avoid the ambiguity, put code in _close_internal() that switches
   * the top level brick to Missing when needed.
   *
   * <li>
   * If levels 1..nlods-1 is surely all constant and level 0 is surely all
   * constant but a different constant, lowres is stale or corrupt. Tell
   * the application to use LodMode::Rebuild. This is an obscure case,
   * probably not worth testing for. Just assume garbage in, garbage out is
   * acceptable.
   *
   * </ol>
   */
  static void
  checkValidLodMode(const InternalZGY::InternalLodMode& lodmode_in, const InternalZGY::ZgyInternalMeta* meta)
  {
    using InternalZGY::LookupTable;
    const LodMode mode = getUserLodMode(lodmode_in);

    if (mode == LodMode::Rebuild)
      return; // Always legal. Unless the sdms protests.

    if (meta == nullptr)
      return; // Not open for write.

    const LookupTable::FinalizedStatus status =
      LookupTable::getFinalizedStatus
      (meta->ih().lodsizes(),
       meta->ih().brickoffsets(),
       meta->blup().lup(),
       meta->blup().lupend());

    switch (status) {
    case LookupTable::FinalizedStatus::NotFinalized:
      if (mode != LodMode::Never && mode != LodMode::Rebuild)
        throw OpenZGY::Errors::ZgyUserError
          ("LodMode::Rebuild or Never needed when the file has no lowres.");
      break;

    case LookupTable::FinalizedStatus::IsFinalized:
      if (mode == LodMode::Never)
        throw OpenZGY::Errors::ZgyUserError
          ("LodMode::Never not allowed on a file that already has lowres.");
      break;

    case LookupTable::FinalizedStatus::SingleBrick:
    case LookupTable::FinalizedStatus::ConstantFile:
    default:
      break; // No restrictions.
    }
  }

  /**
   * \callgraph \callergraph
   * \copydoc OpenZGY::ZgyReader::read(const size3i_t&,const size3i_t&,float*,int)const
   *
   * This method allows reading lod 0 from a file opened for update.
   * Note that IZgyWriter currently does not inherit IZgyReader so the
   * read methods in those two are not the same. The signature also
   * differs because it is only possible to read full resolution data
   * while the file is written. So the lod parameter has been removed.
   *
   * \internal The reason for not inheriting IZgyReader and ZgyReader
   * is to allow the lod parameter to be removed. Also avoid multiple
   * inheritance. The downside is a small amount of code duplication.
   */
  void read(const size3i_t& start, const size3i_t& size, float* data) const override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_readable();
    std::shared_ptr<float> fakeshared = fake_shared(data);
    auto databuffer = std::make_shared<InternalZGY::DataBufferNd<float,3>>(fakeshared, size);
    accessor_rw->readToExistingBuffer(databuffer, start, /*lod*/0, true);
    databuffer.reset();
    if (fakeshared.use_count() != 1)
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "read from write accessor");
  }

  /**
   * \copydoc read(const size3i_t&,const size3i_t&,float*)const
   */
  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_readable();
    std::shared_ptr<std::int16_t> fakeshared = fake_shared(data);
    auto databuffer = std::make_shared<InternalZGY::DataBufferNd<std::int16_t,3>>(fakeshared, size);
    accessor_rw->readToExistingBuffer(databuffer, start, /*lod*/0, false);
    databuffer.reset();
    if (fakeshared.use_count() != 1)
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "read from write accessor");
  }

  /**
   * \copydoc read(const size3i_t&,const size3i_t&,float*)const
   */
  void read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_readable();
    std::shared_ptr<std::int8_t> fakeshared = fake_shared(data);
    auto databuffer = std::make_shared<InternalZGY::DataBufferNd<std::int8_t,3>>(fakeshared, size);
    accessor_rw->readToExistingBuffer(databuffer, start, /*lod*/0, false);
    databuffer.reset();
    if (fakeshared.use_count() != 1)
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "read from write accessor");
  }

  /**
   * \callgraph \callergraph
   * \copydoc ZgyReader::readconst
   *
   * This method allows reading lod 0 from a file opened for update.
   * Note that IZgyWriter currently does not inherit IZgyReader so the
   * read methods in those two are not the same. The signature also
   * differs because it is only possible to read full resolution data
   * while the file is written. So the lod parameter has been removed.
   *
   * \internal The reason for not inheriting IZgyReader and ZgyReader
   * is to allow the lod parameter to be removed. Also avoid multiple
   * inheritance. The downside is a small amount of code duplication.
   */
  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float) const override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_readable();
    auto result = accessor_rw->readConstantValue(start, size, /*lod*/0, as_float);
    ZombieCheck::checkNonUniquePtr(accessor_rw, "read from write accessor");
    return result;
  }

  /**
   * \brief Write an arbitrary region.
   *
   * \details \callgraph \callergraph
   * This will apply conversion float -> storage if needed.
   *
   * Data is ordered inline(slowest), crossline, vertical(fastest).
   *
   * A read/modify/write will be done if the region's start and size
   * doesn't align with bricksize. When writing to the cloud this
   * read/modify/write may incur performance and size penalties. So
   * do write brick aligned data if possible. The same applies to
   * writing compressed data where r/m/w can cause a severe
   * loss of quality.
   *
   * The start position refers to the specified lod level.
   * At lod 0 start + data.size can be up to the survey size.
   * At lod 1 the maximum is just half that, rounded up.
   */
  void write(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();

    // TODO-Worry: The buffer is supposed to be copied at least once
    // before being sent down to the lowest levels. The longer I wait
    // before I do that, the higher the risk is that it isn't done.
    // Currently there is code in _writeAlignedRegion() that ensures
    // the buffer was copied. But if I do the copy too early then it
    // might be wasted effort since r/m/w or value type conversion
    // would copy it anyway.
    //
    // Why I need a copy:
    //
    //   - I need an ugly const_cast because DataBuffer isn't fully const
    //     aware. The cast is particularly dangerous because the lower
    //     levels might be tempted to do padding and byteswapping in place.
    //     If the compiler doesn't stop this due to a "const" declaration
    //     we get subtle bugs.
    //
    //   - The user's buffer is not reference counted. Which means I need
    //     to use a fake shared_ptr. If the lower levels implement delayed
    //     write and are "smart" about not always copying buffers then
    //     the buffer might be accessed after the function returns.
    std::shared_ptr<float> fakeshared = fake_shared(const_cast<float*>(data));
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<float,3>(fakeshared, size));
    accessor_rw->writeRegion(buffer, start, 0, false, _compressor);
    buffer.reset();
    if (fakeshared.use_count() != 1) // Actually a fatal error.
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write an arbitrary region with no conversion.
   *
   * As the write overload with a float buffer but only works for files with
   * SampleDataType::int16 and does not scale the samples.
   */
  void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();
    std::shared_ptr<std::int16_t> fakeshared = fake_shared(const_cast<std::int16_t*>(data));
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<int16_t,3>(fakeshared, size));
    accessor_rw->writeRegion(buffer, start, 0, true, _compressor);
    buffer.reset();
    if (fakeshared.use_count() != 1) // Actually a fatal error.
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write an arbitrary region with no conversion.
   *
   * As the write overload with a float buffer but only works for files with
   * SampleDataType::int8 and does not scale the samples.
   */
  void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();
    std::shared_ptr<std::int8_t> fakeshared = fake_shared(const_cast<std::int8_t*>(data));
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<std::int8_t,3>(fakeshared, size));
    accessor_rw->writeRegion(buffer, start, 0, true, _compressor);
    buffer.reset();
    if (fakeshared.use_count() != 1) // Actually a fatal error.
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write all-constant data.
   *
   * \details \callgraph \callergraph
   * Works as the corresponding write but the entire region is set
   * to the same value. So the provided data buffer needs just one
   * value, or alternatively can be passed as &scalar_value.
   *
   * Calling this method is faster than filling a buffer with constant
   * values and calling write. But it produces the exact same
   * result. This is because write will automatically detect whether
   * the input buffer is all constant.
   */
  void writeconst(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<float,3>(*data, size));
    accessor_rw->writeRegion(buffer, start, 0, false, _compressor);
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write an arbitrary region with no conversion.
   *
   * As the writeconst overload with a float buffer but only works for files with
   * SampleDataType::int16 and does not scale the samples.
   */
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t* data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<std::int16_t,3>(*data, size));
    accessor_rw->writeRegion(buffer, start, 0, true, _compressor);
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write an arbitrary region with no conversion.
   *
   * As the write overload with a float buffer but only works for files with
   * SampleDataType::int8 and does not scale the samples.
   */
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    auto accessor_rw = _accessor_rw;
    throw_if_not_writable();
    std::shared_ptr<InternalZGY::DataBuffer> buffer(new InternalZGY::DataBufferNd<std::int8_t,3>(*data, size));
    accessor_rw->writeRegion(buffer, start, 0, true, _compressor);
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * Copy the on the fly statistics and histogram to the file header.
   * This is part of the new code.
   */
  void _store_stats_histo_new_version()
  {
    std::shared_ptr<InternalZGY::StatisticData> stats;
    std::shared_ptr<InternalZGY::HistogramData> histo;
    std::tie(stats, histo) = _accessor_rw->newTrackedChanges();
    if (!this->_accessor_rw->isStatHistGood()) {
      // Updating a file that had bad or missing stat/histo makes a
      // complete mess if bricks were updated, and is not really
      // useful even it bricks were only appended. It is better to
      // leave them unset on the bulk file. This should only happen
      // for files written by older versions of OpenZGY that were
      // not finalized. So, another approach would be to not allow
      // opening such files for update in the first place.
      _logger(0, "WARNING: Pre-existing statistics and histogram not found.");
      _logger(0, "WARNING: Statistics and histogram were not written.");
      return;
    }
    if (!stats || !histo) {
      // Currently this shouldn't happen.
      _logger(0, "WARNING: Statistics and histogram were not collected.");
      return;
    }

    InternalZGY::ExpandableBuilder hb(*histo, *stats, /*copy=*/true);
    const std::array<double,2> factors = this->_meta->ih().storagetofloat();
    hb.scale(0, 1, factors[1], factors[0] + factors[1]);
    InternalZGY::StatisticData scaled_stats = hb.getstats();
    InternalZGY::HistogramData scaled_histo = hb.finalize(256);
    this->_meta_rw->ih().setstats
      (scaled_stats.getcnt(), scaled_stats.getsum(), scaled_stats.getssq(),
       scaled_stats.getmin(), scaled_stats.getmax());
    this->_meta_rw->hh().sethisto(scaled_histo.getmin(), scaled_histo.getmax(),
                                  scaled_histo.getbins(), scaled_histo.getsize());
    if (_logger(2, "")) {
      std::stringstream ss;
      ss << "[New code] storing histogram Limits "
         << histo->getmin() << " " << histo->getmax()
         << " Converted "
         << scaled_histo.getmin() << " " << scaled_histo.getmax()
         << " counts s=" << scaled_stats.getcnt()
         << " h=" << scaled_histo.getcount();
      _logger(2, ss.str());
    }
  }

  std::string
  _formatListOfAlgo(const std::vector<DecimationType>& list)
  {
    if (list.empty())
      return "<empty>";
    std::stringstream ss;
    for (DecimationType algo : list)
      ss << "," << Formatters::enumToString(algo);
    return ss.str().substr(1);
  }

  /**
   * Override LOD algorithms from environment variables.
   * The envirionment also has precedence over algorithms
   * explicitly specified in the call to finalize().
   * Uses numbers not strings. Because this is just
   * a temporary hack anyway.
   */
  std::vector<DecimationType>
  _lodAlgoFromEnvironment(const std::vector<DecimationType>& decimation_in)
  {
    std::vector<DecimationType> decimation(decimation_in);

    if (decimation.size() < 1)
      decimation.push_back(DecimationType::LowPass);
    if (decimation.size() < 2)
      decimation.push_back(DecimationType::WeightedAverage);
    if (decimation.size() < 3)
      decimation.push_back(decimation.back());

    std::vector<DecimationType> result;
    result.push_back(static_cast<DecimationType>
                     (InternalZGY::Environment::getNumericEnv
                      ("OPENZGY_LODALGO_1", (int)decimation[0])));
    result.push_back(static_cast<DecimationType>
                     (InternalZGY::Environment::getNumericEnv
                      ("OPENZGY_LODALGO_2", (int)decimation[1])));
    result.push_back(static_cast<DecimationType>
                     (InternalZGY::Environment::getNumericEnv
                      ("OPENZGY_LODALGO_N", (int)decimation[2])));

    if (decimation[0] == result[0] &&
        decimation[1] == result[1] &&
        decimation[2] == result[2])
    {
      return decimation_in;
    }
    else {
      _logger(0, "LOD algorithms changed from " +
              _formatListOfAlgo(decimation_in) +
              " to " + _formatListOfAlgo(result));
      return result;
    }
  }

  /**
   * \brief Maybe generate low resolution data, statistics, and histogram.
   *
   * \details \callgraph \callergraph
   * This method will be called automatically from close(), but in
   * that case it is not possible to request a progress callback.
   *
   * \param progress:   Function(done, total) called to report progress.
   *                    If it returns False the computation is aborted.
   *                    The callback implementation should not invoke any
   *                    OpenZGY methods except for any calls expressly
   *                    documented as safe for this purpose and deadlock-free.
   *                    Will be called at least one, even if there is no work.
   *
   * \param decimation: Ignored.
   * \param action:     Ignored.
   * \param force:      Ignored.
   */
  void finalize(const std::vector<DecimationType>& /*decimation_in*/,
                const std::function<bool(std::int64_t,std::int64_t)>& progress,
                FinalizeAction /*action*/, bool /*force*/) override
  {
    // Called with error flag set is actually valid here.
    //throw_if_not_writable();
    if (!_fd || !_accessor_rw || !_meta_rw)
      throw Errors::ZgyUserError("ZGY file is not open for write.");
    _accessor_rw->newFinalize(progress);
    _store_stats_histo_new_version();
  }

  /**
   * \brief Flush the file to disk and close it.
   *
   * The application can request that a file should not contain
   * low resolution data. This setting has been moved to the
   * constructor.
   *
   * An incomplete file is still usable for reading full resolution
   * data as long as it is accessed using OpenZGY and not the old
   * ZGY-Public. OpenZGY reports nlods==1 for incomplete files
   * so applications can deal with it. ZGY-Public will not be allowed
   * to open incomplete files. Version is set >3 to prevent this.
   */
  void _close_internal()
  {
    _logger(1, "Close W " + myname() + std::string((_accessor_rw || _fd) ? "" : " (no-op)"));
    if (!this->_fd || !this->_accessor_rw || !this->_meta_rw) {
      _locked_in_process.reset(); // Just in case.
      return; // looks like _close_internal alrady called.
    }

    finalize(std::vector<DecimationType>{}, nullptr,
             FinalizeAction::BuildDefault, false);

    // Prevent ZGY-Public from opening the file if appropriate.
    this->_meta_rw->fh().set_version(_meta_rw->can_old_library_read() ? 3 : 4);

    {
      // If non-constant bulk exists, and if app didn't want to make
      // lowres (LodMode::Never), and the most decimated lodN brick
      // has type Constant, the file doesn't have valid lowres but
      // readers might think it has. See LookupTable::hasBrickLOD()
      // for an explanation. Switch that brick to Missing to avoid the
      // problem. Make the test as specific as possible, and do the
      // cheaper tests first. Note that the trick here means we need
      // an extra test in checkValidLodMode(). Sigh.
      if (this->_meta_rw->ih().lodsizes().size() > 1) {
        if (!this->_accessor_rw->wasLowResRequested()) {
          using namespace InternalZGY;
          const BrickStatus top_status =
            LookupTable::getBrickFilePosition
            (0, 0, 0, // i,j,k
             this->_meta_rw->ih().lodsizes().size()-1, // last lod
             this->_meta_rw->ih().lodsizes(),
             this->_meta_rw->ih().brickoffsets(),
             this->_meta_rw->blup().lup(),
             this->_meta_rw->blup().lupend(),
             /*bytesperbrick=*/0).status;
          if (top_status == BrickStatus::Constant) {
            const bool lod0_allsame = LookupTable::hasAllBricksSameValue
              (0, 0, this->_meta_rw->ih().brickoffsets(),
               this->_meta_rw->blup().lup(),
               this->_meta_rw->blup().lupend());
            if (!lod0_allsame) {
              LookupTable::setNoBrickLOD
                (this->_meta_rw->ih().lodsizes(),
                 this->_meta_rw->ih().brickoffsets(),
                 &this->_meta_rw->blup().lup(),
                 &this->_meta_rw->blup().lupend());
              _logger(1, "lodN switched to Missing to avoid confusion");
            }
          }
        }
      }
    }

    if (!this->errorflag()) {
      this->_meta_rw->flushMeta(this->_fd);
    }
    else {
      // Write the headers anyway, but don't complain if it didn't work.
      try {
        this->set_errorflag(false);
        this->_meta_rw->flushMeta(this->_fd);
        this->set_errorflag(true);
      }
      catch (std::exception&) {
        this->set_errorflag(true);
      }
    }

    // TODO-Low it would be more consistent if the final write was of the
    // updated headers, in case the app crashes before xx_close. This
    // is true for local file access. But for cloud access the last
    // bulk data segment might be written on xx_close() because it is
    // buffered. While the headers can be written at once. Difficult
    // to change without complicating the internal IFileADT api.

    // Closing the local or cloud handle is always needed even if
    // there was an unrecoverable error. There might be resources that
    // need to be cleaned up. TODO-Low if any of the above writes raised
    // am exception I should probably still try to close.

    // Clearing _fd is also needed to ensure that a subsequent close()
    // that might be triggered from a destructor becomes a no-op.
#if 0
    // TODO, also drop the reference to the accessor?
    // This seems to be an oversight, but risky to change.
    auto accessor_rw = _accessor_rw;
    if (accessor_rw) {
      _accessor_rw.reset();
      ZombieCheck::checkUniquePtr(accessor_rw, "ZgyWriter close");
      accessor_rw.reset();
    }
#endif

    this->_fd->xx_close();
    this->_fd.reset();
    this->_locked_in_process.reset();
    // ZombieCheck::checkUniquePtr(victim, "ZgyWriter close file");

    // Kludge for performance measurements,
    // The timers in DataBuffer are global. For some experiments
    // it makes sense to print and reset then just before finalize
    // and also after closing each file.
    // The same applies to the histogram builder, and others.
    // If multiple unrelated files are open at the same time,
    // these measurements might not be useful.
    InternalZGY::DataBuffer::cleartimers(true);
    InternalZGY::HistogramBuilder::cleartimers();

    // If errorflag() is set and the file is new or has been
    // successfully written to at least once then the client code is
    // strongly advised to delete the file.
    // TODO-Low: Later: this is in the bells & whistles category.
    //if not self._precious_set_on_open and was_written_to
    //    self._fd.xx_delete_on_close(); self._fd.xx_close()
    //    self._fd.xx_close_if_needed_and_delete()
    //    ZgyUtils(saved_iocontext).delete(self._filename)
  }

  /**
   * \brief Flush the file to disk and close it.
   *
   * There is no longer any difference between close_incomplete()
   * and close(). Statistics, histogram, and lowres are all computed
   * while the fullres is written. So there is not much that can be
   * saved by not writing them out.
   *
   * TODO-WIP-BrickedAPI: A possible new feature could be an abandon()
   * function that corrupts the file by closing it immediately. Useful
   * if the application has decided that the file is useless and
   * doesn't want to wait. Do NOT re-purpose close_incomplete() to do
   * this. Existing applications might be surprised.
   */
  void close_incomplete() override
  {
    _close_internal();
  }

  /**
   * \brief Flush the file to disk and close it.
   *
   * If the file has been written to, the application is encouraged to
   * call finalize() before close(). This gives more control over the
   * process and allows using a progress callback to track generation
   * of any remaining low resolution data.
   *
   * If the application fails to call finalize() it will be done
   * automatically here.
   *
   * The function won't bother with statistics, histogram, lowres if
   * there has been an unrecoverable error. The headers might still be
   * written out in case somebody wants to try some forensics.
   *
   * The ZgyWriter destructor will call close() if not done already,
   * but that will catch and swallow any exception. Relying on the
   * destructor to close the file is strongly discouraged.
   */
  void close() override
  {
    // TODO-@@@: If the file has never been written to and the error
    // flag is set then discard everyhing and do NOT write any data.
    // This can in some cases avoid corrupting a file that was opened
    // for write and then has an error thrown.
    // The same logic may be needed in _close_internal.
    _close_internal();
  }

  /**
   * Return true if this open file has encountered an unrecoverable error.
   * The error should previously have caused an exception to be thrown.
   * If this flag is set, no further writes will be allowed.
   *
   * Application code might check this flag if they are considering
   * trying to recover from an error. Internally the flag is also checked
   * and if set it will (mostly) prevent other writes from being done.
   *
   * Implementation note: Currently the ZgyInternalMeta and ZgyInternalBulk
   * instances each contains an _is_bad member. The reader or writer is
   * considered bad if either of those are set. This scheme improves
   * isolation somewhat, but TODO-Low it might backfire. If writing metadata
   * to file failed, the bulk accessor should probably behave as if it
   * also has seen an error.
   *
   * Currently only the ZgyWriter uses this mechanism. It might not make
   * that much sense in ZgyReader, because as long as opening the file
   * succeded no operation should manage to corrupt it.
   *
   * Implementation note: Unlike most other ZgyWriter members, this one is
   * declared const. But it can still access the mutable _accessor_rw and
   * _meta_rw data members. You'll just have to trust me when I tell you
   * they aren't being modified here. Or refactor somehow to allow the
   * compiler (or at least a runtime check) catch it.
   */
  bool errorflag() const override
  {
    return _accessor_rw->errorflag() || _meta_rw->errorflag();
  }

  /**
   * Force the error flag for this open file to true or false.
   * This should normally be done only for testing.
   */
  void set_errorflag(bool value) override
  {
    _accessor_rw->set_errorflag(value);
    _meta_rw->set_errorflag(value);
  }

  /**
   * Get the file statistics of a file currently opened for write.
   * There is no caching because the result will change whenever data
   * is written. If this is a problem then implement a cache as done
   * in ZgyReader and arrange for it to be cleared every time the
   * metadata is touched.
   */
  std::shared_ptr<const FileStatistics> filestats() const override
  {
    std::shared_ptr<FileStatistics> result
      (new FileStatistics(*filestats_nocache()));
    // The base class has no _fd member so I need to set the size here.
    result->_file_size = _fd->xx_eof();
    result->_segment_sizes = _fd->xx_segments(false);
    result->_data_start = std::min(result->_data_start, result->_file_size);
    return result;
  }

private:
  /**
   * Test for common user errors as early as possible.
   *
   * Checking that
   * the file is still open for write is particularly important. Some
   * write operations are deferred, possibly because they only affect
   * the lookup table. An error in the application would be caught
   * eventually by IFileADT::xx_write() but reported way too late.
   *
   * Also, a read/modify/write might cause the error to be reported in
   * the read step. Which is confusing.
   */
  void throw_if_not_writable()
  {
    if (errorflag())
      throw Errors::ZgyCorruptedFile("Cannot continue due to previous errors.");
    if (!_fd || !_accessor_rw || !_meta_rw)
      throw Errors::ZgyUserError("ZGY file is not open for write.");
  }

  void throw_if_not_readable() const
  {
    if (!_fd || !_accessor_rw || !_meta_rw)
      throw Errors::ZgyUserError("ZGY file is not open for read.");
  }
};

/**
 * \brief Concrete implementation of IZgyUtils.
 *
 * \details Thread safety: Depends on the function being executed.
 * Currently implemented methods just forward the requests to SDAPI
 * or some other cloud back end. So thread safety depends on the
 * cloud provider.
 *
 * \internal TODO-Low: One utility class per backend plug-in.
 * This will be a pain to maintain though; as the
 * Python wrapper will need to do the same.
 */
class ZgyUtils : public IZgyUtils
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

private:
  std::shared_ptr<InternalZGY::IFileBase> _fd;
  LoggerFn _logger;

public:
  /**
   * \copydoc IZgyUtils::utils()
   */
  ZgyUtils(const std::string& prefix, const IOContext* iocontext)
  {
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(iocontext);
    _logger(1, "open utility class");
    _fd = InternalZGY::FileFactory::instance().create
      (prefix, InternalZGY::OpenMode::Closed, ctxt.get());
  }

  void deletefile(const std::string& filename, bool missing_ok) override
  {
    _fd->deleteFile(filename, missing_ok);
  }

  std::string alturl(const std::string& filename) override
  {
    return _fd->altUrl(filename);
  }

  std::string idtoken() override
  {
    return _fd->idToken();
  }

  std::map<std::string,std::string> basicinfo(const std::string& filename) override
  {
    return _fd->getBasicInfo(filename);
  }
};

#if 0

// TODO-Low: Refactor to a cleaner way of choosing a compressor.
// Since the factories have variable argument lists I might
// not be able to encapsulate as much as I do in the Python version.
// Ummm... didn't I fix that by using a string list?

//def ZgyCompressFactory(name, *args, **kwargs):
//    return _internal.CompressFactoryImpl.factory(name, *args, **kwargs)
//

std::vector<std::string>
ZgyKnownCompressors()
{
  throw std::runtime_error("Not implemented: ZgyKnownCompressors()");
  //return _internal.CompressFactoryImpl.knownCompressors()
}

std::vector<std::string>
ZgyKnownDecompressors()
{
  throw std::runtime_error("Not implemented: ZgyKnownDecompressors()");
  //return _internal.CompressFactoryImpl.knownDecompressors()
}

#endif

} // namespace OpenZGY::Impl
// Remain in namespace OpenZGY

/**
 * \brief Create a fixed range empty histogram.
 *
 * \internal
 *
 * - Could, but doesn't, simplify by using the internal ExtendableBuilder
 *   also when fixed-width is requested. That way of using ExtendableBuilder
 *   should work. But it has not been tested as much as the rest.
 *
 * - Could, but doesn't, use the fact that the internal ExtendableBuilder
 *   inherits HistogramBuilder. Using that would be fragile. Because
 *   add(), operator+=, and operator-= are overridden but not virtual.
 */
SampleHistogramBuilder::SampleHistogramBuilder(size_type nbins, double min, double max)
  : fixed_(std::make_shared<InternalZGY::HistogramBuilder>(nbins, min, max))
  , dynamic_(nullptr)
{
}

SampleHistogramBuilder::SampleHistogramBuilder(size_type nbins)
  : fixed_(nullptr)
  , dynamic_(std::make_shared<InternalZGY::ExpandableBuilder>(nbins))
{
}

SampleHistogramBuilder&
SampleHistogramBuilder::operator+=(const SampleHistogramBuilder& other)
{
  if (this->dynamic_) {
    if (other.dynamic_)
      *this->dynamic_ += *other.dynamic_;
    else
      *this->dynamic_ += *other.fixed_;
  }
  else {
    if (other.dynamic_)
      *this->fixed_ += *other.dynamic_;
    else
      *this->fixed_ += *other.fixed_;
  }
  return *this;
}

SampleHistogramBuilder&
SampleHistogramBuilder::operator-=(const SampleHistogramBuilder& other)
{
  if (this->dynamic_) {
    if (other.dynamic_)
      *this->dynamic_ -= *other.dynamic_;
    else
      *this->dynamic_ -= *other.fixed_;
  }
  else {
    if (other.dynamic_)
      *this->fixed_ -= *other.dynamic_;
    else
      *this->fixed_ -= *other.fixed_;
  }
  return *this;
}

SampleHistogramBuilder&
SampleHistogramBuilder::operator*=(count_type factor)
{
  if (dynamic_)
    *dynamic_ *= factor;
  else
    *fixed_ *= factor;
  return *this;
}

SampleStatistics
SampleHistogramBuilder::getstats() const
{
  const InternalZGY::StatisticData& s =
    dynamic_ ? dynamic_->getstats() : fixed_->getstats();
  return SampleStatistics
    (s.getcnt(), s.getsum(), s.getssq(), s.getmin(), s.getmax());
}

SampleHistogram
SampleHistogramBuilder::gethisto() const
{
  const InternalZGY::HistogramData& hh =
    dynamic_ ? dynamic_->gethisto() : fixed_->gethisto();
  if (hh.getcount() > 0 && hh.getbins() != nullptr)
    return SampleHistogram
      (hh.getcount(), hh.getmin(), hh.getmax(),
       std::vector<std::int64_t>(hh.getbins(),
                                 hh.getbins() + hh.getsize()));
  else
    return SampleHistogram();
}

SampleHistogram
SampleHistogramBuilder::finalize(size_type want_nbins) const
{
  if (dynamic_) {
    InternalZGY::HistogramData hh = dynamic_->finalize(want_nbins);
    if (hh.getcount() > 0 && hh.getbins() != nullptr) {
      return SampleHistogram
        (hh.getcount(), hh.getmin(), hh.getmax(),
         std::vector<std::int64_t>(hh.getbins(),
                                   hh.getbins() + hh.getsize()));
    }
  }
  return gethisto();
}

void
SampleHistogramBuilder::add(const std::int8_t  *begin, const std::int8_t  *end)
{
  if (dynamic_)
    dynamic_->add(begin, end);
  else
    fixed_->add(begin, end);
}

void
SampleHistogramBuilder::add(const std::int16_t *begin, const std::int16_t *end)
{
  if (dynamic_)
    dynamic_->add(begin, end);
  else
    fixed_->add(begin, end);
}

void
SampleHistogramBuilder::add(const float        *begin, const float        *end)
{
  if (dynamic_)
    dynamic_->add(begin, end);
  else
    fixed_->add(begin, end);
}

void
SampleHistogramBuilder::scale(double oldmin, double oldmax, double newmin, double newmax)
{
  if (dynamic_)
    dynamic_->scale(oldmin, oldmax, newmin, newmax);
  else
    fixed_->scale(oldmin, oldmax, newmin, newmax);
}

class ProgressWithDots::Impl
{
public:
  int dots_printed_;
  int length_;
  std::ostream& outfile_;
  std::mutex mutex_;
  Impl(int length, std::ostream& outfile)
  : dots_printed_(0)
  , length_(length)
  , outfile_(outfile)
  , mutex_()
  {}
};

ProgressWithDots::ProgressWithDots(int length, std::ostream& outfile)
  : pimpl_(std::make_shared<Impl>(length, outfile))
{
}

bool
ProgressWithDots::operator()(std::int64_t done, std::int64_t total)
{
  if (pimpl_->length_ < 1)
    return true;
  std::lock_guard<std::mutex> lk(pimpl_->mutex_);
  if (pimpl_->dots_printed_ == 0) {
    pimpl_->outfile_ << "[" + std::string(pimpl_->length_, ' ') + "]\r[" << std::flush;
  }
  std::int64_t needed = (total <= 0) ? 1 : 1 + ((done * (pimpl_->length_-1)) / total);
  if (needed > pimpl_->dots_printed_) {
    while (needed > pimpl_->dots_printed_) {
      pimpl_->outfile_ << '.';
      pimpl_->dots_printed_ += 1;
    }
    pimpl_->outfile_ << std::flush;
  }
  if (done == total)
    pimpl_->outfile_ << "\n" << std::flush;
  return true;
}

class FancyProgressWithDots::Impl
{
public:
  const int     length_;
  std::ostream& outfile_;
  std::mutex    mutex_;
  double        starttime_;
  double        prev_starttime_;
  int           prev_printed_;
  std::int64_t  prev_done_;
  int           output_count_;

  Impl(int length, std::ostream& outfile)
    : length_(length)
    , outfile_(outfile)
    , mutex_()
    , starttime_(timestamp())
    , prev_starttime_(starttime_)
    , prev_printed_(0)
    , prev_done_(0)
    , output_count_(0)
  {}
};

FancyProgressWithDots::FancyProgressWithDots(int length, std::ostream& outfile)
  : pimpl_(std::make_shared<Impl>(length, outfile))
{
}

static void
to_hm(double seconds, std::stringstream& ss, bool do_seconds = false)
{
  int sec = do_seconds ? (int)ceil(seconds) : (int)ceil(seconds/60.0)*60;
  int min = sec / 60;
  sec -= min * 60;
  int hour = min / 60;
  min -= hour * 60;
  if (seconds < 0)
    ss << " --:--";
  else if (!do_seconds)
    ss << " " << std::setw(2) << std::setfill(' ') << hour
       << ":" << std::setw(2) << std::setfill('0') << min
       << std::setfill(' ');
  else if (hour == 0)
    ss << " " << std::setw(2) << std::setfill(' ') << min
       << ":" << std::setw(2) << std::setfill('0') << sec
       << std::setfill(' ');
  else
      ss << " " << std::setw(2) << std::setfill(' ') << hour
         << ":" << std::setw(2) << std::setfill('0') << min
         << ":" << std::setw(2) << std::setfill('0') << sec
         << std::setfill(' ');
}

bool
FancyProgressWithDots::operator()(std::int64_t done, std::int64_t total)
{
  Impl& p(*this->pimpl_);
  if (p.length_ < 1)
    return true;
  std::lock_guard<std::mutex> lk(p.mutex_);
  const int needed =
    (total <= 0) ? 1 : 1 + (int)((done * (p.length_-1)) / total);
  if (needed == p.prev_printed_ && timestamp() - p.prev_starttime_ < 20.0)
    return true;

  // Timing
  if (done == 0)
    p.starttime_ = timestamp();
  //const std::int64_t this_done  = done - p.prev_done_;
  const double this_endtime     = timestamp();
  //const double this_elapsed     = this_endtime - p.prev_starttime_;
  const double all_elapsed      = this_endtime - p.starttime_;
  //const double average_all_step = done<=0 ? 0 : all_elapsed / done;
  //const double average_per_step = this_done<=0 ? 0 : this_elapsed / this_done;
  const double est_totaltime    = done <= 0 ? 0 : all_elapsed * ((double)total / (double)done);
  const double est_remain       = est_totaltime - all_elapsed;
  p.prev_starttime_ = this_endtime;
  p.prev_printed_   = needed;
  p.prev_done_      = done;

  // Progress bar
  std::stringstream ss;
  ss << '\r' << '[';
  for (int ii = 0; ii < p.length_; ++ii)
    ss << (ii >= needed ? ' ' : ii % 5 == 0 ? '+' : '.');
  ss << ']';

  // Timing
  // What is most useful?
  //  - Total elapsed time.
  //  - Estimated time left (ETL) based on average over entire run.
  //  - Estimated total time based on average over entire run.
  //  - Time for this 2% step vs. the average time for the same period.
  //  - Fancier estimation based on moving average (not implemnented)
  //  - Wall clock ETA.
  //  - A single char showing something happening every 10 seconds or so,
  //ss << " Run";
  to_hm(all_elapsed, ss, true);
  //ss << " ETL ";
  to_hm(est_remain, ss, true);
  //ss << " Est. total";
  //to_hm (est_totaltime, ss);
  //ss << " average " << average_all_step << " " << average_per_step;
  ss << " ";
  //ss << ((p.output_count_ % 2) == 0 ? '*' : '+');
  //ss << " ";
  //ss << "done " << done << "/" << total << " ";
  if (done == total)
    ss << '\n';
  p.outfile_ << ss.str() << std::flush;
  ++p.output_count_;
  return true;
}

double
FancyProgressWithDots::timestamp()
{
  // Remember the good old days with ::time(nullptr) ?
  typedef std::chrono::high_resolution_clock clock;
  static constexpr double factor =
    static_cast<double>(clock::duration::period::num) /
    static_cast<double>(clock::duration::period::den);
  return clock::now().time_since_epoch().count() * factor;
}

// Dummy destructors for the interface types.
// Some compilers like to have at least one non-abstract member in
// each class, implemented outside the header file, as this gives
// the compiler an obvious place to put the vtbl.

IZgyMeta::~IZgyMeta() {}
IZgyTools::~IZgyTools() {}
IZgyReader::~IZgyReader() {}
IZgyWriter::~IZgyWriter() {}
IZgyUtils::~IZgyUtils() {}

// Not-quite factory methods for opening a file and getting a handle
// for it. Lexically scoped inside the pure interface classes but
// don't really belong there. Free functions would also work.

std::shared_ptr<IZgyReader>
IZgyReader::open(const std::string& filename, const IOContext* iocontext)
{
  std::shared_ptr<IZgyReader> reader =
    std::shared_ptr<IZgyReader>(new Impl::ZgyReader(filename, iocontext));
  reader = Telemetry::ApiReaderRecorder::inject(reader, filename);
  return reader;
}

std::shared_ptr<IZgyWriter>
IZgyWriter::open(const ZgyWriterArgsV3& args)
{
  auto writer = std::shared_ptr<IZgyWriter>(new Impl::ZgyWriter(args));
  writer = std::shared_ptr<IZgyWriter>(new Impl::ZgySafeWriter(writer));
  writer = Telemetry::ApiWriterRecorder::inject(writer, args.impl().filename_);
  return writer;
}

std::shared_ptr<IZgyWriter>
IZgyWriter::open(const ZgyWriterArgsV2& oldargs)
{
  ZgyWriterArgsV3 args(oldargs, ZgyWriterArgsV3::FromOldWriterArgsV2{});
  return open(args);
}

std::shared_ptr<IZgyWriter>
IZgyWriter::open(const ZgyWriterArgs& oldargs)
{
  ZgyWriterArgsV3 args(oldargs, ZgyWriterArgsV3::FromOldWriterArgs{});
  return open(args);
}

/**
 * \brief Open an existing ZGY file for writing.
 *
 * \details
 * There are several restrictions on using this feature.
 * Some of those might be relaxed in the future.
 *
 * \li The file must have the latest version. Currently not enforced.
 *
 * \li size, bricksize, datatype, and datarange are not allowed to change.
 *
 * \li filename, iocontext, compressor, and lodcompressor are not stored
 *     in the file itself. So the caller needs to re-specify them.
 *
 * \li ilstart, ilinc, xlstart, xlinc, zstart, zinc, corners
 *     can all be changed but note that this can only be done in the
 *     actual call to reopen().
 *
 * \li zunit and hunit cannot be changed because the string table that
 *     holds the unit names is still read only. I don't think this will
 *     be a problem.
 *
 * \li If the file to be updated is stored on a cloud back-end then it must
 *     have been created by OpenZGY and written directly to the cloud.
 *     The reason: Even though the ZGY-Public and OpenZGY file formats are
 *     identical, the data on the cloud can be split into multiple chunks
 *     which are logically treated as a single byte stream. The choice
 *     of where to split differs between OpenZGY, ZGY-Public, and uploads
 *     done by sdutil.
 *
 * \li If the file does contain bulk data
 *     then it is strongly advised to not change from compressed to
 *     uncompressed or vice versa. Even though this might not be checked.
 *     If you ignore this advice then he new compression mode will only
 *     apply to subsequent writes.
 *
 * \li If the file does not contain bulk data but is flagged as v4
 *     anyway, it is ok to reopen it as uncompressed. In that case
 *     it is unspecifed whether the version is changed back from 4 to 3
 *     (which is the current behavior) or not.
 *     Caveat: If the file was originally created as compressed then
 *     bricks are not guaranteed to be stored aligned to the brick size.
 *
 * \li If the file does contain bulk data
 *     then it is strongly advised to not update any already existing data.
 *     Including overwrite due to a read/modify/write cycle. With a cloud
 *     backend this could waste storage and might also be refused at runtime.
 *
 * \li The operation can leave the file unusable or even deleted if some
 *     error occurs.
 */
std::shared_ptr<IZgyWriter>
IZgyWriter::reopen(const ZgyWriterArgsV3& args)
{
  auto writer = std::shared_ptr<IZgyWriter>(new Impl::ZgyWriter(args, Impl::ZgyWriter::OpenForUpdate{}));
  writer = std::shared_ptr<IZgyWriter>(new Impl::ZgySafeWriter(writer));
  writer = Telemetry::ApiWriterRecorder::inject(writer, args.impl().filename_);
  return writer;
}

std::shared_ptr<IZgyWriter>
IZgyWriter::reopen(const ZgyWriterArgsV2& oldargs)
{
  ZgyWriterArgsV3 args(oldargs, ZgyWriterArgsV3::FromOldWriterArgsV2{});
  return reopen(args);
}

std::shared_ptr<IZgyWriter>
IZgyWriter::reopen(const ZgyWriterArgs& oldargs)
{
  ZgyWriterArgsV3 args(oldargs, ZgyWriterArgsV3::FromOldWriterArgs{});
  return reopen(args);
}

std::shared_ptr<IZgyUtils>
IZgyUtils::utils(const std::string& prefix, const IOContext* iocontext)
{
  return std::shared_ptr<IZgyUtils>(new Impl::ZgyUtils(prefix, iocontext));
}

ZgyWriterArgs&
ZgyWriterArgs::metafrom(const std::shared_ptr<OpenZGY::IZgyReader>& reader)
{
  OpenZGY::IZgyReader *r = reader.get();
  size(r->size()[0], r->size()[1], r->size()[2]);
  bricksize(r->bricksize()[0], r->bricksize()[1], r->bricksize()[2]);
  datatype(r->datatype());
  datarange(r->datarange()[0], r->datarange()[1]);
  zunit(r->zunitdim(), r->zunitname(), r->zunitfactor());
  hunit(r->hunitdim(), r->hunitname(), r->hunitfactor());
  ilstart(r->annotstart()[0]);
  ilinc(r->annotinc()[0]);
  xlstart(r->annotstart()[1]);
  xlinc(r->annotinc()[1]);
  zstart(r->zstart());
  zinc(r->zinc());
  corners(r->corners());
  return *this;
}

ZgyWriterArgs&
ZgyWriterArgs::merge(const ZgyWriterArgs& other)
{
  if (other._have_size)
    size(other._size[0], other._size[1], other._size[2]);
  if (other._have_bricksize)
    bricksize(other._bricksize[0], other._bricksize[1], other._bricksize[2]);
  if (other._have_datatype)
    datatype(other._datatype);
  if (other._have_datarange)
    datarange(other._datarange[0], other._datarange[1]);
  if (other._have_zunit)
    zunit(other._zunitdim, other._zunitname, other._zunitfactor);
  if (other._have_hunit)
    hunit(other._hunitdim, other._hunitname, other._hunitfactor);
  if (other._have_ilstart)
    ilstart(other._annotstart[0]);
  if (other._have_ilinc)
    ilinc(other._annotinc[0]);
  if (other._have_xlstart)
    xlstart(other._annotstart[1]);
  if (other._have_xlinc)
    xlinc(other._annotinc[1]);
  if (other._have_zstart)
    zstart(other._zstart);
  if (other._have_zinc)
    zinc(other._zinc);
  if (other._have_corners)
    corners(other._corners);
  return *this;
}

ZgyWriterArgs&
ZgyWriterArgs::compressor(const std::string& name, const std::vector<std::string>& args)
{
  return compressor(InternalZGY::CompressFactoryImpl::getCompressor(name,args));
}

ZgyWriterArgs& ZgyWriterArgs::zfp_compressor(float snr)
{
  return compressor("ZFP", std::vector<std::string>{"snr", std::to_string(snr)});
}

ZgyWriterArgs& ZgyWriterArgs::lodcompressor(const std::string& name, const std::vector<std::string>& args)
{
  return lodcompressor(InternalZGY::CompressFactoryImpl::getCompressor(name,args));
}

ZgyWriterArgs&
ZgyWriterArgs::zfp_lodcompressor(float snr)
{
  return lodcompressor("ZFP", std::vector<std::string>{"snr", std::to_string(snr)});
}

ZgyWriterArgs&
ZgyWriterArgs::iocontext(const IOContext *value)
{
  _iocontext = value ? value->clone() : nullptr;
  return *this;
}

namespace Formatters {
  /**
   * \brief Return the string representation of the input enum type.
   */
  std::string enumToString(SampleDataType value)
  {
    switch (value) {
    case SampleDataType::unknown: return "SampleDataType::unknown";
    case SampleDataType::int8:    return "SampleDataType::int8";
    case SampleDataType::int16:   return "SampleDataType::int16";
    case SampleDataType::float32: return "SampleDataType::float32";
    default: return "SampleDataType::" + std::to_string((int)value);
    }
  }
  /**
   * \brief Return the string representation of the input enum type.
   */
  std::string enumToString(UnitDimension value)
  {
    switch (value) {
    case UnitDimension::unknown:  return "UnitDimension::unknown";
    case UnitDimension::time:     return "UnitDimension::time";
    case UnitDimension::length:   return "UnitDimension::length";
    case UnitDimension::arcangle: return "UnitDimension::arcangle";
    default: return "UnitDimension::" + std::to_string((int)value);
    }
  }
  /**
   * \brief Return the string representation of the input enum type.
   */
  std::string enumToString(DecimationType value)
  {
    switch (value) {
    case DecimationType::LowPass:          return "DecimationType::LowPass";
    case DecimationType::LowPassNew:       return "DecimationType::LowPassNew";
    case DecimationType::WeightedAverage:  return "DecimationType::WeightedAverage";
    case DecimationType::Average:          return "DecimationType::Average";
    case DecimationType::Median:           return "DecimationType::Median";
    case DecimationType::Minimum:          return "DecimationType::Minimum";
    case DecimationType::Maximum:          return "DecimationType::Maximum";
    //case DecimationType::MinMax:         return "DecimationType::MinMax";
    case DecimationType::Decimate:         return "DecimationType::Decimate";
    case DecimationType::DecimateSkipNaN:  return "DecimationType::DecimateSkipNaN";
    //case DecimationType::DecimateRandom: return "DecimationType::DecimateRandom";
    case DecimationType::AllZero:          return "DecimationType::AllZero";
    //case DecimationType::WhiteNoise:     return "DecimationType::WhiteNoise";
    case DecimationType::MostFrequent:     return "DecimationType::MostFrequent";
    case DecimationType::MostFrequentNon0: return "DecimationType::MostFrequentNon0";
    case DecimationType::AverageNon0:      return "DecimationType::AverageNon0";
    default: return "DecimationType::" + std::to_string((int)value);
    }
  }
  /**
   * \brief Return the string representation of the input enum type.
   */
  std::string enumToString(FinalizeAction value)
  {
    switch (value) {
    case FinalizeAction::Delete:           return "FinalizeAction::Delete";
    case FinalizeAction::Keep:             return "FinalizeAction::Keep";
    case FinalizeAction::BuildIncremental: return "FinalizeAction::BuildIncremental";
    case FinalizeAction::BuildFull:        return "FinalizeAction::BuildFull";
    case FinalizeAction::BuildNoHistogram: return "FinalizeAction::BuildNoHistogram";
    default: return "FinalizeAction::" + std::to_string((int)value);
    }
  }
  /**
   * \brief Output the string representation of the input enum type.
   */
  std::ostream& operator<<(std::ostream& os, SampleDataType value) {
    return os << enumToString(value);
  }
  /**
   * \brief Output the string representation of the input enum type.
   */
  std::ostream& operator<<(std::ostream& os, UnitDimension value) {
    return os << enumToString(value);
  }
  /**
   * \brief Output the string representation of the input enum type.
   */
  std::ostream& operator<<(std::ostream& os, DecimationType value) {
    return os << enumToString(value);
  }
  /**
   * \brief Output the string representation of the input enum type.
   */
  std::ostream& operator<<(std::ostream& os, FinalizeAction value) {
    return os << enumToString(value);
  }
} // namespace Formatters
} // namespace OpenZGY

#ifdef _MSC_VER
#pragma warning(pop)
#endif
