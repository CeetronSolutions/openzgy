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
#include "exception.h"
#include "safewriter.h"
#include "impl/enum.h"
#include "impl/file.h"
#include "impl/meta.h"
#include "impl/bulk.h"
#include "impl/databuffer.h"
#include "impl/transform.h"
#include "impl/lodalgo.h"
#include "impl/statisticdata.h"
#include "impl/histogramdata.h"
#include "impl/genlod.h"
#include "impl/compression.h"
#include "impl/guid.h"
#include "impl/logger.h"
#include "impl/environment.h"

#include <tuple>
#include <list>
#include <sstream>

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

namespace OpenZGY {
  class IOContext;
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
    logger_(1, "ERROR: " + msg);
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

/** \cond IMPL */ // Doxygen needs api.cpp, but not all of it.

/**
 * \brief convert between internal and external data types.
 *
 * \details Thread safety:
 * Thread safe because the class only contains static methods with no state.
 *
 * \internal
 * TODO-Test: Move to a separate file to allow it to be unit tested.
 */
class EnumMapper
{
  EnumMapper() = delete;
  EnumMapper(const EnumMapper&) = delete;
  EnumMapper& operator=(const EnumMapper&) = delete;
public:
  static SampleDataType mapRawDataTypeToSampleDataType(InternalZGY::RawDataType);
  static UnitDimension  mapRawHorizontalDimensionToUnitDimension(InternalZGY::RawHorizontalDimension);
  static UnitDimension  mapRawVerticalDimensionToUnitDimension(InternalZGY::RawVerticalDimension);
  static InternalZGY::RawDataType            mapSampleDataTypeToRawDataType(SampleDataType);
  static InternalZGY::RawHorizontalDimension mapUnitDimensionToRawHorizontalDimension(UnitDimension);
  static InternalZGY::RawVerticalDimension   mapUnitDimensionToRawVerticalDimension(UnitDimension);
  static InternalZGY::ZgyInternalWriterArgs mapWriterArgs(const ZgyWriterArgs&);
  static InternalZGY::LodAlgorithm mapDecimationTypeToLodAlgorithm(DecimationType value);
  static std::vector<InternalZGY::LodAlgorithm> mapDecimationTypeToLodAlgorithm(const std::vector<DecimationType>& values);
};

/** \endcond */

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
    // TODO-Low: Performance: Cache usable nlods.
    // Note that a file currently being written to but not finalized yet
    // will see nlods() == 1 (which is entirely correct) but to avoid
    // some chicken and egg problems it will still be permitted to read
    // back higher LOD levels. See ZgyInternalBulk::_validateUserPosition.
    //return _meta->ih().nlods();
    return InternalZGY::LookupTable::usableBrickLOD
      (this->_meta->ih().lodsizes(),
       this->_meta->ih().brickoffsets(),
       this->_meta->blup().lup(),
       this->_meta->blup().lupend());
  }

  // Currently not needed by any client.
  //virtual std::string dataid() const override
  //{
  //  return InternalZGY::GUID(_meta->ih().dataid()).toString();
  //}

  std::string verid() const override
  {
    return InternalZGY::GUID(_meta->ih().verid()).toString();
  }

  // Currently not needed by any client.
  //virtual std::string previd() const override
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
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(iocontext);
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
    _accessor.reset(new InternalZGY::ZgyInternalBulk(_fd, _meta, nullptr, false, _logger));
  }

  ~ZgyReader()
  {
    try {
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
      _logger(0, "ERROR closing a file opened for read: " +
              std::string(ex.what() ? ex.what() : "(null)"));
    }
  }

  /**
   * \brief Read an arbitrary region.
   *
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
  mutable bool _dirty; // if true need LOD, stats, histogram.
  compressor_t _compressor;
  compressor_t _lodcompressor;
  std::shared_ptr<LockedInProcess> _locked_in_process;
  LoggerFn _logger;

public:
  struct OpenForUpdate{};

  /**
   * \copydoc IZgyWriter::open()
   */
  explicit ZgyWriter(const ZgyWriterArgs& args)
    : _fd()
    , _accessor_rw()
    , _dirty(false)
    , _compressor(args._compressor)
    , _lodcompressor(args._lodcompressor ? args._lodcompressor : args._compressor)
  {
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(args._iocontext.get());
    _locked_in_process.reset(new LockedInProcess(args._filename, true, _logger));

    for (int ii=0; ii<3; ++ii) {
      if (args._size[ii] < 1)
        throw Errors::ZgyUserError("Survey size cannot be empty or negative.");
      else if (args._size[ii] >= std::numeric_limits<std::int32_t>::max())
        throw Errors::ZgyUserError("Survey size is too large.");
      }
    const bool compress = args._compressor || args._lodcompressor;
    // This is both pointless and expensive to support.
    if (compress && args._datatype != SampleDataType::float32)
      throw Errors::ZgyUserError("Compressed files need to be stored as float.");

    // Set the protected metadata information in the ZgyMeta base class.
    // Also store a mutable copy in the current class.
    // ZgyInternalMeta does not contain any file descriptor. This means
    // we need to hold on to _fd ourselves and provide it when it is time
    // to flush metadata to disk.
    InternalZGY::ZgyInternalWriterArgs iargs = EnumMapper::mapWriterArgs(args);
    _meta_rw.reset(new InternalZGY::ZgyInternalMeta(iargs, compress, _logger));
    _meta = _meta_rw; // in base class

    // The file creation was deferred until after the consistency checks.
    _fd = InternalZGY::FileFactory::instance().create(
        args._filename, InternalZGY::OpenMode::Truncate, ctxt.get());
    _meta_rw->flushMeta(_fd);
    // Compress or not at this level only controls alignment etc.
    // The actual compression functor is passed to each write.
    _accessor_rw.reset(new InternalZGY::ZgyInternalBulk(_fd, _meta_rw, _meta_rw, compress, _logger));

    // A file that is recently created is dirty because even if nothing
    // was written to it the low resolution bricks should be generated.
    // Otherwise the file will be flagged as incomplete, nlods()==1.
    // If a file is opened for update (not currently supported) *and*
    // the low resolution data is ok then the dirty flag should not be
    // set until the first write. Not a big deal because why open a file
    // for write without writing to it?
    // And after we get proper update support with incremental finalize
    // it shouldn't matter anyway.
    _dirty = true;
  }

  /**
   * \brief Open an existing file for update.
   * \details
   * See the ZgyReader constructor and the "truncate" ZgyWriter constructor.
   * Understandably this nethod contains elements from both.
   */
  ZgyWriter(const ZgyWriterArgs& args, OpenForUpdate)
    : _fd()
    , _accessor_rw()
    , _dirty(false)
    , _compressor(args._compressor)
    , _lodcompressor(args._lodcompressor ? args._lodcompressor : args._compressor)
  {
    std::shared_ptr<IOContext> ctxt;
    std::tie(_logger, ctxt) = setupLogging(args._iocontext.get());
    _logger(1, "open for update: \"" + args._filename + "\"");
    _locked_in_process.reset(new LockedInProcess(args._filename, true, _logger));

    const bool compress = _compressor || _lodcompressor;

    _fd = InternalZGY::FileFactory::instance().create(
        args._filename, InternalZGY::OpenMode::ReadWrite, ctxt.get());

    // Set the protected metadata information in the ZgyMeta base class.
    // The metadata ctor is very different on update. bulkdata is not.
    // ZgyInternalMeta does not retain the file descriptor. The file is
    // only used inside the constructor to populate the headers.
    InternalZGY::ZgyInternalWriterArgs iargs = EnumMapper::mapWriterArgs(args);
    _meta_rw.reset(new InternalZGY::ZgyInternalMeta(_fd, iargs, compress, _logger));
    _meta = _meta_rw;

    if (compress && _meta->ih().datatype() != InternalZGY::RawDataType::Float32)
      throw Errors::ZgyUserError("Compressed files need to be stored as float.");

    // At the implementation level the bulk- and meta access are separate,
    // and the bulk accessor needs some of the meta information to work.
    // The accessor will have its own metadata pointer which it holds on to.
    // It also holds on to the file descriptor. Unlike the metadata reader
    // it is obviously not possible to read everything up front.
    _accessor_rw.reset(new InternalZGY::ZgyInternalBulk(_fd, _meta_rw, _meta_rw, compress, _logger));

    // If low resolution data, statistics, or histogram are not
    // available, mark as dirty because even if there has been no
    // change to the bulk data the code should by default try to
    // create lods and statistics on close.
    if (!_meta->has_finalized_feature())
      this->_dirty = true;

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

  /**
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
    this->_dirty = true;
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
    this->_dirty = true;
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
    this->_dirty = true;
    buffer.reset();
    if (fakeshared.use_count() != 1) // Actually a fatal error.
      throw Errors::ZgyInternalError("A Reference to the user's buffer was retained.");
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Write all-constant data.
   *
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
    this->_dirty = true;
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
    this->_dirty = true;
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
    this->_dirty = true;
    ZombieCheck::checkNonUniquePtr(accessor_rw, "write accessor");
  }

  /**
   * \brief Remove the presumably stale statistics, histogram, and lowres.
   */
  void _delete_derived()
  {
    // Called explicitly from user code, meaning the application
    // doesn't want the expensive finalize step. Possibly because
    // additional data will be written later (new feature). Or
    // because the consumer isn't interested in lowres data.
    //
    // If update is allowed then a few extra steps will be needed.
    // If this is the first write then these steps are no-ops

    // 1. Clear the statistics because it is out of date.
    constexpr double inf{std::numeric_limits<double>::infinity()};
    this->_meta_rw->ih().setstats(0,0,0,inf,-inf);

    // 2. Clear the histogram, except that we have already started
    //    deciding the histogram range to use as the range of all
    //    samples written so far. Store that range from
    //    _written_sample_min/_max so the next write can pick
    //    up that range from where we left off.
    //    For float cubes this is useful. For integral cubes,
    //    not so much. The real histogram range currently ends up
    //    as the possible range, not the codingrange. In case
    //    that changes I'll save the actual range anyway. Including
    //    converting from storage (used to track the range)
    //    to float (expected in the histogram range).
    const std::array<double,2> range=this->_accessor_rw->valueRangeWritten();
    const InternalZGY::IInfoHeaderAccess& ih = this->_meta_rw->ih();
    const double slope     = ih.storagetofloat_slope();
    const double intercept = ih.storagetofloat_intercept();
    this->_meta_rw->hh().sethisto(range[0] * slope + intercept,
                                  range[1] * slope + intercept,
                                  nullptr, 0);

    // 3. Logically erase all low resolution bricks.
    std::vector<std::uint64_t>& blup = _meta_rw->blup().lup();
    std::vector<std::uint64_t>& bend = _meta_rw->blup().lupend();
    //auto code1 = blup[0];
    //auto nl1 = InternalZGY::LookupTable::usableBrickLOD
    //  (ih.lodsizes(), ih.brickoffsets(), blup, bend);
    InternalZGY::LookupTable::setNoBrickLOD
      (ih.lodsizes(), ih.brickoffsets(), &blup, &bend);
    //auto code2 = blup[0];
    //auto nl2 = InternalZGY::LookupTable::usableBrickLOD
    //  (ih.lodsizes(), ih.brickoffsets(), blup, bend);
    //std::cout << "Was nlods " << nl1 << " code " << code1
    //          << " is nlods " << nl2 << " code " << code2 << std::endl;

    // 4. Change the version number to 4 is now done in _close_internal().
  }

  /**
   * \brief Remove the presumably stale statistics, histogram, and lowres.
   */
  void _create_derived(
       const std::vector<DecimationType>& decimation,
       const std::function<bool(std::int64_t,std::int64_t)>& progress,
       FinalizeAction action,
       ZgyWriter::LoggerFn logger)
  {
    InternalZGY::DataBuffer::cleartimers(true);
    this->_dirty = false;
    std::shared_ptr<InternalZGY::StatisticData> stats;
    std::shared_ptr<InternalZGY::HistogramData> histo;

    if (action == FinalizeAction::BuildIncremental) {
      std::tie(stats, histo) = _accessor_rw->trackedChanges();
      if (!stats || !histo) {
        action = FinalizeAction::BuildFull;
        logger(2, "finalize: switched to full because no tracked changes.");
      }
      else if (histo->getcount() <= 0) {
        // Obscure case: The data has changed so much that none of the
        // existing samples, not even defaultvalue, now fit inside the
        // histogram. Be nice to applications (and unit tests) and
        // just fall back to a full rebuild. See test_reopen() in
        // test_reopen.cpp.
        action = FinalizeAction::BuildFull;
        logger(2, "finalize: switched to full because of empty histogram.");
      }
      else {
        logger(2, "finalize: will attempt incremental build.");
      }
    }
    else {
      logger(2, "finalize: will do a full build.");
    }

    // Note: this->_accessor._metadata == this->_meta but the former is
    // private and only used to allow ZgyInternalBulk to use some parts of
    // the metadata. This is why both _accessor and _meta are needed below.
    std::tie(stats, histo) =
      InternalZGY::GenLodC(_accessor_rw, _meta_rw, _lodcompressor,
                           EnumMapper::mapDecimationTypeToLodAlgorithm(decimation),
                           action == FinalizeAction::BuildIncremental,
                           action == FinalizeAction::BuildNoHistogram,
                           progress, _logger).call();

    // If doing an incremental rebuild: The statistics collected by
    // GenLodC are bogus because we don't know which bricks were
    // tallied. Those in _modified_bricks should have been but there
    // might be more. Fortunately, if the file was reopened then the
    // histogram range is already fixed. It is safe to tally changes
    // to statistics and histogram as they occur.
    if (action == FinalizeAction::BuildIncremental) {
      //std::cout << "Scanned changes " << stats->toString()
      //          << " histo " << histo->toString() << std::endl;
      std::tie(stats, histo) = _accessor_rw->trackedChanges();
      //std::cout << "Tracked changes " << stats->toString()
      //          << " histo " << histo->toString() << std::endl;
    }

    _accessor_rw->trackedBricksShowDirty(/*loglevel=*/6);

    if (action == FinalizeAction::BuildNoHistogram || !stats || !histo) {
      // GenLod couldn't collect stats, or it did but we asked it not to.
      logger(4, "As requested, do not store statistics or histogram");
      return;
    }

    // Statistics and histogram are accumulated in storage units
    // and we want to store user units on the file.
    // Technically I could scale in place since stats and histo
    // shouldn't be used after this so it is ok to clobber them.
    // But I will copy just in case.
    InternalZGY::StatisticData scaled_stats(*stats);
    InternalZGY::HistogramData scaled_histo(*histo);
    const std::array<double,2> factors = this->_meta->ih().storagetofloat();
    scaled_stats.scale(0, 1, factors[1], factors[0] + factors[1]);
    scaled_histo.scale(0, 1, factors[1], factors[0] + factors[1]);
    // TODO-Low: if using a 64k temporary histogram. scaled_histo.resize(256);
    // A temp histogram would probably mess up the incremental tracking
    // of histogram information, so take care. Might need to simply
    // disable the histogram for files that were opened for update.
    // TODO-Low: Some kind of static_assert making the following cast safe.
    // Note reason for cast: ZGY only deals with int8, int16, float
    // which means that converting the value range to float will not
    // lose any precision. If we later decide to support int32 then
    // this and many similar places need to be updated.
    this->_meta_rw->ih().setstats
      (scaled_stats.getcnt(), scaled_stats.getsum(), scaled_stats.getssq(),
       scaled_stats.getmin(), scaled_stats.getmax());
    this->_meta_rw->hh().sethisto(scaled_histo.getmin(), scaled_histo.getmax(),
                                  scaled_histo.getbins(), scaled_histo.getsize());
    if (logger(4, "")) {
      std::stringstream ss;
      ss << "Histogram Limits "
         << histo->getmin() << " " << histo->getmax()
         << " Converted "
         << scaled_histo.getmin() << " " << scaled_histo.getmax();
      logger(4, ss.str());
    }
  }

  /**
   * \brief Maybe generate low resolution data, statistics, and histogram.
   *
   * This method will be called automatically from close(), but in
   * that case it is not possible to request a progress callback.
   *
   * It is valid to call finalize() and then continue writing to the
   * file. This is currently not useful. The reason the application
   * might want this is to allow reading low resolution data from
   * a file that is still open for write. The API blocks this today,
   * simply by removing the lod parameter from IZgyWriter::read().
   *
   * If the processing raises an exception the data is still marked
   * as clean. So a second attempt will do nothing unless the
   * caller passes force=True.
   *
   * The C++ code is very different from Python because it needs
   * an entirely different approach to be performant.
   *
   * \param decimation: Optionally override the decimation algorithms by
   *                    passing an array of DecimationType
   *                    with one entry for each level of detail. If the
   *                    array is too short then the last entry is used for
   *                    subsequent levels.
   *
   * \param progress:   Function(done, total) called to report progress.
   *                    If it returns False the computation is aborted.
   *                    The callback implementation should not invoke any
   *                    OpenZGY methods except for any calls expressly
   *                    documented as safe for this purpose and deadlock-free.
   *                    Will be called at least one, even if there is no work.
   *
   * \param action:     Whether to compute decimated bricks, statistics,
   *                    and histogram.
   *                    - Delete:
   *                      Remove any existing information. If force=true
   *                      this is unconditinal. Otherwise it changes to
   *                      "Keep" if the information already exists
   *                      and is not stale. With force=off, Delete and Keep
   *                      actually do the same thing.
   *                    - Keep:
   *                      Keep any existing information. If force=true this
   *                      will keep the old information even if it is looks
   *                      stale. Otherwise the mode changes to "Delete"
   *                      if OpenZGY cannot guarantee that the information
   *                      is still correct.
   *                    - BuildIncremental:
   *                      Do the minimum amount of work needed to bring the
   *                      derived data up to date. May cause the statistics
   *                      and histogram to be less accurate. Changes to "Keep"
   *                      if the information is already up to date. changes
   *                      to "BuildFull" if information does not exist already
   *                      or if the file is compressed or if incremental
   *                      cannot work for other reasons. Setting force=true
   *                      might do a little more work, to be defined.
   *                    - BuildFull:
   *                      Delete and re-create derived information from
   *                      scratch. On the cloud this might waste disk space.
   *                      If force=true the information is re-calculated even
   *                      if already good. That may be useful if you want to
   *                      change the decimation algorithm. Otherwise change
   *                      to "Keep" if the information is already up to date.
   *                      Never changes to BuildIncremental.
   *                    - BuildNoHistogram:
   *                      As BuildFull, but do not collect or store histogram
   *                      and statistics. This might be useful to speed up
   *                      writes that we know will only be read by Petrel.
   *                      Or other apps that don't strictly need them.
   *                      CAVEAT: BuildNoHistogram is incompatible with
   *                      DecimationType::WeightedAverage. So, indirectly
   *                      this option will lower the quality of LoD data,
   *
   * \param force:      If true, use the "action" parameter exactly as given.
   *                    Use with caution.
   */
  void finalize(const std::vector<DecimationType>& decimation_in,
                const std::function<bool(std::int64_t,std::int64_t)>& progress,
                FinalizeAction action = FinalizeAction::BuildDefault,
                bool force = false) override
  {
    std::vector<DecimationType> decimation(decimation_in);
#if 1
    // UGLY KLUDGE for testing a questionable feature. If possible, set
    // this programmatically instead. Or, better still, don't use it at all.
    const std::string kludge("OPENZGY_KLUDGE_WRITE_NO_HISTOGRAM");
    if (InternalZGY::Environment::getNumericEnv(kludge.c_str(), -1) > 0) {
      if (action  != FinalizeAction::BuildDefault) {
        _logger(0, "Ignoring " + kludge + " because action is not BuildFull");
      }
      else {
        action = FinalizeAction::BuildNoHistogram;
        if (decimation.empty()) {
          decimation = std::vector<DecimationType>
            {DecimationType::LowPass, DecimationType::Average};
        }
        else {
          for (std::size_t ii=0; ii<decimation.size(); ++ii)
            if (decimation[ii] == DecimationType::WeightedAverage)
              decimation[ii] = DecimationType::Average;
        }
        _logger(0, "DANGER: " + kludge + " is enabled.");
      }
    }
#endif // UGLY KLUDGE END

    // Called with error flag set is actually valid here.
    //throw_if_not_writable();
    if (!_fd || !_accessor_rw || !_meta_rw)
      throw Errors::ZgyUserError("ZGY file is not open for write.");

    if (!force) {
      if (!this->_dirty)
        action = FinalizeAction::Keep;
      else if (action == FinalizeAction::Keep)
        action = FinalizeAction::Delete;
    }

    if (this->errorflag() &&
        (action == FinalizeAction::BuildIncremental ||
         action == FinalizeAction::BuildNoHistogram ||
         action == FinalizeAction::BuildFull))
    {
      // Don't do anything that requires more bulk data to be written.
      // Attempting to update the headers, to clear exisisting data,
      // is ok to try.
      action = FinalizeAction::Delete;
    }

    switch (action) {

    case FinalizeAction::Delete:
      _delete_derived();
      if (progress)
        progress(0, 0);
      // Next finalize must be a full one, so stop tracking.
      _accessor_rw->trackedBricksTryEnable(false);
      break;

    case FinalizeAction::Keep:
      if (progress)
        progress(0, 0);
      break;

    case FinalizeAction::BuildIncremental:
    case FinalizeAction::BuildNoHistogram:
    case FinalizeAction::BuildFull:
      try {
        _create_derived(decimation, progress, action, _logger);
        // Reset and start over tracking changes.
        _accessor_rw->trackedBricksTryEnable(true);
      }
      catch(...) {
        // Next finalize must be a full one.
        _accessor_rw->trackedBricksTryEnable(false);
        throw;
      }
      break;
    }

    // If the application makes an explicit call to finalize() to
    // suppress the low resolution generation then it shouldn't need
    // to also call close_incomplete() instead of close() to prevent
    // the lowres from be computed anyway.

    // So: Always mark the file as clean at this point even when it
    // provably isn't. Lowres can be missing or even stale.
    // Nothing will be done about that. A finalize with force=true
    // will be needed if the application changes its mind.

    this->_dirty = false;
  }

  /**
   * \brief Flush the file to disk and close it.
   *
   * This version of close() will not calculate statistics and low
   * resolution bricks if thay are missing. It will delete them if
   * they exist and are stale.
   *
   * An incomplete file is still usable for reading full resolution
   * data as long as it is accessed using OpenZGY and not the old
   * ZGY-Public. OpenZGY reports nlods==1 for incomplete files
   * so applications can deal with it. ZGY-Public will not be allowed
   * to open incomplete files. Version is set >3 to prevent this.
   */
  void _close_internal()
  {
    if (!this->_fd || !this->_accessor_rw || !this->_meta_rw) {
      _locked_in_process.reset(); // Just in case.
      return; // looks like _close_internal alrady called.
    }

    // Prevent ZGY-Public from opening the file if appropriate.
    this->_meta_rw->fh().set_version(_meta_rw->can_old_library_read() ? 3 : 4);

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
    InternalZGY::DataBuffer::cleartimers(true);

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
   * This version of close() will not calculate statistics and low
   * resolution bricks if thay are missing. It will delete them if
   * they exist and are stale.
   *
   * Calling this method is equivalent to calling finalize() with
   * the "Keep" action and then calling the regular close().
   *
   * An incomplete file is still usable for reading full resolution
   * data as long as it is accessed using OpenZGY and not the old
   * ZGY-Public. OpenZGY reports nlods==1 for incomplete files
   * so applications can deal with it. ZGY-Public will not be allowed
   * to open incomplete files. Version is set >3 to prevent this.
   */
  void close_incomplete() override
  {
    if (_fd) {
      // Delete instead of keep if the derived information is stale.
      // Always a no-op if an explicit call was made to finalize().
      finalize(std::vector<DecimationType>{}, nullptr,
               FinalizeAction::Keep, false);
      _close_internal();
    }
  }

  /**
   * \brief Flush the file to disk and close it.
   *
   * If the file has been written to, the application is encouraged to
   * call finalize() before close(). This gives more control over the
   * process and allows using a progress callback to track generation
   * of low resolution data.
   *
   * If the application fails to call finalize() it will be done
   * automatically here. If the application really wants to skip
   * the finalize then call close_incomplete() instead.
   * Or call finalize() with the "Keep" action.
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
    if (_fd) {
      finalize(std::vector<DecimationType>
               {DecimationType::LowPass, DecimationType::WeightedAverage},
               nullptr,
               FinalizeAction::BuildDefault,
               false);
      _close_internal();
    }
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

  void deletefile(const std::string& filename, bool missing_ok)
  {
    _fd->deleteFile(filename, missing_ok);
  }

  std::string alturl(const std::string& filename)
  {
    return _fd->altUrl(filename);
  }

  std::string idtoken()
  {
    return _fd->idToken();
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

/** \cond IMPL */ // Doxygen needs api.cpp, but not all of it.

/**
 * Map between enums used in the public API and the internal ones that
 * might change without notice and might be used to define the actual
 * numbers written to file.
 *
 * As a general rule, if an invalid enum tag is encountered while
 * converting from public to internal then this will throw an exceptiom
 * because it would be a user error. When converting from internal to
 * public the error is probably bad data encountered on the file.
 * The mapping function might just return "unknown", leaving to the
 * caller to decide whether this should be silently ignored.
 */
SampleDataType
EnumMapper::mapRawDataTypeToSampleDataType(InternalZGY::RawDataType value)
{
  using InternalZGY::RawDataType;
  switch (value) {
  case RawDataType::SignedInt8:    return SampleDataType::int8;
  case RawDataType::UnsignedInt8:  return SampleDataType::unknown;
  case RawDataType::SignedInt16:   return SampleDataType::int16;
  case RawDataType::UnsignedInt16: return SampleDataType::unknown;
  case RawDataType::SignedInt32:   return SampleDataType::unknown;
  case RawDataType::UnsignedInt32: return SampleDataType::unknown;
  case RawDataType::Float32:       return SampleDataType::float32;
  case RawDataType::IbmFloat32:    return SampleDataType::unknown;
  default:                         return SampleDataType::unknown;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
InternalZGY::RawDataType
EnumMapper::mapSampleDataTypeToRawDataType(SampleDataType value)
{
  using InternalZGY::RawDataType;
  switch (value) {
  case SampleDataType::int8:    return RawDataType::SignedInt8;
  case SampleDataType::int16:   return RawDataType::SignedInt16;
  case SampleDataType::float32: return RawDataType::Float32;
  case SampleDataType::unknown:
    throw Errors::ZgyUserError("SampleDataType::unknown is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for SampleDataType");
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
UnitDimension
EnumMapper::mapRawHorizontalDimensionToUnitDimension(InternalZGY::RawHorizontalDimension value)
{
  using InternalZGY::RawHorizontalDimension;
  switch (value) {
  default:
  case RawHorizontalDimension::Unknown:  return UnitDimension::unknown;
  case RawHorizontalDimension::Length:   return UnitDimension::length;
  case RawHorizontalDimension::ArcAngle: return UnitDimension::arcangle;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * The distinction between TWT and OWT, which is stored on the file, is lost.
 * To my knowledge there is no code that recognizes OWD in ZGY files anyway,
 * so it is better to not expose this to the API.
 */
UnitDimension
EnumMapper::mapRawVerticalDimensionToUnitDimension(InternalZGY::RawVerticalDimension value)
{
  using InternalZGY::RawVerticalDimension;
  switch (value) {
  default:
  case RawVerticalDimension::Unknown:    return UnitDimension::unknown;
  case RawVerticalDimension::Depth:      return UnitDimension::length;
  case RawVerticalDimension::SeismicTWT: return UnitDimension::time;
  case RawVerticalDimension::SeismicOWT: return UnitDimension::time;
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * An explicit "unknown" is allowed, but passing e.g. "time" in this
 * context (which is supposedly horozontal) will raise an exception.
 */
InternalZGY::RawHorizontalDimension
EnumMapper::mapUnitDimensionToRawHorizontalDimension(UnitDimension value)
{
  using InternalZGY::RawHorizontalDimension;
  switch (value) {
  case UnitDimension::unknown:  return RawHorizontalDimension::Unknown;
  case UnitDimension::length:   return RawHorizontalDimension::Length;
  case UnitDimension::arcangle: return RawHorizontalDimension::ArcAngle;
  case UnitDimension::time:
    throw Errors::ZgyUserError("UnitDimension::time is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for UnitDimension");
  }
};

/**
 * See mapRawDataTypeToSampleType() for comments.
 *
 * An explicit "unknown" is allowed, but passing e.g. "arcangle" in this
 * context (which is supposedly vertical) will raise an exception.
 */
InternalZGY::RawVerticalDimension
EnumMapper::mapUnitDimensionToRawVerticalDimension(UnitDimension value)
{
  using InternalZGY::RawVerticalDimension;
  switch (value) {
  case UnitDimension::unknown: return RawVerticalDimension::Unknown;
  case UnitDimension::time:    return RawVerticalDimension::SeismicTWT;
  case UnitDimension::length:  return RawVerticalDimension::Depth;
  case UnitDimension::arcangle:
    throw Errors::ZgyUserError("UnitDimension::arcangle is not allowed here.");
  default:
    throw Errors::ZgyUserError("Invalid enum tag for UnitDimension");
  }
}

InternalZGY::ZgyInternalWriterArgs
EnumMapper::mapWriterArgs(const ZgyWriterArgs& in)
{
  InternalZGY::ZgyInternalWriterArgs result;
  result.filename    = in._filename;
  result.size        = in._size;
  result.bricksize   = in._bricksize;
  result.datatype    = mapSampleDataTypeToRawDataType(in._datatype);
  result.datarange   = in._datarange;
  result.zunitdim    = mapUnitDimensionToRawVerticalDimension(in._zunitdim);
  result.hunitdim    = mapUnitDimensionToRawHorizontalDimension(in._hunitdim);
  result.zunitname   = in._zunitname;
  result.hunitname   = in._hunitname;
  result.zunitfactor = in._zunitfactor;
  result.hunitfactor = in._hunitfactor;
  result.zstart      = in._zstart;
  result.zinc        = in._zinc;
  result.annotstart  = in._annotstart;
  result.annotinc    = in._annotinc;
  result.corners     = in._corners;
  result.have_size      = in._have_size;
  result.have_bricksize = in._have_bricksize;
  result.have_datatype  = in._have_datatype;
  result.have_datarange = in._have_datarange;
  result.have_zunit     = in._have_zunit;
  result.have_hunit     = in._have_hunit;
  result.have_ilstart   = in._have_ilstart;
  result.have_ilinc     = in._have_ilinc;
  result.have_xlstart   = in._have_xlstart;
  result.have_xlinc     = in._have_xlinc;
  result.have_zstart    = in._have_zstart;
  result.have_zinc      = in._have_zinc;
  result.have_corners   = in._have_corners;
  return result;
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
InternalZGY::LodAlgorithm
EnumMapper::mapDecimationTypeToLodAlgorithm(DecimationType value)
{
  using InternalZGY::LodAlgorithm;
  switch (value) {
  case DecimationType::LowPass:          return LodAlgorithm::LowPass;
  case DecimationType::WeightedAverage:  return LodAlgorithm::WeightedAverage;
  case DecimationType::Average:          return LodAlgorithm::Average;
  case DecimationType::Median:           return LodAlgorithm::Median;
  case DecimationType::Minimum:          return LodAlgorithm::Minimum;
  case DecimationType::Maximum:          return LodAlgorithm::Maximum;
  //case DecimationType::MinMax:         return LodAlgorithm::MinMax;
  case DecimationType::Decimate:         return LodAlgorithm::Decimate;
  case DecimationType::DecimateSkipNaN:  return LodAlgorithm::DecimateSkipNaN;
  //case DecimationType::DecimateRandom: return LodAlgorithm::DecimateRandom;
  case DecimationType::AllZero:          return LodAlgorithm::AllZero;
  //case DecimationType::WhiteNoise:     return LodAlgorithm::WhiteNoise;
  case DecimationType::MostFrequent:     return LodAlgorithm::MostFrequent;
  case DecimationType::MostFrequentNon0: return LodAlgorithm::MostFrequentNon0;
  case DecimationType::AverageNon0:      return LodAlgorithm::AverageNon0;
  default:
    throw Errors::ZgyUserError("Invalid enum tag for DecimationType");
  }
}

/**
 * See mapRawDataTypeToSampleType() for comments.
 */
std::vector<InternalZGY::LodAlgorithm>
EnumMapper::mapDecimationTypeToLodAlgorithm(const std::vector<DecimationType>& values)
{
  using InternalZGY::LodAlgorithm;
  if (!values.size())
    return std::vector<InternalZGY::LodAlgorithm>
      {LodAlgorithm::LowPass, LodAlgorithm::WeightedAverage};
  std::vector<InternalZGY::LodAlgorithm> result;
  for (DecimationType value : values)
    result.push_back(mapDecimationTypeToLodAlgorithm(value));
  return result;
}

/** \endcond */

} // namespace OpenZGY::Impl
// Remain in namespace OpenZGY

ProgressWithDots::ProgressWithDots(int length, std::ostream& outfile)
  : _dots_printed(0)
  , _length(length)
  , _outfile(outfile)
  , _mutex()
{
}

bool
ProgressWithDots::operator()(std::int64_t done, std::int64_t total)
{
  if (_length < 1)
    return true;
  std::lock_guard<std::mutex> lk(_mutex);
  if (_dots_printed == 0) {
    _outfile << "[" + std::string(_length, ' ') + "]\r[" << std::flush;
  }
  std::int64_t needed = (total <= 0) ? 1 : 1 + ((done * (_length-1)) / total);
  if (needed > _dots_printed) {
    while (needed > _dots_printed) {
      _outfile << '.';
      _dots_printed += 1;
    }
    _outfile << std::flush;
  }
  if (done == total)
    _outfile << "\n" << std::flush;
  return true;
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
  return std::shared_ptr<IZgyReader>(new Impl::ZgyReader(filename, iocontext));
}

std::shared_ptr<IZgyWriter>
IZgyWriter::open(const ZgyWriterArgs& args)
{
  auto unsafe = std::shared_ptr<IZgyWriter>(new Impl::ZgyWriter(args));
  return std::shared_ptr<IZgyWriter>(new Impl::ZgySafeWriter(unsafe));
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
IZgyWriter::reopen(const ZgyWriterArgs& args)
{
  auto unsafe = std::shared_ptr<IZgyWriter>(new Impl::ZgyWriter(args, Impl::ZgyWriter::OpenForUpdate{}));
  return std::shared_ptr<IZgyWriter>(new Impl::ZgySafeWriter(unsafe));
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
