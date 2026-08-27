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

/** file: impl/bulk.cpp
 *
 * This file contains a number of classes and free functions dealing with
 * bulk data access. There is a similar file impl.meta that deals with
 * metadata. Nothing in this file should be directly visible to users of
 * the public API.
 *
 *     impl.bulk.ZgyInternalBulk:
 *
 *         * Internal methods to read and write bulk.
 *         * Needs access to much of the meta data information.
 *           Currently this is handled by giving this class references to
 *           the ZgyInternalMeta instance. This gives access to more
 *           metadata than the class actually needs.
 *         * Refactoring notes: There are a couple of sections in this
 *           class that ought to have been separated out. As it is now
 *           the class does too much.
 *
 *             - Write support
 *             - Lod generation
 *             - Histogram generation
 *
 *     Code flow for writing:
 *
 *     The main logic is a couple of methods where each function might have
 *     ended with calling the next level down. But to make multi threading
 *     easier to implement (and obfuscate the code a bit),
 *     prepareWriteOneBrick() and several others return
 *     to the calling writeAlignedBrickList() with an argument package that
 *     should be sent to the next method in the chain.
 *
 *     The order of calls is:
 *
 *     writeBricksInternal(list of bricks)
 *         Used only by genlod and the future IZgyWriter::writebricks().
 *         Input validation.
 *         Assumes caller is sending storage value type, does NOT convert.
 *
 *     writeRegion (arbitrary region, user's value_type)
 *         Used only from the general API: IZgyWriter::write() and writeconst()
 *         Input validation.
 *         Apply conversion float -> storage if needed.
 *         Apply splitting into bricks.
 *
 *     writeAlignedBrickList(list of bricks)
 *         From here to the end used by both the general and the brick API.
 *         ST:  Handle the reset-all shortcut.
 *         ST*: Call to read old brick contents (cheap if never written before).
 *         MT:  Update statistics, histogram, and dirty brick information.
 *         MT:  Copy data if read/modify/write.
 *         ST*: Incremental low resolution generation.
 *         ST#: Calls to to the actual write.
 *         Legend:
 *         MT = multi-threaded. ST = Single threaded.
 *         ST* = single call on list, callee might do multi-threading.
 *         ST# = single threaded loop, callee can buffer data and write MT.
 *
 *     readBricksToNewBuffer(list of bricks)
 *         Used for r/m/w if misaligned access.
 *         Used for incremental tracking of statistics and histogram.
 *
 *     prepareWriteOneBrick (single brick, native vt) [Called MT]
 *
 *         Might convert from a scalar to a regular buffer or vice versa.
 *         Might decide whether to update an existing brick, and where.
 *         Might decide to veto the compression.
 *         The next step depends on whether the buffer ended up scalar
 *         (writeOneConstantBrick) or not (prepareWriteOneNormalBrick
 *          eventually followed by writeOneNormalBrickWithRetry)
 *
 *     prepareWriteOneNormalBrick [Called MT]
 *
 *         Apply padding, compression, byte swapping, and convert from
 *         DataBuffer to void* + size suitable for the file i/o layer.
 *
 *     writeOneConstantBrick
 *         Apply conversion of the constant value to a lookup table entry,
 *         Pass on the write request to the function that updates the lut.
 *
 *     writeOneRealBrickWithRetry
 *
 *        Allocate space on the bulk file as needed.
 *        Handle retries needed if the segment (on the cloud) was closed.
 *        Update the lookup table if a new block was allocated.
 *
 *     Code flow for reading (Brick API, from genlod and while writing)
 *
 *     readBricksToNewBuffers
 *        readConstValueBricks [ST]
 *           readConstantValue
 *        readBricksInternal
 *           readBricksToExistingBuffer
 *              logReadingSingleBrick
 *              startReadingSingleBrick
 *              actualReadingAllBricks
 *              createDeliverance
 *              -via functor-
 *              xx_readv
 *              deliverOneBrick
 *
 *     Code flow for reading (General API, i.e. read())
 *
 *        readToExistingBuffer
 *        (calls same as readBricksToExistingBuffer)
 *
 */

#include "enum.h"
#include "types.h"
#include "bulk.h"
#include "meta.h"
#include "file.h"
#include "databuffer.h"
#include "arrayops.h"
#include "structaccess.h"
#include "subtiling.h"
#include "logger.h"
#include "compression.h"
#include "environment.h"
#include "fancy_timers.h"
#include "workorder.h"
#include "histogrambuilder.h"
#include "expandablebuilder.h"
#include "tracktouched.h"
#include "statisticdata.h"
#include "histogramdata.h"
#include "../exception.h"
#include "genlod.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <sstream>
#include <atomic>
#include <omp.h>
#include <mutex>

namespace InternalZGY {
#if 0
}
#endif

using namespace InternalZGY::ArrayOps;
using namespace InternalZGY::Formatters;

namespace {
  /**
   * Return non-zero if there should be special case handling of the trivial
   * case of reading exactly one thread. The feature is on by default.
   */
  static bool
  expedited_read()
  {
    static int enable = Environment::getNumericEnv("OPENZGY_EXPEDITED_READ", 1);
    return enable > 0;
  }

  /**
   * Return non-zero if there should be special case handling of the trivial
   * case of writing exactly one brick.
   */
  static bool
  expedited_write()
  {
    static int enable = Environment::getNumericEnv("OPENZGY_EXPEDITED_WRITE", 1);
    return enable > 0;
  }

  /**
   * Strictly for debugging. Count the number of negative bin counts,
   * indicating a rounding error or an outright bug. The caller is
   * responsible for any required locking.
   */
  static std::int64_t
  negativeBins(const HistogramData& histo)
  {
    std::int64_t result{0};
    const std::int64_t* ptr = histo.getbins();
    const std::int64_t* end = ptr + histo.getsize();
    for (; ptr < end; ++ptr)
      if (*ptr < 0)
        result -= *ptr;
    return result;
  }

#if 0
  /**
   * Add to or subtract from the statistics and histogram provided
   *
   * Invokes ExpandableBuilder::staticAddOrSub() with some logging.
   * Caveat: Should not be used from inside a "parallel for".
   * Caveat: Not useful if updating a temporary variable.
   *
   * So, for old code (only) it can be invoked instead of staticAddOrSub()
   * in trackChangesT(). For new code (only) it can be invoked instead of
   * staticAddOrSub() near the end of writeAlignedBrickList() where the
   * "real" statistics get updated.
   */
  static void
  addChanges(
       StatisticData *stats,
       HistogramData *histo,
       const StatisticData& stats_change,
       const HistogramData& histo_change,
       bool add,
       // The rest is just for debug logging.
       const index3_t& size, // data->size3d()
       const index3_t& valid, // range.second
       bool new_code)
  {
    // Note that counting negative bins is somewhat useless if multi-threaded.
    // Which is why I don't bother locking, either.
    std::int64_t neg1 = negativeBins(*histo);
    std::cout << "track changes"
              << " size " << size << " valid " << valid
              << "\n"
              << (new_code ? "OLD" : "old")
              << " stats " << stats->toString()
              << " histo " << histo->toString(false)
              << "\n"
              << (new_code ? add ? "ADD" : "SUB" : add ? "add" : "sub")
              << " stats " << stats_change.toString()
              << " histo " << histo_change.toString(false)
              << "\n"
              << std::flush;
    ExpandableBuilder::staticAddOrSub
      (stats, histo, stats_change, histo_change, add);
    std::int64_t neg2 = negativeBins(*histo);
    std::cout << (new_code ? "OUT" : "out")
              << " stats " << stats->toString()
              << " histo " << histo->toString(false)
              << " negative " << neg2 - neg1 << " (now " << neg2 << ")"
              << "\n"
              << std::flush;
    //if (neg2 != neg1)
    //  BREAK();
  }
#endif

  /**
   * \brief Add or subtract this buffer's samples in the statistics.
   *
   * \details
   *
   * This is used to keep track of statistics while data is being updated.
   * Used for both the old and the new api.
   *
   * - New API, in ZgyInternalBulk::writeAlignedBrickList(),
   *   collect changes in the brick list before and after the r/m/w
   *   processing. The part of r/m/w that reads in old data and
   *   updates statistics is always active. In contrast to the old code.
   *   The new code Works on a temporary because the change should
   *   not be applied if there is an exception.
   *   Locking is needed because multiple bricks might be processed
   *   in parallel.
   *
   * - Old API, in ZgyInternalBulk::writeRegion(), collect changes
   *   only if the file is opened for update. For practical reasons
   *   changes will also be tracked for newly created files when
   *   actual r/m/w is needed. In that case the statistics get
   *   discarded later. Locking is not needed because the entire
   *   buffer provided by the application will be processed at once.
   *   Exceptions (such as trying r/m/w on a compressed file) might
   *   cause the histogram to get out of sync. This is a bug. But,
   *   unlike the new code where the bug has been fixed, the old code
   *   will only see the problem on update. Which is rately used.
   *   There is no unit test to catch this bug, so I am not 100% sure.
   *
   * The code tries to keep the statistical range somewhat correct by
   * trimming it to the first and last histogram bin with samples.
   * See HistogramBuildet::operator+=() and StatisticData::trimRange().
   * TODO-WIP-BrickedAPI: the statistics min/max might end up too
   * narrow if samples in the histogram are clipped. Which might
   * happen when a file is open for update, where histogram range
   * cannot grow. statistics will also be wrong. For consistency,
   * consider to not even try. Just document that the range will never
   * shrink as long as the derived information is generated
   * incrementally.
   *
   * It would work for the int8 case, or any other type where the number
   * of bins equals the number of discrete values in the valuetype.
   *
   * The method could have been moved inside class DataBuffer which
   * would remove the inelegant template kludge. It would also avoid
   * duplicating the functionality here and in GenLodImpl::
   * _accumulateT(). But class DataBuffer is bloated enough. I don't
   * want it aware of HistogramBuilder etc.
   *
   * Note that tracking needs to be done in storage values in case
   * sample values are clipped while converting float to integral.
   *
   * Subtle issue: Prefer (stats - old + new) to (stats + (new - old))
   * because of an obscure issue with count == 0 where the code might
   * incorrectly assume that the ssq etc. should not be updated.
   *
   * Thread safety:
   * The function knows that the input stats and histo pointers might
   * refer to shared state. It will place a lock, using the optionally
   * supplied mutex, while accessing those. Arguably this is not very
   * elegant. It might be better to return the histogram builder to
   * the caller and have the caller lock.
   */
  template <typename T>
  static void
  trackChangesT(const std::shared_ptr<const DataBuffer>& data_in,
                const std::pair<index3_t,index3_t>& range,
                StatisticData *stats, HistogramData *histo,
                bool add, std::mutex& mutex,
                BulkAdHocTimers& timers)
  {
    SimpleTimerEx t10(timers.mt10);
    const std::int64_t len = range.second[0]*range.second[1]*range.second[2];
    if (!data_in || !stats || !histo || len == 0)
      return;                   // Not tracking changes.
    std::shared_ptr<const DataBufferNd<T,3>> data =
      std::dynamic_pointer_cast<const DataBufferNd<T,3>>(data_in);
    if (!data)
      return;                   // Maybe throw in this case?
    // Make an empty builder with identical histogram range and isfixed.
    std::unique_lock<std::mutex> lck(mutex);
    ExpandableBuilder hb(*histo, StatisticData(), /*copy=*/false);
    lck.unlock();
    t10.stop();
    if (data->isScalar()) {
      SimpleTimerEx t11(timers.mt11);
      hb.add(data->data(), data->data() + 1);
      hb *= len;
    }
    else if (range.first == index3_t{0,0,0} && range.second == data->size3d()) {
      SimpleTimerEx t12(timers.mt12);
      hb.add(data->data(), data->data() + len); // HIGH CPU CONSUMPTION
    }
    else {
      SimpleTimerEx t13(timers.mt13);
      // Assumes data->stride3d() == index3_t{size[1]*size[2], size[2], 1}
      // TODO: Too much code assumes default stride. Change DataBuffer to
      // remove its half-hearted attempt of supporting arbitrary strides?
      // TODO-WIP-BrickedAPI: Can use the InternalZGY::EdgeBrick iterator
      // to *maybe* speed up processing, reduce locks, etc. Or make a
      // templated or overloaded add(data, size, roi) which might be
      // more efficient. Needs proper benchmarking.
      const std::array<std::int64_t,3> size = data->size3d();
      for (std::int64_t ii=range.first[0]; ii<range.first[0]+range.second[0]; ++ii) {
        for (std::int64_t jj=range.first[1]; jj<range.first[1]+range.second[1]; ++jj) {
          const T* ptr = data->data() + ii*size[1]*size[2] + jj*size[2];
          hb.add(ptr + range.first[2], ptr + range.first[2] + range.second[2]);
        }
      }
    }
    SimpleTimerEx t14(timers.mt14);
    std::lock_guard<std::mutex> lk(mutex);
    t14.stop();
    SimpleTimerEx t15(timers.mt15);
    ExpandableBuilder::staticAddOrSub
      (stats, histo, hb.getstats(), hb.gethisto(), add);
    // the addChanges wrapper is not useful here, because stats and histo
    // will be temporaries. Also, no locking.
    // addChanges(stats, histo, hb.getstats(), hb.gethisto(), add, data->size3d(), range.second, true);
  }

  /**
   * \copydoc TrackChangesT
   * This method has the usual boilerplate for not-quite templated code.
   */
  static void
  trackChanges(const std::shared_ptr<const DataBuffer>& data_in,
               const std::pair<index3_t,index3_t>& range,
               StatisticData *stats, HistogramData *histo,
               bool add, std::mutex& mutex, BulkAdHocTimers& timers)
  {
    switch (data_in->datatype()) {
    case RawDataType::SignedInt8:
      trackChangesT<std::int8_t>(data_in, range, stats, histo, add, mutex, timers);
      break;
    case RawDataType::SignedInt16:
      trackChangesT<std::int16_t>(data_in, range, stats, histo, add, mutex, timers);
      break;
    case RawDataType::Float32:
      trackChangesT<float>(data_in, range, stats, histo, add, mutex, timers);
      break;
    //case RawDataType::UnsignedInt8:
    //case RawDataType::UnsignedInt16:
    //case RawDataType::UnsignedInt32:
    //case RawDataType::SignedInt32:
    //case RawDataType::IbmFloat32:
    default:
      throw OpenZGY::Errors::ZgyInternalError("Unrecognized valuetype.");
    }
  }

  /**
   * Return true if cube a is equal to or larger than cube b in all dimensions.
   */
  static bool
  covers(const std::array<int64_t,3>& astart,
         const std::array<int64_t,3>& asize,
         const std::array<int64_t,3>& bstart,
         const std::array<int64_t,3>& bsize)
  {
    bool ok = true;
    for (int ii=0; ii<3; ++ii)
      ok = ok && (astart[ii] <= bstart[ii] &&
                  astart[ii] + asize[ii] >= bstart[ii] + bsize[ii]);
    return ok;
  }

#if 0 // Yagni?
  /**
   * Return true if cubes a and b have at least one sample in common.
   */
  static bool
  overlaps(const std::array<int64_t,3>& astart,
           const std::array<int64_t,3>& asize,
           const std::array<int64_t,3>& bstart,
           const std::array<int64_t,3>& bsize)
  {
    bool ok = true;
    for (int ii=0; ii<3; ++ii)
      ok = ok && (astart[ii] < bstart[ii] + bsize[ii] &&
                  bstart[ii] < astart[ii] + asize[ii]);
    return ok;
  }

  /**
   * Return true if cubes a and b have at least one sample in common.
   * Return the actual intersection if it exists. Return unspecified
   * start and zero size if the function return is false.
   */
  static bool
  overlaps(
       const std::array<int64_t,3>& astart,
       const std::array<int64_t,3>& asize,
       const std::array<int64_t,3>& bstart,
       const std::array<int64_t,3>& bsize,
       std::array<int64_t,3>& rstart,
       std::array<int64_t,3>& rsize)
  {
    bool ok = true;
    for (int ii=0; ii<3; ++ii) {
      rstart[ii] = std::max(astart[ii], bstart[ii]);
      rsize[ii] = (std::min(astart[ii] + asize[ii],
                            bstart[ii] + bsize[ii]) -
                   rstart[ii]);
      ok = ok && rsize[ii] > 0;
    }
    if (!ok) {
      rstart = std::array<int64_t,3>{0,0,0};
      rsize = std::array<int64_t,3>{0,0,0};
    }
    return ok;
  }
#endif
}

struct IJK
{
  std::int64_t i0;
  std::int64_t j0;
  std::int64_t k0;
  std::int64_t ni;
  std::int64_t nj;
  std::int64_t nk;
  std::ostream& dump(std::ostream &os) const {
    os << "start (" << i0 << ", " << j0 << ", " << k0 << ")"
       << " size (" << ni << ", " << nj << ", " << nk << ")";
    return os;
  }
  friend std::ostream& operator<<(std::ostream& os, const IJK& ijk) {
    return ijk.dump(os);
  }
};

/**
 * Argument package used by prepareWriteOneBrick(), prepareWriteOneNormalBrick(),
 * writeOneConstantBrick(). One size fits all. Which is in general
 * not a god idea. But there is a compromise between readablilty
 * and ease of implementing multi-threading.
 *
 * The original implementation processed a single brick at a time.
 * Using this class it becomes reasonably easy to pass around
 * lists of bricks to be handled.
 *
 *   * brickpos, lod: Always set and not modified by layers below.
 *   * data:          Actual bulk data to output. May be changed.
 *   * compressor:    Might be cleared. N/A if brick ends up constant.
 *   * fileoffset:    Only when known. Never known yet in prepareWriteOneBrick(),
 *                    and N/A for constant bricks.
 */
struct WriteBrickArgPack
{
  std::array<std::int64_t,3>  brickpos;
  std::int32_t                lod;
  std::shared_ptr<DataBuffer> data;
  compressor_t                compressor; // TODO-Low, caller use std::ref?
  //bool                        use_compressor; // or this?
  std::int64_t                fileoffset;
  WriteBrickArgPack(const std::array<std::int64_t,3>& brickpos_in,
                    std::int32_t lod_in,
                    const std::shared_ptr<DataBuffer>& data_in,
                    const compressor_t& compressor_in,
                    std::int64_t fileoffset_in)
    : brickpos(brickpos_in)
    , lod(lod_in)
    , data(data_in)
    , compressor(compressor_in)
    , fileoffset(fileoffset_in)
  {
  }
  std::string toString() const
  {
    std::stringstream ss;
    ss << "pos=" << brickpos
       << ", lod=" << lod
       << ", size=" << data->size3d();
    if (data->isScalar())
      ss << ", scalar=" << data->scalarAsDouble();
    if (compressor)
      ss << ", compressor=*";
    if (fileoffset)
      ss << ", fileoffset=" << std::hex << fileoffset << std::dec;
    return ss.str();
  }
};

/**
 * Argument package used by writeOneNormalBrickWithRetry() which is at a lower level
 * than prepareWriteOneBrick() etc. and it got too awkward to use the
 * one size fits all rule.
 *
 * There are new fields rawdata (replaces data) and brickstatus.
 * Also the compressor is no longer needed as compression is now done.
 * But we need to know whether compression was done because that might
 * affect alignment whan a new brick is allocated.
 */
struct WriteNowArgPack
{
  std::array<std::int64_t,3>  brickpos;
  std::int32_t                lod;
  std::int64_t                fileoffset;
  rawdata_t                   rawdata;
  BrickStatus                 brickstatus;
  bool                        align;
  WriteNowArgPack(const std::array<std::int64_t,3>& brickpos_in,
                  std::int32_t lod_in,
                  std::int64_t fileoffset_in,
                  const rawdata_t rawdata_in,
                  BrickStatus brickstatus_in,
                  bool align_in)
    : brickpos(brickpos_in)
    , lod(lod_in)
    , fileoffset(fileoffset_in)
    , rawdata(rawdata_in)
    , brickstatus(brickstatus_in)
    , align(align_in)
  {
  }
  std::string toString() const
  {
    std::stringstream ss;
    ss << "pos=" << brickpos
       << ", lod=" << lod
       << ", size=" << rawdata.second;
    if (fileoffset)
      ss << ", fileoffset=" << std::hex << fileoffset << std::dec;
    // Should use symbolic names for enums, but this is just for verbose logs.
    ss << ", brickstatus=" << (int)brickstatus;
    if (align)
      ss << ", alignment needed";
    return ss.str();
  }
};

/**
 * TODO-Low might want to fold this into LookupTable::LutInfo.
 *
 * Add position in samples. LookupTable::getBrickFilePosition()
 * cannot easily store this because it only has the brick position
 * and would need to know the brick size to get sample position.
 *
 * Add constvalue after decoding.
 */
struct LutInfoEx : public LookupTable::LutInfo
{
  // Inherited: status, offset_in_file, size_in_file, raw_constant;
  IJK    survey_position;
  double double_constvalue;

  LutInfoEx(const LookupTable::LutInfo& info, const IJK& pos_in, double constvalue_in)
    : LookupTable::LutInfo(info)
    , survey_position(pos_in)
    , double_constvalue(constvalue_in)
  {
  }
};

/**
 * Duplicated between impl/bulk.cpp and impl/meta.cpp but sets
 * different flags.
 *
 * Start a critical section where any exception means that the
 * owner class should be permanently flagged with _is_bad = True.
 * Typically this is used to prevent secondary errors after a
 * write failure that has most likely corrupted the entire file.
 * The exception itself will not be caught.
 *
 * The _is_bad flag normally means that any further attempts
 * to access this class, at least for writing, will raise a
 * ZgyCorruptedFile exception. Regardless of what the exception
 * was that caused the flag to be set.
 *
 * C++ note: Unlike Python it isn't trivial (or portable) to check
 * whether a destructor is being called due to leaving scope normally
 * or due to unwinding an exception. So in the C++ version the code
 * should explicitly call disarm() at the end of the critical section.
 *
 * Thread safety:
 * The class itself is not thread safe but this should not be an issue,
 * because instances are meant to be short lived. Typically used inside
 * a method and not acceesible outside.
 */
class ZgyInternalBulk::ErrorsWillCorruptFile
{
  ZgyInternalBulk *_owner;
public:
  explicit ErrorsWillCorruptFile(ZgyInternalBulk* owner) : _owner(owner)
  {
  }
  ~ErrorsWillCorruptFile()
  {
    if (_owner)
      _owner->set_errorflag(true);
    _owner = nullptr;
  }
  void disarm()
  {
    _owner = nullptr;
  }
};

/**
 * \details \callgraph \callergraph
 *
 * ZgyInternalBulk is used both from ZgyReader, which is in general not
 * allowed to change any state to help make it thread safe, and from the
 * non threadsafe ZgyWriter.
 *
 * To mitigate the risk of accidentally modifying data using a ZgyReader
 * there is both a const and a mutable pointer to the underlying
 * ZgyInternalMeta. When instanciated from a ZgyReader the mutable pointer
 * will be empty. When instanciated from a ZgyWriter the two pointers
 * will be identical. This mitigation by itself will cause a null pointer
 * exception if trying to modify data that shouldn't be change.
 * This of course assumes that ZgyInternalMeta is const-correct so that
 * no state can be changed by using a const pointer to a ZgyInternalMeta
 * instance.
 *
 * An additional mitigation is to not use the _metadata_rw directly but
 * instead call a get_metadata_rw() method. If ZgyInternalBulk is
 * const-correct then you get a compile time error if trying to call
 * get_metadata_rw() from inside a method declared as const. If not,
 * get_metadata_rw() is still preferable because it can raise a proper
 * error message instead of the null pointer exception.
 */
ZgyInternalBulk::ZgyInternalBulk(
     const std::shared_ptr<IFileADT>& file,
     const std::shared_ptr<const ZgyInternalMeta>& metadata,
     const std::shared_ptr<ZgyInternalMeta>& metadata_rw,
     bool compressed_write,
     const std::array<float, 2>& force_histo_range,
     const LoggerFn& logger)
  : _file(file)
  , _metadata(metadata)
  , _metadata_rw(metadata_rw)
  , _update_mode(compressed_write? UpdateMode::Constant : UpdateMode::Always)
  , _compressed_write(compressed_write)
  , _is_bad(false)
  , _loggerfn(logger ? logger : LoggerBase::standardCallback(LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"), "openzgy-bulk: ", ""))
  , _new_mutex()
  , _new_modified_bricks()
  , _new_modified_stats()
  , _new_modified_histo()
  , _new_stathist_good(true)
  , _new_genlod()
  , _new_internal_lod_mode()
  , _force_histo_range(force_histo_range)
  , _timers()
{
  this->newTrackedBricksTryEnable();
}

/**
 * \brief Get hint about all constant region.
 *
 * \details \callgraph \callergraph
 *
 * Check to see if the specified region is known to have all samples
 * set to the same value. Returns a pair of (is_const, const_value).
 *
 * The function only makes inexpensive checks so it might return
 * is_const=false even if the region was in fact constant. It will not
 * make the opposite mistake. This method is only intended as a hint
 * to improve performance.
 *
 * For int8 and int16 files the caller may specify whether to scale
 * the values or not. Even if unscaled the function returns the value
 * as a double.
 */
std::pair<bool,double>
ZgyInternalBulk::readConstantValue(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod, bool as_float) const
{
  validateUserLod(lod);
  validateUserPosition(start, size, lod);
  const double missing_val = missingValue(/*as_float=*/false);
  std::vector<LutInfoEx> bricklist = partsNeeded(start, size, lod);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  double result = nan;
  bool first = true;
  for (const LutInfoEx& brick : bricklist) {
    if (blogger(5))
      blogger(5, std::stringstream()
              << "brick " << brick.survey_position
              << " " << (int)brick.status << "\n");
    switch (brick.status) {
    case BrickStatus::Constant:
      // Two NaN values should compare equal in this test.
      if (!first && result != brick.double_constvalue &&
          !(std::isnan(result) && std::isnan(brick.double_constvalue)))
        return std::make_pair(false,nan);
      result = brick.double_constvalue;
      break;
    case BrickStatus::Missing:
      if (!first && result != missing_val)
        return std::make_pair(false,nan);
      result = missing_val;
      break;
    default:
      return std::make_pair(false,nan);
    }
    first = false;
  }
  if (as_float && std::isfinite(result)) {
    result = (result * this->_metadata->ih().storagetofloat_slope() +
              this->_metadata->ih().storagetofloat_intercept());
  }

  // Only logging success, because on failure the application will
  // need to read the actual data. Which will trigger the logging
  // in readToExistingBuffer.
  if (blogger(2))
    blogger(2, std::stringstream()
            << "read(start="
            << "(" << start[0]
            << "," << start[1]
            << "," << start[2]
            << ")"
            << ", size="
            << "(" << size[0]
            << "," << size[1]
            << "," << size[2]
            << ")"
            << ", lod=" << lod
            << std::boolalpha << ", as_float=" << as_float << ")"
            << " => constant " << result);

  return std::make_pair(true,result);
}

/**
 * Create a functor that will accept data read from the lower levels
 * and pass it on to deliverOneBrick() for eventual delivery to the user.
 */
ReadRequest::delivery_t
ZgyInternalBulk::createDeliverance(
     const std::shared_ptr<DataBuffer>& result,
     const std::array<std::int64_t,3>& start,
     bool as_float,
     const LutInfoEx& brick) const
{
  // The deliverance functor is called when a single chunk of data
  // has been read by the low level xx_readv(). Only a void* pointer
  // and a size is passed in. The lambda captures the output buffer
  // and the positions of the user's read request and where in the
  // output buffer the data being delivered now should be placed.
  // The functor also knows whether the data, when it arrives, needs
  // to be decompressed. If used from the brick-based api, the data
  // will always be stored at the beginning of one of the the
  // provided brick-sized buffer.

  // result, start, and as_float come straight from user code.
  // b.position and b.status vary per brick which is why I copy them
  // to make sure the up to date contents are captured. raw and
  // rawsize (the arguments to the lambda) is what was read from
  // storage.

  // TODO-WIP-BrickAPI: In brick mode, more shortcuts might in
  // theory be possible by letting xx_readv() write directly into
  // the user's buffer(s). Not using deliverance(). The messiness
  // might not be worth it.

  // In brick mode, there is no reshaping bricks read from physical
  // file into the user's arbitrary buffer. "start" differs for each
  // deliverance, and will be the same as "bpos". Which will be
  // start_in[ii] since we don't do any splitting or combining of
  // bricks. There should (TODO-WIP-BrickAPI check) be code that
  // does a shortcut and just does a memcpy into the target.
  // "result" also differs for each deliverance, and will be the
  // buffer for this part of the request.
  const std::array<std::int64_t,3> bpos
    {brick.survey_position.i0,
     brick.survey_position.j0,
     brick.survey_position.k0};
  const BrickStatus bstatus = brick.status;
  return [this,result,start,as_float,bpos,bstatus]
    (ReadRequest::data_t raw, std::int64_t rawsize) {
           this->deliverOneBrick
             (result, start, bpos, raw, rawsize, bstatus, as_float);
         };
}

namespace {
/**
 * Debug code. Log what we are going to do with this single brick.
 * Used only in readToExistingBuffer(); refactored out to avoid clutter.
 */
static void
logReadingSingleBrick(
     const LutInfoEx& brick,
     const ZgyInternalBulk::LoggerFn& logger)
{
  if (!logger(2, ""))
    return;
  std::stringstream ss;
  ss << "  Read ";
  switch (brick.status) {
  case BrickStatus::Missing:
    ss << "missing brick at " << brick.survey_position
       << ", use \"defaultstorage\"\n";
    break;
  case BrickStatus::Constant:
    ss << "constant brick at " << brick.survey_position << "\n";
    break;
  case BrickStatus::Normal:
    ss << "normal brick at " << brick.survey_position
       << " from file offset " << std::hex
       << (std::intptr_t)brick.offset_in_file << std::dec << "\n";
    break;
  case BrickStatus::Compressed:
    ss << "compressed brick at " << brick.survey_position
       << " from file offset " << std::hex
       << (std::intptr_t)brick.offset_in_file << std::dec
       << " size " << brick.size_in_file << "\n";
    break;
  default:
    ss << "brick of invalid type " << (int)brick.status << "\n";
    break;
  }
  logger(2, ss.str());
}

/**
 * Process, but possibly defer, a read request for a single brick.
 * The "deliverance" functor will be invoked directly for simple
 * cases such as constant bricks. Otherwise it gets wrapped into
 * a ReadRequest and returned to the caller. To be handled in
 * actualReadingAllBricks()
 */
static void
startReadingSingleBrick(
     const LutInfoEx& brick,
     const ReadRequest::delivery_t& deliverance,
     double missing_val,
     std::int64_t bytesperbrick,
     std::vector<ReadRequest> *requests)
{
  switch (brick.status) {
  case BrickStatus::Missing:
    deliverance
      (std::make_shared<double>(missing_val), sizeof(double));
    break;
  case BrickStatus::Constant:
    deliverance
      (std::make_shared<double>(brick.double_constvalue), sizeof(double));
    break;
  case BrickStatus::Normal:
    requests->push_back
      (ReadRequest{brick.offset_in_file, bytesperbrick, deliverance});
    break;
  case BrickStatus::Compressed:
    // TODO-Worry obscure corner case, might need to re-try if we didn't
    // get enough data. I try to prevent this with the rule that
    // no compressed brick may be larger than the uncompressed version.
    requests->push_back
      (ReadRequest{brick.offset_in_file, brick.size_in_file, deliverance});
    break;
  default:
    throw OpenZGY::Errors::ZgyInternalError("Internal error, bad brick status");
  }
}

/**
 * Process real brick read requests that were deferred and queued up
 * by startReadingSingleBrick. The main point of this convoluted
 * control flow is to be able to call the low level scatter/gather
 * read function with as many bricks as possible. When reading from
 * the cloud, neighboring brick reads might also be consolidated.
 */
static void
actualReadingAllBricks(
     const std::vector<ReadRequest>& requests,
     const std::shared_ptr<IFileADT>& file,
     const ZgyInternalBulk::LoggerFn& logger)
{
  if (logger(2, ""))
    logger(2, std::to_string(requests.size()) + " read requests are queued\n");

  if (!requests.empty())
    file->xx_readv(requests, true, false, true, UsageHint::Data);

  // Note-Performance: If passing true in the second argument above, this
  // could help performance a lot. Especially for reading compressed files
  // where the user sends large requests without multithreading. Also
  // when finalizing compressed files. parallel_ok=true will cause the
  // decompression step to be multi-threaded. Also the conversion to
  // float (if needed) and the copy-out to the applicaton's buffer will
  // be multi-threaded. But there are caveats:
  //
  // * deliverOneBrick() must be thread safe.
  //
  // * The cloud backend doesn't honor the parallel_ok argument.
  //   While this would be a very good idea it is also rather difficult
  //   to implement.
  //
  // * There is a risk of creating too many threads if the application
  //   is doing its own multi threading. Ideally the user should
  //   be able to configure this.
  //
  // * Ditto for the low resolution generator. It can probably speed
  //   up by reading with 4 (but only 4) threads. So this isn't as
  //   efficient as setting parallel_ok=true here with respect to
  //   speeding up compression. But it might help I/O against the
  //   cloud. Which the setting here won't.
  //
  // * See commants in LocalFileLinux::xx_readv() and deliverOneBrick().
  //   And GenLodImpl::_calculate().
}
} // anonymous namespace

/**
 * \brief Read arbitrary region to a buffer owned by the caller.
 *
 * \details \callgraph \callergraph
 *
 * Read bulk data starting at "start" in index space and store the
 * result in the provided DataBuffer. Start should be in the range
 * (0,0,0) to Size-1. The count of samples to read is implied by the
 * size of the DataBuffer that is passed in. The valid data types for
 * the result are float32 (in which case samples stored as int8 or
 * int16 will be scaled) or the files's storage value type (in which
 * case there is no scaling). It is valid to pass a size that
 * includes the padding area between the survey and the end of the
 * current brick. But not past that point.
 */
void
ZgyInternalBulk::readToExistingBuffer(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      int32_t lod, bool as_float) const
{
  validateUserLod(lod);
  validateUserPosition(start, result->size3d(), lod);
  RawDataType expected_dtype = as_float ? RawDataType::Float32 : _metadata->ih().datatype();
  // When called from ZgyReader.read(), the check on expected_dtype represents
  // a user error trying to convert to an intergral type while the others
  // should not be possible to trigger.
  if (!result || result->isScalar() || !result->voidData())
    throw OpenZGY::Errors::ZgyInternalError("Reading to missing or wrong buffer type.");
  if (result->datatype() != expected_dtype)
    throw OpenZGY::Errors::ZgyUserError("Requested data type not supported for this file.");

  if (blogger(2))
    blogger(2, std::stringstream()
            << "read(start="
            << "(" << start[0]
            << "," << start[1]
            << "," << start[2]
            << "), size="
            << "(" << result->size3d()[0]
            << "," << result->size3d()[1]
            << "," << result->size3d()[2]
            << "), lod=" << lod
            << std::boolalpha << ", as_float=" << as_float << ")");

  // Make a separate pass to gather all the bricks we need to read.
  // The lower levels might fetch some of them in parallel and we might be
  // able to combine bricks to read larger blocks at a time. Both changes
  // can have a dramatic performance impact effect on cloud access.
  // In brick mode we obviosly have that list of bricks already.

  // In generic read mode, call partsNeeded() to split the region into
  // bricks and then get some information about those bricks. In brick
  // mode, only the last step, which is makeLutInfoEx(), is needed.
  std::vector<LutInfoEx> bricklist = partsNeeded(start, result->size3d(), lod);

  // In generic read mode:
  // After all bricks have been processed, the padding past the
  // end if the survey might still not have been touched. Just in
  // case the request did in fact include such samples we will
  // initialize the entire result buffer to the default value.
  // TODO-Performance can I refactor the code so this isn't needed?
  // Or can I make it fill only the actual padding area?
  // Note, the long term goal is to make the generic read call
  // the brick api internally. This moves the issue. See below.
  // Since the outside samples in the user buffer are not touched
  // at all, the value that is filled needs to be converted.

  // In brick mode, there is no need to clear the buffers first.
  // Trust that the read will set every sample of the buffers.
  result->fill(missingValue(as_float));

  std::vector<ReadRequest> requests;
  for (const LutInfoEx& brick : bricklist) {
    auto deliverance = createDeliverance(result, start, as_float, brick);
    logReadingSingleBrick(brick, this->_loggerfn);
    startReadingSingleBrick(brick, deliverance, missingValue(false),
                            this->_metadata->ih().bytesperbrick(),
                            &requests);
  }
  actualReadingAllBricks(requests, this->_file, this->_loggerfn);
}

/**
 * \brief Read bricks to a list of buffers owned by the caller.
 *
 * \details \callgraph \callergraph
 *
 * Read one or more bricks of bulk data given the starting point in
 * index space and the existing brick-sized buffers to store the
 * result in. The function is very similar to readToExistingBuffer(),
 * but works on a list of bricks instead of an arbitrary region.
 *
 * Attempts to read to a missing or scalar buffer are quietly skipped.
 *
 * Attempts to read outside the survey has already been handled by
 * readBricksInternal(). Ditto for expedited single brick reads.
 * That function sets the already handled buffers to null or scalar.
 */
void
ZgyInternalBulk::readBricksToExistingBuffer(
     const std::vector<std::shared_ptr<DataBuffer>>& result,
     const std::vector<std::array<std::int64_t,3>>& start,
     int32_t lod, bool as_float) const
{
  validateUserLod(lod);
  checkBricksInternal(start, result, as_float, /*strict=*/false);
  std::vector<ReadRequest> requests;
  for (std::size_t ii = 0; ii < result.size() && ii < start.size(); ++ii) {
    if (!result[ii] || result[ii]->isScalar() || !result[ii]->voidData())
      continue;
    const LutInfoEx brick = makeLutInfoEx(start[ii], lod);
    auto deliverance = createDeliverance(result[ii], start[ii], as_float, brick);
    logReadingSingleBrick(brick, this->_loggerfn);
    startReadingSingleBrick(brick, deliverance, missingValue(false),
                            this->_metadata->ih().bytesperbrick(),
                            &requests);
  }
  actualReadingAllBricks(requests, this->_file, this->_loggerfn);
}

/*
 * TODO-WIP-BrickedAPI: The old ZgyInternalBulk::readToNewBuffer() did:
 * Pass check_constant=true to check extra hard for all-constant data.
 * A region written with all samples identical, as opposed to a region
 * flagged as constant without taking up space, will also be detected.
 */

/**
 * \brief Return true if the write is updating all samples in the survey.
 *
 * \details
 * Whether the caller also tries to update samples outside the survey
 * is irrelevant.
 *
 * When called from (not yet implemented) writeConstBricksInternal,
 * it would be nice if that function could give us a hint that the
 * write covers the whole survey. Instead of figuring it out locally.
 *
 * Not done: The function might check the actual brick positions
 * instead of assuming that the same brick is not included more than
 * once. Which would be invalid. But defensive coding is nice.
 */
bool /*static*/
ZgyInternalBulk::coversEntireSurvey(
     const std::vector<std::array<int64_t,3>>& position,
     const std::array<int64_t,3>& surveysize,
     const std::array<int64_t,3>& bricksize)
{
  // Survey size in bricks.
  const std::array<std::int64_t,3> ssize
    {(surveysize[0] + bricksize[0] - 1) / bricksize[0],
     (surveysize[1] + bricksize[1] - 1) / bricksize[1],
     (surveysize[2] + bricksize[2] - 1) / bricksize[2]};
  const std::int64_t survey_bricks = ssize[0] * ssize[1] * ssize[2];

  if ((std::int64_t)position.size() < survey_bricks)
    return false; // Shortcut

  // Count number of bricks inside the survey, including those with padding.
  // Already know that the start position is brick-aligned, but check this.
  std::int64_t provided_bricks{0};
  for (const auto& pos : position) {
    bool inside{true};
    for (int dim = 0; dim < 3; ++dim) {
      if ((pos[dim] % bricksize[dim]) != 0)
        throw OpenZGY::Errors::ZgyInternalError
          ("Misaligned brickpos in coversEntireSurvey()");
      if (pos[dim] < 0 || pos[dim] >= surveysize[dim])
        inside = false;
    }
    if (inside)
      ++provided_bricks;
  }

  // Note that if there are both duplicated bricks and missing bricks
  // then we might not detect this, and might also return a false positive.
  if (provided_bricks > survey_bricks)
    throw OpenZGY::Errors::ZgyInternalError
      ("Duplicated brickpos in coversEntireSurvey()");

  return provided_bricks == survey_bricks;
}

/**
 * \brief Return true if the write is setting all samples in the survey to
 * the same sample value.
 *
 * \details
 * Missing buffers, i.e. finding an empty shared_ptr, will return false.
 *
 * Not done: Samples in the padding area and bricks completely outside
 * the survey should be ignored. In practice, this thorough check
 * might not be needed.
 */
std::pair<bool, double> /*static*/
ZgyInternalBulk::isSameScalarValue(
     const std::vector<std::shared_ptr<DataBuffer>>& data)
{
  if (data.empty() || !data.front() || !data.front()->isScalar())
    return std::pair<bool, double>(false, 0);
  const double value = data.front()->scalarAsDouble();
  if (std::isnan(value)) {
    for (const auto& buffer : data)
      if (!buffer || !buffer->isScalar() || !std::isnan(buffer->scalarAsDouble()))
        return std::pair<bool, double>(false, 0);
  }
  else {
    for (const auto& buffer : data)
      if (!buffer || !buffer->isScalar() || buffer->scalarAsDouble() != value)
        return std::pair<bool, double>(false, 0);
  }
  return std::pair<bool, double>(true, value);
}

/**
 * \brief Return true if the write is setting all samples in the survey to
 * the same sample value. Which will be: data.front()->scalarAsDouble().
 *
 * \details
 * Not done: Samples in the padding area and samples completely outside
 * the survey should be ignored. In practice, this thorough check
 * might not be needed. Just see if all the buffers are scalar and
 * have the same value.
 *
 * When true, the new version of the statistics and histogram reader
 * can be completely reset after this write has completed. This allows
 * the histogram range and the statistical min/max to shrink. It also
 * allows to switch back to an expandable histogram if applicable.
 * And it fixes any bad histogram e.g. if it was read from an existing
 * ZGY file which didn't have a good histogram. If none of these
 * features are needed, the special handling could have been skipped.
 *
 * When true, the now redundant tracking before and after the
 * write can, but doesn't need to, be skipped. This would slightly
 * improve performance. Probably insignificant, so don't do this if it
 * causes additional spaghettification of the code.
 *
 * Not done: If all samples are being updated but not all to the same
 * value, it is still possible to get the benefits listed above. That
 * requires more code. It is also a very unlikely scenario, because
 * even a reasonably sized survey would be too large to write in a
 * single call. If we want it anyway, consider detecting this case
 * and inserting a dummy writeconst of the entire survey
 */
std::pair<bool, double> /*static*/
ZgyInternalBulk::isResettingEntireSurvey(
     const std::vector<std::array<std::int64_t,3>>& position,
     const std::array<std::int64_t,3>& surveysize,
     const std::array<std::int64_t,3>& bricksize,
     const std::vector<std::shared_ptr<DataBuffer>>& data)
{
  // Educated guess: This test is probably the cheaper one.
  std::pair<bool, double> isconst = isSameScalarValue(data);
  if (!isconst.first)
    return isconst;
  if (!coversEntireSurvey(position, surveysize, bricksize))
    return std::pair<bool, double>(false, 0);
  return isconst;
}

/**
 * \brief Discard the old statistics and histogram and initialize them
 * to reflect an all-constant survey.
 *
 * \details
 * This is where the in-memory histogram size is chosen. And whether
 * the histogram should be expandable or not.
 *
 * This function needs to be called explicitly after creating a new
 * file. The existing histogram will in that case be incorrect. It is
 * beneficial to also call it in other all-const cases. Such as
 * opening an existing file that has no data yet. Or after an explicit
 * writeconst() of the entire survey called from the application code.
 * See isResettingEntireSurvey() for more details.
 *
 * Might be called indirectly from the constructor. Must not be virtual.
 */
void
ZgyInternalBulk::newTrackedBricksSetAllConst(double value_in)
{
  const IInfoHeaderAccess& ih = _metadata_rw->ih();
  switch (ih.datatype()) {
  default:
  case RawDataType::Float32:
    {
      if (_force_histo_range[0] < _force_histo_range[1]) {
        // This user provided histogram range is honored on initial create,
        // and also if there is a full reset of the survey to constant-value.
        HistogramBuilder hb(256, _force_histo_range[0], _force_histo_range[1]);
        float fvalue[1]{static_cast<float>(value_in)};
        hb.add(&fvalue[0], &fvalue[1]);
        hb *= (ih.size()[0] *ih.size()[1] *ih.size()[2]);
        _new_modified_stats = std::make_shared<StatisticData>(hb.getstats());
        _new_modified_histo = std::make_shared<HistogramData>(hb.gethisto());
      }
      else {
        ExpandableBuilder hb(4096); // Careful! Not polymorphic.
        float fvalue[1]{static_cast<float>(value_in)};
        hb.add(&fvalue[0], &fvalue[1]);
        hb *= (ih.size()[0] *ih.size()[1] *ih.size()[2]);
        _new_modified_stats = std::make_shared<StatisticData>(hb.getstats());
        _new_modified_histo = std::make_shared<HistogramData>(hb.gethisto());
      }
      break;
    }
  case RawDataType::SignedInt16:
    {
      // TODO-WIP-BrickedAPI: Consider a dynamic histogram here.
      // Or even a fixed, 64k internal histogram. In the latter
      // case, allow finalize to squeeze and resample also
      // previously fixed range histogram. Some caveats regarding
      // zero-centric. The user's coding range will most likely
      // be zero-centric, but we need anti-zero-centric.
      HistogramBuilder hb(256, -32768, +32767);
      std::int16_t ivalue[1]{static_cast<std::int16_t>(value_in)};
      hb.add(&ivalue[0], &ivalue[1]);
      hb *= (ih.size()[0] *ih.size()[1] *ih.size()[2]);
      _new_modified_stats = std::make_shared<StatisticData>(hb.getstats());
      _new_modified_histo = std::make_shared<HistogramData>(hb.gethisto());
      break;
    }
  case RawDataType::SignedInt8:
    {
      HistogramBuilder hb(256, -128, +127);
      std::int8_t ivalue[1]{static_cast<std::int8_t>(value_in)};
      hb.add(&ivalue[0], &ivalue[1]);
      hb *= (ih.size()[0] *ih.size()[1] *ih.size()[2]);
      _new_modified_stats = std::make_shared<StatisticData>(hb.getstats());
      _new_modified_histo = std::make_shared<HistogramData>(hb.gethisto());
      break;
    }
  }
  if (blogger(1))
    blogger(1, std::stringstream()
            << "Will Track new changes. Survey is now all const. "
            << _new_modified_stats->toString()
            << " histo " << _new_modified_histo->toString());
}

/**
 * \brief Read existing histogram and statistics from a file being
 * opened for update.
 *
 * \details
 *
 * Will be called from the constructor. Must not be virtual.
 *
 * Try to detect the case where the file to be updated is still empty
 * or all-const, and call newTrackedBricksSetAllConst instead.
 */
void
ZgyInternalBulk::newTrackedBricksInitFromFile()
{
  const IHistHeaderAccess& hh = _metadata_rw->hh();
  const IInfoHeaderAccess& ih = _metadata_rw->ih();
  // Trust that the statistics and histogram read from file are good.
  // CAVEAT: If a user has created and not finalized a file in an older
  // version of the library then this probably won't work.
  // The same applies to a completely empty file, as it will have
  // empty statistics instead of counting all the defaultvalue.
  // Actually, updating *any* file written by old ZgyPublic will
  // likely end up with bogus histogram and statistics.
  // When initialized from file, always create a fixed range histogram
  // with the size (currently always 256) that it has on file.
  _new_modified_stats = std::make_shared<StatisticData>
    (ih.scnt(), /*inf=*/0, ih.ssum(), ih.sssq(), ih.smin(), ih.smax());
  _new_modified_histo = std::make_shared<HistogramData>
    (hh.bins(), (int)hh.bincount(), hh.minvalue(), hh.maxvalue());
  if (blogger(2))
    blogger(2, std::stringstream()
            << "Will track new changes. Initial stats (real)    "
            << _new_modified_stats->toString()
            << " histo " << _new_modified_histo->toString());
  // Statistics saved in the ZGY file uses float data, but internal
  // calculation needs to be in storage because otherwise any clipped
  // values would mess things up.
  const std::array<double,2> factors = ih.storagetofloat();
  _new_modified_stats->scale(factors[1], factors[0] + factors[1], 0, 1);
  _new_modified_histo->scale(factors[1], factors[0] + factors[1], 0, 1);
  if (blogger(1))
    blogger(1, std::stringstream()
            << "Will track new changes. Initial stats (storage) "
            << _new_modified_stats->toString()
            << " histo " << _new_modified_histo->toString());
}

/**
 * In the new mechanism, incremental tracking of statistics and histogram
 * is always enabled when a file is open for write.
 *
 * Call this function immediately after creating a file or writing a
 * constant value to the entire survey. Also when opening a file for update.
 *
 * Technically there are still a few situations where tracking is
 * redundant. But my gut feeling says to leave it in.
 * Example: Compressed files need a full lowres build anyway (plan C).
 * The code might have skipped reading old data (if not r/m/w).
 * But compressed bricks cannot be updated. So reading old data is cheap.
 * Similarly, if an inconsistency is added that trips _new_stathist_good
 * (such as discovering that the statistical min/max are too narrow)
 * then we probably end up with a full rebuild. But maybe the caller
 * didn't want lowres and decided to trust the stats and histo anyway?
 *
 * TODO-WIP-BrickedAPI: Special case for all-const files that used to
 * contain real data. I.e. it contains inflated all-const bricks. See
 *
 *  - ZgyInternalBulk::newTrackedBricksTryEnable() // needs check_constant
 *  - ZgyInternalBulk::writeAlignedBrickList()     // refuses to leak
 *  - GenLodC::_readbricks()                       // check_constant=false
 *
 * In this case the statistics and histogram could have been reset, and
 * low resolution computation could have been speeded up. What happens
 * now is that the histogram ends up useless with min~=max.
 * See unit test reopen.alignment and possibly others.
 *
 * This method must be nonvirtual because it is also called from the
 * constructor.
 */
void
ZgyInternalBulk::newTrackedBricksTryEnable()
{
  _new_modified_bricks.reset();
  _new_modified_stats.reset();
  _new_modified_histo.reset();
  _new_stathist_good = true;
  if (!_metadata_rw) {
    blogger(1, "Will not track new changes. Not open for write.");
  }
  else {
    const IInfoHeaderAccess& ih = _metadata_rw->ih();
    const IHistHeaderAccess& hh = _metadata_rw->hh();
    std::pair <bool,double> global_value = readConstantValue
      (std::array<std::int64_t,3>{0,0,0},
       ih.size(), /*lod*/0, /*as_float*/false);

    // Keep track of where LOD bricks are needed.
    _new_modified_bricks = std::make_shared<TrackTouched>
      (ih.lodsizes(), ih.brickoffsets(), _loggerfn);

    if (global_value.first) {
      if (blogger(007), "")
        blogger(007, std::stringstream().flush()
                << "newTrackedBricksTryEnable:"
                << " all const " << global_value.second);
      newTrackedBricksSetAllConst(global_value.second);
    }
    else if (hh.samplecount() != 0 && ih.scnt() != 0 && ih.smin()<=ih.smax()) {
      // The file is being updated and has good histogram and statistics.
      // Note, if either is missing, both are assumed to be bad.
      // Note, statistical and/or histogram range might still be good.
      // For the purpose of tracking, this does not help.
      blogger(007, "newTrackedBricksTryEnable: read existing histogram");
      newTrackedBricksInitFromFile();
    }
    else {
      // Problem. SINGLEPASS always does incremental update of
      // statistics and histogram, so the values collected so far need
      // to be available. Maybe this file was written by an older
      // version of OpenZGY that allowed skipping this information.
      // Version 1 files are missing that information as well, but
      // they aren't updateable in the first place.
      //
      // Setting _new_stathist_good=false will:
      //
      // - Only track additions, not subtractions. Stats/histo will be
      //   wrong but still somewhat useful for the WeightedAverage
      //   algorithm.
      //
      // - Only use the stats/histo for decimation. Don't save it on close.
      //
      // - An alternative, not implemented, is to change any WeightedAverage
      //   to Average and do no tracking at all.
      //
      // To test this by hand: Patch this function to not see any good
      // stats/histo in the existing file being opened. Patch
      // check_incrlod() in test_api.cpp to set kludge=true. Patch
      // somewhere to save stats/histo in spite of them being bad. Run
      // Run api.incrlod3 and see that it succeeds. Can also run
      // api.incrlod2 but I haven't added expected results for that one.
      blogger(0, "Statistics and/or histogram are missing and won't be added.");
      newTrackedBricksSetAllConst(0);
      (*_new_modified_stats) *= 0;
      (*_new_modified_histo) *= 0;
      _new_stathist_good = false;
    }
  }
}

namespace {
template<typename T>
static void
fillMeT(void* data, const std::array<std::int64_t,3>& size, double value)
{
  std::fill(static_cast<T*>(data),
            static_cast<T*>(data) + size[0]*size[1]*size[2],
            static_cast<T>(value));
}

static void
fillMe(void* data, const std::array<std::int64_t,3>& size, double value, RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::Float32:    fillMeT<float>(data, size, value);        break;
  case RawDataType::SignedInt16:fillMeT<std::int16_t>(data, size, value); break;
  case RawDataType::SignedInt8: fillMeT<std::int8_t>(data, size, value);  break;
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized valuetype.");
  }
}
}

/**
 * \brief Short cut to read exactly one brick.
 *
 * \details \callgraph \callergraph
 * More efficient read of a region that happes to be exactly one brick,
 * not compressed, not converted int to float, and possibly other limitatons.
 * This can bypass a lot of logic and will elide at least one memory
 * allocation. Return false if the conditions are not met. The caller can
 * invoke this function unconditionally as long as readToExistingBuffer()
 * is used as a fallback.
 */
bool
ZgyInternalBulk::expeditedRead(const std::array<std::int64_t,3>& start, const std::array<std::int64_t,3>& size, void* data, int lod, RawDataType result_type) const
{
  // TODO-Performance: Used in ZgyWriter::read(), not just ZgyReader::read().
  // This ought to work but will require more testing.
  if (!expedited_read())
    return false;

  // Not posible due to 8x8 subtiling.
  if (this->_metadata->fh().version() == 1)
    return false;

  // TODO-Performance: Might implement in-place integral to float conversion.
  if (result_type != this->_metadata->ih().datatype())
    return false;
  const std::array<std::int64_t,3> bs  = this->_metadata->ih().bricksize();
  if ((size[0] != bs[0]) ||
      (size[1] != bs[1]) ||
      (size[2] != bs[2]) ||
      (start[0] % bs[0]) != 0 ||
      (start[1] % bs[1]) != 0 ||
      (start[2] % bs[2]) != 0)
    return false;

  // TODO-Performance: Could be simplified.
  validateUserLod(lod);
  validateUserPosition(start, size, lod); // Could be simplified.

  // TODO-Performance: Could be simplified.
  std::vector<LutInfoEx> bricklist = partsNeeded(start, size, lod);
  if (bricklist.size() != 1)
    throw OpenZGY::Errors::ZgyInternalError("expeditedRead messed up.");
  const LutInfoEx& brick = bricklist.front();

  switch (brick.status) {
  case BrickStatus::Missing:
    if (blogger(2, ""))
      blogger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of missing brick"
              << " at lod=" << lod);
    // TODO-WIP-BrickAPI: The user might have established a different
    // default value. Bricks inside the survey has already been stored
    // with that constant. Reading bricks completely outside the survey
    // ought to return the same value. GenLod might care.
    fillMe(data, size, missingValue(/*as_float=*/false), result_type);
    break;

  case BrickStatus::Constant:
    if (blogger(2, ""))
      blogger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of constant brick"
              << " at lod=" << lod);
    fillMe(data, size, brick.double_constvalue, result_type);
    break;

  case BrickStatus::Normal:
    if (brick.size_in_file != size[0] * size[1] * size[2] *
        static_cast<std::int64_t>(RawDataTypeDetails(result_type).size)) {
      throw OpenZGY::Errors::ZgyInternalError("Bad size in expeditedRead.");
    }
    if (blogger(2, ""))
      blogger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of regular brick"
              << " at lod=" << lod);
    this->_file->xx_read(data, brick.offset_in_file, brick.size_in_file, UsageHint::Data);
    // TODO-Low: If deliverOneBrick() clears padding samples, do likewise here.
    break;

  case BrickStatus::Compressed:
    // TODO-Performance: Also support decompression.
    if (blogger(2, ""))
      blogger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of compressed brick (not implemented)");
    return false;

  default:
    throw OpenZGY::Errors::ZgyInternalError("Internal error, bad brick status");
  }
  return true;
}

/**
 * Mitigate the problem that ZgyInternalBulk is used both for reading
 * and writing, making it difficult to keep track of thread safety.
 * The function below forces the code to be more explicit read vs write.
 *
 * See ZgyInternalBulk::ZgyInternalBulk for more details.
 */
std::shared_ptr<ZgyInternalMeta>
ZgyInternalBulk::get_metadata_rw()
{
  if (!_metadata_rw)
    throw OpenZGY::Errors::ZgyUserError("ZGY file not open for write");
  return _metadata_rw;
}

/**
 * Convenience for invoking _loggerfn with a simple message.
 * This isn't much different from invoking the callback directly.
 * But it makes debugging slightly simpler to have an easy place
 * to set a breakpoint. It also adds more symmetry with respect
 * to the stringstream version, which does add value.
 */
bool
ZgyInternalBulk::blogger(int priority, const std::string& message) const
{
  return _loggerfn(priority, message);
}

/**
 * Convenience for invoking _loggerfn with a stringstream.
 * Due to a somewhat naughty cast, the function can be caller as:
 *
 *   if(blogger(pr1))
 *    blogger(pri, std::stringstream() << some << data << here);
 *
 * The first line is optional. It just prevents the expression in
 * the second line from being evaluated if debugging is disabled.
 *
 * CAVEAT: Modifying a temporary in this manner is not very clean.
 * Newer versions of C++ forbid assigning a non-const reference to
 * a temporary, which can result in some confusing build errors.
 * And it might work for operator<< defined as members of ostream,
 * but not user defined operator<< outside the class. So the result
 * depends on both compiler version and the ordering of the data
 * being output.
 *
 * There is yet another trick to work around the issues:
 * Use std::stringstream().flush() instead of just std::stringstream()
 * when calling the logger. flush() on a new sstream is a no-op, but
 * it returns a non-const reference. This fools the compiler, as it
 * will no longer realize it shouldn't allow that.
 *
 * Stackoverflow has more details: https://stackoverflow.com/\
 * questions/7979215/c-stringstream-to-ostream-to-string
 */
bool
ZgyInternalBulk::blogger(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return blogger(priority, sstream ? sstream->str() : std::string());
}

/**
 * Check that lod is inside the valid bounds. With additional checks
 * for whether this lod can be accessed at this point in time.
 *
 * Note that a file currently being written to but not finalized yet
 * will see nlods() == 1 (which is entirely correct) but to avoid
 * some chicken and egg problems it will still be permitted to read
 * back higher LOD levels. See ZgyInternalBulk::validateUserPosition.
 * When an incomplete file is opened read only, attempts to read low
 * resolution data will throw an exception.
 *
 * This check should be done before calling validateUserPosition().
 * If calling that function in a loop with different positions, the
 * lod only needs to be checked once.
 */
void
ZgyInternalBulk::validateUserLod(int32_t lod) const
{
  if (lod == 0)
    return;
  const bool open_for_write = this->_metadata_rw != nullptr;
  const bool has_lods = open_for_write ? true :
    this->_metadata->cached_hasBrickLOD();
  if (lod > 0 && !open_for_write && !has_lods) {
    blogger(1, "This file does not contain low resolution data.\n");
    throw OpenZGY::Errors::ZgyUserError
      ("This file does not contain low resolution data.");
  }
  const std::int32_t nlods = !has_lods ? 1 :
    static_cast<std::int32_t>(this->_metadata->ih().lodsizes().size());
  // Better error message for the most obvious error.
  if (lod > 0 && !open_for_write && !has_lods) {
    blogger(1, "This file does not contain low resolution data.\n");
    throw OpenZGY::Errors::ZgyUserError
      ("This file does not contain low resolution data.");
  }
  if (lod < 0 || lod >= nlods) {
    std::stringstream ss;
    ss << "Requested lod " << lod <<
      " is outside the valid range 0 to " << nlods - 1 << " inclusive";
    blogger(1, ss.str() + "\n");
    throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
}

/**
 * Check that i,j,k,lod are all inside valid bounds. Throw if not.
 * If used to validate an alpha tile then k should be passed as 0.
 * This is similar to LookupTable::_validatePosition(). The input
 * is in trace numbers, not brick numbers.
 *
 * The supplied lod is trivially checked to make it safe to access
 * ih.lodsizes(). Caller ought to have called validateUserLod *before*
 * calling validateUserPosition(). So, the test here should never fail.
 *
 * Note on test coverage: For defensive coding these checks are done
 * both here at a high level and in lookuptable. This hurts coverage
 * since it might not be possible to get those other checks to fail.
 * And some errors may be caught even earlier e.g. if trying to make
 * a DataBuffer with negative size.
 */
void
ZgyInternalBulk::validateUserPosition(
    const std::array<std::int64_t,3>& start,
    const std::array<std::int64_t,3>& size,
    std::int32_t lod) const
{
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const std::array<std::int64_t,3> one{1,1,1} ;
  const std::array<std::int64_t,3> bs  = this->_metadata->ih().bricksize();
  const std::array<std::int64_t,3> end = ((start + size + bs - one) / bs) * bs;
  if (lod < 0 || lod >= (std::int32_t)ih.lodsizes().size())
    throw OpenZGY::Errors::ZgyUserError("Requested lod is outside range");
  const std::array<std::int64_t,3> ssize = ih.lodsizes()[lod] * bs;
  if (size[0] <= 0 || size[1] <= 0 || size[2] <= 0) {
    blogger(1, "Requested region is empty\n");
    throw OpenZGY::Errors::ZgyUserError("Requested region is empty");
  }
  if (start[0] < 0 || end[0] > ssize[0] ||
      start[1] < 0 || end[1] > ssize[1] ||
      start[2] < 0 || end[2] > ssize[2]) {
    std::stringstream ss;
    ss << "Requested region"
       << " from (" << start[0] << ", " << start[1] << ", " << start[2] << ")"
       << " to (" << end[0] << ", " << end[1] << ", " << end[2] << ")"
       << " lod " << lod
       << " is outside the valid range"
       << " (0, 0, 0)"
       << " to (" << ssize[0] << ", " << ssize[1] << ", " << ssize[2] << ")"
       ;
    blogger(1, ss.str() + "\n");
    throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
}

bool
ZgyInternalBulk::singleBrickOutsideSurvey(
    const std::array<std::int64_t,3>& start,
    const std::array<std::int64_t,3>& size,
    int32_t lod) const
{
  const std::array<std::int64_t,3> bs  =   _metadata->ih().bricksize();
  const std::array<std::int64_t,3> ssize = _metadata->ih().lodsizes()[lod] * bs;
  int inside{0};
  for (int dim = 0; dim < 3; ++dim) {
    if ((start[dim] % bs[dim]) != 0)
      return false;
    if (size[dim] != bs[dim])
      return false;
    if (start[dim] >= 0 && start[dim] + bs[dim] <= ssize[dim])
      ++inside;
  }
  if (inside == 3)
    return false;
  return true;
}

/**
 * Scale integral data from storage according to the coding range.
 * The operation is a no-op if file_dtype is float. The data must
 * be known to be in storage domain, we don't try to guess based on
 * the DataBuffer's valuetype.
 */
std::shared_ptr<DataBuffer>
ZgyInternalBulk::scaleDataToFloat(
      const std::shared_ptr<DataBuffer>& data) const
{
  // The data buffer type should match file_dtype but I won't bother to check.
  // Also, if result ends up empty this means no conversion is needed.
  std::array<double,2> factors = this->_metadata->ih().storagetofloat();
  std::shared_ptr<DataBuffer> result = data->scaleToFloat(factors);
  return result ? result : data;
}

std::shared_ptr<DataBuffer>
ZgyInternalBulk::scaleDataToStorage(
      const std::shared_ptr<DataBuffer>& data) const
{
  std::array<double,2> factors = this->_metadata->ih().storagetofloat();
  std::shared_ptr<DataBuffer> result =
    data->scaleToStorage(factors, _metadata->ih().datatype());
  return result ? result : data;
}

namespace {
  template<typename T> double decodeConstantT(std::uint32_t in)
  {
    // LookupTable::getBrickFilePositionFromIndex() has already
    // handled the v1 special case of 0x01 meaning all-zero. Checking
    // for it would NOT WORK HERE, because the upper bits have already
    // been stripped off. The old "all zero" and the new "all 1" would
    // have looked the same.
    T value;
    byteswapT(&in);
    memcpy(&value, &in, sizeof(T));
    byteswapT(&value);
    return static_cast<double>(value);
  }

  template<typename T> std::uint32_t encodeConstantT(double in)
  {
    if (std::numeric_limits<T>::is_integer)
      in = std::max((double)std::numeric_limits<T>::lowest(),
                    std::min((double)std::numeric_limits<T>::max(),
                             std::round(in)));
    T value = static_cast<T>(in);
    std::uint32_t result{0};
    byteswapT(&value);
    memcpy(&result, &value, sizeof(T));
    byteswapT(&result);
    return result;
  }
}

/**
 * Decode the constant brick value from the re-purposed file offset
 * stored in the lookup table in the ZGY file. Return as a double,
 * which should have enough precision to precisely represent all
 * supported value types.
 *
 * In ZGY v2 and later the actual constant value is stored on the file
 * as little-endian in the first 1, 2, or 4 bytes in the 64-bit slot
 * that normally holds the file offset. The v1 format does not store
 * arbitrary constants. If it had, it would be subject to the swapping
 * of low and high 32-bit words.
 *
 * Interpreting that value in C++ safely and portably is tricky.
 *
 * The input to this method has already been byteswapped to machine
 * byte order and contains the least significant 32 bits of the stored
 * file offset. So the correct handling is to mask out the least
 * significant bytes (however many we need) and return that. It may be
 * simpler to temporarily byte swap back to little endian so the code
 * can copy "first N bytes" instead. Because, memcpy seems to be the
 * safest approach to avoid casts that are technically illegal.
 *
 * The same reasoning applies to the inverse encodeConstantT().
 * The returned result will be byteswapped before being written out,
 * so the constant needs to be stored in the least significant part
 * of the result and not the first N bytes of the result.
 *
 * TODO-Low: Idea: Could I convert a 0x80xxx pointer to the pointer's address
 * and set the size to one sample? This might actually work...
 * The "1" pointer from the v1 format might even be treated the same way;
 * the upper 23 bits of the pointer will be all zero so we could
 * point to that. Just be careful to not overwrite it. And caller
 * may need the not-converted version to test whether this is a constant.
 * Or can it simply test for size == 1 sample?
 *
 * TODO-Low: Idea: passing the pointer to the first byte might actually
 * be easier to deal with.
 */
/*static*/ double
ZgyInternalBulk::decodeConstant(std::uint32_t in, RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::SignedInt8:    return decodeConstantT<std::int8_t>(in);
  case RawDataType::UnsignedInt8:  return decodeConstantT<std::uint8_t>(in);
  case RawDataType::SignedInt16:   return decodeConstantT<std::int16_t>(in);
  case RawDataType::UnsignedInt16: return decodeConstantT<std::uint16_t>(in);
  case RawDataType::SignedInt32:   return decodeConstantT<std::int32_t>(in);
  case RawDataType::UnsignedInt32: return decodeConstantT<std::uint32_t>(in);
  case RawDataType::Float32:       return decodeConstantT<float>(in);
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Inverse of decodeConstant(). See that function for details.
 */
/*static*/ std::int32_t
ZgyInternalBulk::encodeConstant(double in, RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::SignedInt8:    return encodeConstantT<std::int8_t>(in);
  case RawDataType::UnsignedInt8:  return encodeConstantT<std::uint8_t>(in);
  case RawDataType::SignedInt16:   return encodeConstantT<std::int16_t>(in);
  case RawDataType::UnsignedInt16: return encodeConstantT<std::uint16_t>(in);
  case RawDataType::SignedInt32:   return encodeConstantT<std::int32_t>(in);
  case RawDataType::UnsignedInt32: return encodeConstantT<std::uint32_t>(in);
  case RawDataType::Float32:       return encodeConstantT<float>(in);
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Return the sample value to be used where no real data is used. Such as:
 *
 * - A brick inside the survey is missing on read.
 * - A requested brick completely outside the survey is attempted read.
 * - A sample is read from the padding area in an edge brick.
 * - A sample is written to the padding area in an edge brick.
 *
 * Normally, the chosen storage value is chosen heuristically as that
 * sample value which, if converted to float, ends up as close as
 * possible to zero. The function returns that storage value cast to
 * double. In most cases you do want the storage value. Because, if
 * the user asked for integer storage to float user values, this
 * conversion happens at a higher level.
 *
 * Set as_float=true if you really want the value converted here. You
 * might need this if filling a buffer to be used for read, where you
 * know that missing samples will not be delivered, i.e. never touched
 * at all.
 *
 * TODO-WIP-BrickAPI: The default value needs to be added to the
 * physical file format. The current implementation has many caveats.
 *
 * Caveat: The behavior does not match the old ZGY-Public accessor.
 * The behavior in that case was even weirder, so it made no sense to
 * try to replicate.
 *
 * Caveat: Some places might call IInfoHeaderAccess::defaultvalue()
 * directly. If the rules are going to change, it might be better to
 * update that lower level function to make sure all the code uses the
 * same default. See GenLodBase::_wa_defaultstorage in GenLodBase.
 * Which... isn't quite used yet ;-(
 *
 * Caveat: The heuristic might not be appropriate for samples that
 * don't represent seismic.
 *
 * Caveat: If the application fills the newly created survey with its
 * choice of default value (which is the recommended practice), the
 * value is established when the file is written. There will then be
 * no missing bricks inside the survey. If it doesn't, the value is
 * established when the file is read. This difference violates the
 * principle of least surprise.
 *
 * Caveat: If the application fills the newly created survey with its
 * choice of default value, the heuristic default value will still be
 * used when reading a brick completely outside the survey. Which
 * genlod might do. This is not critical. But it might cause low
 * resolution bricks to end up non-constant when they need not be.
 * Additional heuristics in genlod might solve this.
 *
 * Caveat: The value to be used for writing padding samples to an
 * otherwise const brick should be that const value instead of the
 * heuristic default. Already handled (I think) at higher levels.
 *
 * Caveat: In brick mode, after reading is completed, the padding area
 * inside each brick (if any) should be filled with deterministic
 * data. Failure to do this is not critical but might cause
 * heissenbugs. Note, is this already covered? Can probably use
 * setPaddingSamples().
 */
double
ZgyInternalBulk::missingValue(bool as_float) const
{
  return (as_float ?
          this->_metadata->ih().defaultvalue() :
          this->_metadata->ih().defaultstorage());
}

/**
 * Get information about a single brick to be read.
 * status, offset_in_file, size_in_file, raw_constant, double_constvalue.
 * The input "start" wil be copied to survey_position.{i0,j0,k0}.
 */
LutInfoEx
ZgyInternalBulk::makeLutInfoEx(
      const std::array<std::int64_t,3>& start,
      int32_t lod) const
{
  const IInfoHeaderAccess& ih         = this->_metadata->ih();
  const std::array<std::int64_t,3> bs = this->_metadata->ih().bricksize();
  const ILookupTableAccess& blup      = this->_metadata->blup();
  const IJK pos{start[0], start[1], start[2], bs[0], bs[1], bs[2]};
  const LookupTable::LutInfo info = LookupTable::getBrickFilePosition
    (start[0]/bs[0], start[1]/bs[1], start[2]/bs[2], lod,
     /*needed?*/ih.lodsizes(), ih.brickoffsets(), blup.lup(), blup.lupend(),
     ih.bytesperbrick());
  return LutInfoEx(info, pos, decodeConstant(info.raw_constant, ih.datatype()));
}

/**
 * Return the list of bricks needed to cover the entire area given by
 * start and size. The resulting list gives the sample position
 * relative to the start of the survey. To get the brick position you
 * need to divide by brick size. Also compute file offsets etc. from
 * the lookup tables. This will indirectly validate start and size.

 * Implementation note: Don't collect the brick list in one loop and
 * then compute lookup info later. If requested size is ridiculously
 * large this may cause us to run out of memory before getting to the
 * validation stage.
 */
std::vector<LutInfoEx>
ZgyInternalBulk::partsNeeded(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const
{
  std::vector<LutInfoEx> bricklist;
  const std::array<std::int64_t,3> bs = this->_metadata->ih().bricksize();
  const std::array<std::int64_t,3>
    brick0{(start[0] / bs[0]) * bs[0],
           (start[1] / bs[1]) * bs[1],
           (start[2] / bs[2]) * bs[2]};
  const std::array<std::int64_t,3>
    brickN{((start[0] + size[0] - 1) / bs[0]) * bs[0],
           ((start[1] + size[1] - 1) / bs[1]) * bs[1],
           ((start[2] + size[2] - 1) / bs[2]) * bs[2]};
  std::array<std::int64_t,3> iter;
  for (iter[0] = brick0[0]; iter[0] <= brickN[0]; iter[0] += bs[0]) {
    for (iter[1] = brick0[1]; iter[1] <= brickN[1]; iter[1] += bs[1]) {
      for (iter[2] = brick0[2]; iter[2] <= brickN[2]; iter[2] += bs[2]) {
        bricklist.push_back(makeLutInfoEx(iter, lod));
      }
    }
  }
  return bricklist;
}

/**
 * \brief This is the final step in readToExistingBuffer().
 *
 * \details \callgraph \callergraph
 *
 * This is the final step in readToExistingBuffer(). The data has
 * been read from storage, so now it needs to be copied back to
 * the user. This function may be invoked multiple times if data
 * was needed from more than one brick.
 *
 * Arguments:
 *    result   -- user's buffer which was passed to readToExistingBuffer().
 *    start    -- user's supplied position as a 3-tuple of index values.
 *    as_float -- user's choice about converting data.
 *    startpos -- position of the start of this brick.
 *    raw      -- bulk data for this brick as it exists on file.
 *    rawsize  -- size in BYTES of "raw".
 *    brickstatus -- normal, compressed, constant, ...
 *
 * This low level function technically deals with raw bytes as input and
 * and the higher level DataBuffer as output.
 *
 * The function will always be called to deliver one full brick. If the raw
 * data only contains a single sample then this is a constant-value brick.
 * otherwise there should be exactly bricksize-in-bytes available.
 * After any decompression has been done.
 * If the brick size is (1,1,1) (not even possible right now) then the
 * distinction between constant-value and regular goes away.
 *
 * TODO-Worry: beware this kludge: If passing sizeof(double) data then this also
 * represents a constant-value brick where the constant needs to be cast
 * to the correct type. There are some highly theoretical cases where this
 * becomes ambiguous.
 *
 * Caveat: If the data pointer is not correctly aligned to the type of data
 * then the returned DataBuffer will also have a misaligned pointer.
 * Another caveat: Even though result is reference counted, its buffer
 * might end up set to "raw" which means it must not be used after the
 * delivery function returns. Mitigated by adding a clear() member to
 * DataBuffer. That doesn't help if somebody already has a raw pointer though.
 *
 * Thread safety: Yes, but TODO-Worry this is not trivial to verify.
 *
 * Multithreading by having multiple read requests from the API layer is
 * safe. Those requests don't share any mutable data. Specifically,
 * "result" will be different for each read request.
 *
 * Multithreading by having the low level xx_readv() deliver the
 * requested bricks in parallel is also supposed to be thread safe.
 * The result parameter is shared between threads but different
 * threads will be writing to different parts of the buffer due to
 * "bpos" being different. (Note 1,2). The raw pointer might also be
 * referring to shared data. E.g. a cache. (Note 3). But since the
 * data will be destined for different parts of the result buffer it
 * shouldn't be possible for the same bytes to be sent to different
 * threads.
 *
 * Note 1:
 *
 * No two threads should be trying to write to the same location in
 * result. But there are no guarantees of alignment. E.g, thread#1
 * might write the first 7 bytes, thread#2 the next 7 bytes. The
 * assumption is that if the hardware needs to do a read/modify/write
 * due to the memory bus or cache line being wider than 8 bits (which
 * they definitely are) then the hardware or the compiler is
 * responsible for preventing a race condition. For C++11 this appears
 * to be the case. C++11 1.7 [intro.memory]/p2, irrelevant note
 * omitted:
 *
 *    A memory location is either an object of scalar type or a
 *    maximal sequence of adjacent bit-fields all having non-zero
 *    width. Two or more threads of execution (1.10) can update
 *    and access separate memory locations without interfering
 *    with each other.
 *
 * Note (*2):
 *
 * In the lower layers there are situations where the same block can
 * be requested twice in a single request. The reason for that is that
 * the caller needs two regions in the file and padding causes the two
 * to overlap. Once we get here the padding shouldn't be visible.
 *
 * Note (*3):
 *
 * Besides, if xx_readv() was called with immutable_ok=true then there
 * will be no race condition because "raw" will only be read from. And
 * if immutable_ok=false e.g. because of byteswapping then "raw" needs
 * to be a private copy of the data anyway.
 */
void
ZgyInternalBulk::deliverOneBrick(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& startpos,
      const ReadRequest::data_t& raw, std::int64_t rawsize,
      BrickStatus brickstatus, bool as_float) const
{
  if (blogger(5))
    blogger(5, std::stringstream()
            << "  Special delivery! type " << (int)brickstatus
            << " user buffer start "
            <<"(" << start[0] << ", " << start[1] << ", " << start[2] << ")"
            << " this brick start "
            <<"(" << startpos[0] <<
            ", " << startpos[1]
            << ", " << startpos[2] << ")"
            << " got " << rawsize << " bytes.\n");

  const RawDataType dtype = _metadata->ih().datatype();
  const std::array<std::int64_t,3> bricksize = _metadata->ih().bricksize();

  // TODO-Low const-correctness? Might not be worth the trouble.
  std::shared_ptr<DataBuffer> data;

  // TODO-Low this might not be the right place.
  // Note that if the file contained compressed integral data
  // and the application asked for float data back, the decompressor
  // would be responsible for the int to float conversion but NOT the
  // scaling from storage values to real values. Currently this is
  // N/A becase compressed files are always float. There is almost
  // nothing to gain by compressing integral data.
  // The one case I can think of is lossless compression of
  // discrete data such as classification. Not supported but might
  // be supported in the future. In that case retrieving float
  // does not make sense and must be disallowed.

  // TODO-Low, to make it clear what is going on I should probaby
  // rename Compressed to CompressedFloat.

  switch (brickstatus) {
  case BrickStatus::Compressed: {
    rawdata_t rawdata(raw, rawsize);
    rawdata_t tmp = CompressFactoryImpl::decompress(rawdata, brickstatus, bricksize);
    if (!tmp.first || !tmp.second)
      throw OpenZGY::Errors::ZgyFormatError("Compression type not recognized");
    if (tmp.second != bricksize[0]*bricksize[1]*bricksize[2]*(int)sizeof(float))
      throw OpenZGY::Errors::ZgyInternalError("Wrong decompressed size");
    auto floatdata = std::static_pointer_cast<const float>(tmp.first);
    // TODO-Low: Const-correctness of DataBuffer, at least at runtime.
    // Maybe just an overloaded constructor with const arg that will
    // set a readonly flag. And have both data() and rwdata() methods.
    auto floatunconst = std::const_pointer_cast<float>(floatdata);
    data.reset(new DataBufferNd<float,3>(floatunconst, bricksize));
    break;
  }
  case BrickStatus::Normal: {
    // TODO-Low: byteswap here. maybe clone first if not owndata.
    // TODO-Low: Casting away the const is ugly. DataBuffer not const-correct.
    data = DataBuffer::makeDataBuffer3d(std::const_pointer_cast<void>(raw), rawsize, bricksize, dtype);
    if (_metadata->fh().version() == 1) {
      // Rare case. Does not need to be performant.
      if (rawsize != bricksize[0]*bricksize[1]*bricksize[2]*data->itemsize())
        throw OpenZGY::Errors::ZgyInternalError("Wrong size of subtiled data");
      data = data->clone();
      subtiling(bricksize, data->itemsize(), data->voidData().get(), raw.get(), true);
    }
    break;
  }
  case BrickStatus::Missing:
  case BrickStatus::Constant:
    // TODO-Low: byteswap? Or did the constvalue decode do that?
    if (!raw || rawsize != static_cast<std::int64_t>(sizeof(double)))
      throw OpenZGY::Errors::ZgyInternalError("Expected 8 bytes of scalar value");
    data = DataBuffer::makeScalarBuffer3d(*static_cast<const double*>(raw.get()), bricksize, dtype);
    break;
  default:
    throw OpenZGY::Errors::ZgyInternalError("Unknown brickstatus");
  }

  if (as_float && dtype != RawDataType::Float32)
    data = scaleDataToFloat(data);

  // Note, we might pass in the survey range as well to force
  // all padding bytes to be set to the default value. Less
  // surprises for the caller. It may look a bit odd if the user
  // does a flood-fill of the entire survey to a given value and
  // later sees that the content is different in the padding area.
  // But, the caller should ignore the padding.
  // TODO-Worry: Keep this decision in sync with what is done in
  // the shortcut in expeditedRead().
  //
  // On write the padding samples should also be forced to contain
  // the same value. If nothing else, to help compression. But for
  // efficiency reasons the value is not specified. Typically it
  // will be the default "absent" value but if the rest of the
  // brick has a const value and no allocated disk space then
  // they will inherit that constant.

  result->copyFrom(data.get(), startpos.data(), start.data(), nullptr, nullptr);

  // Belt and suspenders mode...
  // The bulk data inside the DataBuffer I created will now go out of scope;
  // it is only valid for the duration of this function. The data smart pointer
  // should also be going out of scope now. But in case it doesn't I will
  // zero out the buffer pointer.

  data->clear();
}

// write support //

/**
 * Compute the size of the used (inside survey) area of a data buffer
 * with size "size". size will probably always be equal to bricksize
 * but the function should still work if it isn't.
 */
std::array<std::int64_t,3>
ZgyInternalBulk::usedPartOfBrick(
      const std::array<std::int64_t,3>& size,
      const std::array<std::int64_t,3>& brickpos_in,
      std::int32_t lod) const
{
  std::array<std::int64_t,3> used{0,0,0};
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  for (int dim=0; dim<3; ++dim) {
    const std::int64_t bricksize     = ih.bricksize()[dim];
    const std::int64_t surveysize0   = ih.size()[dim];
    const std::int64_t brickpos      = brickpos_in[dim];
    const std::int64_t lodfactor     = static_cast<std::int64_t>(1) << lod;
    const std::int64_t surveysizelod = (surveysize0 + lodfactor - 1) / lodfactor;
    const std::int64_t available     = surveysizelod - (brickpos * bricksize);
    used[dim] = std::max((std::int64_t)0, std::min(available, size[dim]));
  }
  return used;
}

/**
 * Return True if all useful samples in this brick have the same value.
 * Padding samples outside the survey are not useful and should not
 * be checked since they could easily have been set to a padding value
 * that doesn't mach the interesting stuff.
 */
bool
ZgyInternalBulk::isUsedPartOfBrickAllConstant(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& brickpos,
      int32_t lod) const
{
  const std::array<std::int64_t,3> used =
    this->usedPartOfBrick(data->size3d(), brickpos, lod);
  return data->isAllSame(used.data());
}

/**
 * Pad unused parts of the data buffer by replicating the last samples,
 * but only up to a multiple of 'modulo' samples. Handles just one
 * dimension, so caller will typically invoke us three times.
 */
/*static*/ void
ZgyInternalBulk::setPaddingToEdge(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      int modulo, int dim)
{
  std::int64_t slice = used[dim];
  while (slice > 0 && slice < data->size3d()[dim] && (slice % modulo) != 0) {
    std::array<std::int64_t,3> srcorig{0,0,0};
    std::array<std::int64_t,3> dstorig{0,0,0};
    std::array<std::int64_t,3> cpyorig{0,0,0};
    std::array<std::int64_t,3> cpysize=data->size3d();
    srcorig[dim] = 1;
    cpyorig[dim] = slice;
    cpysize[dim] = 1;
    // Had I set srcorig[dim] = 0, this would have been a no-op copying one
    // 2d slice, at one past the last valid slice of samples, to itself.
    // By setting srcorig[dim] = 1, I "trick" the copier to read from the
    // previous slice instead. Which is what is needed to get padding.
    // CAVEAT: copyFrom() does not handle overlapping data. As long as the
    // loop only copies one slice at a time, this should not be an issue.
    data->copyFrom(data.get(), srcorig.data(), dstorig.data(), cpyorig.data(), cpysize.data());
    ++slice;
  }
}

/**
 * Pad unused parts of the data buffer with a constant.
 * Handles just one dimension, so caller should invoke us three times.
 */
/*static*/ void
ZgyInternalBulk::setPaddingToConst(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missing_val, int dim)
{
  if (used[dim] > 0 && used[dim] < data->size3d()[dim]) {
    std::array<std::int64_t,3> srcorig{0,0,0};
    std::array<std::int64_t,3> dstorig{0,0,0};
    std::array<std::int64_t,3> cpyorig{0,0,0};
    std::array<std::int64_t,3> cpysize=data->size3d();
    cpyorig[dim] = used[dim];
    cpysize[dim] -= used[dim];

  std::shared_ptr<DataBuffer> filler =
    DataBuffer::makeScalarBuffer3d(missing_val, data->size3d(), data->datatype());
  data->copyFrom(filler.get(), srcorig.data(), dstorig.data(), cpyorig.data(), cpysize.data());
  }
}

/**
 * Make sure the contants of the padding area, if any, is
 * deterministic instead of whatever garbage the buffer holds.
 * The buffer is updated in place.
 *
 * TODO-Low: Different padding might be indicated depending on
 * compression algorithm. Ideally the compressor itself should be
 * allowed to choose. The current implementation is tuned for ZFP.
 * The algorithm also works well for low resolution bricks, whether
 * compressed or not.
 *
 * Possible implementations:
 *
 *   1 - pass the compressor_t into this function, and invoke a new
 *       .pad() method that lets the compressor have full control
 *       over how the padding is done.
 *   2 - When eventually invoking the compressor, pass the "used"
 *       information as the second argument. The compressor can use
 *       the information to add padding, or as a hint to the actual
 *       compression algorithm.
 *   3 - Include "used" as part of the DataBuffer type.
 *
 * Current strategy:
 *
 *   ZFP:
 *      Replicate to make size a multiple of 4 samples, pad with
 *      storage-zero outside that.
 *   Lowres:
 *      Same as ZFP. We only need a multiple of 2 since most LOD
 *      algorithms operate on a 2x2x2 input, so 4 will also work.
 *      Technically the lowpass algorithm could use more padding
 *      vertically but it is unlikely that anybody will notice.
 *   Old ZGY (if we ever support that)
 *      Might work best padding all unused samples with the edge
 *      value. Or perhaps only in the inline direction which for
 *      the old files is the fastest varying.
 */
/*static*/ void
ZgyInternalBulk::setPaddingSamples(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missing_val, const compressor_t& compressor)
{
  // Currently not possible because "compressor" is a
  // std::function which means it has just one entry point.
  //if (compressor)
  //  return compressor.pad(data, used);

  for (int dim=0; dim<3; ++dim)
    setPaddingToConst(data, used, missing_val, dim);

  for (int dim=0; dim<3; ++dim)
    setPaddingToEdge(data, used, 4, dim);
}

/**
 * \brief Write a single brick.
 *
 * \details \callgraph \callergraph
 *
 * Done in prepareWriteOneNormalBrick() which is called first: Apply padding,
 * compression, byte swapping, and convert from DataBuffer to void*.
 * Done here: write the brick to the physical file or to the cloud.
 * Update the lookup table if a new block was allocated.
 *
 * The operation is currently not thread safe. For cloud writes it would
 * not help much to allow multi threading since most writes end up being
 * buffered anyway. Supporting parallel writes in that scenario would
 * require changes at a lower level. Also, SegmentClosedException would
 * get really problematic to handle.
 *
 * For local writes, temporarily dropping a write lock while executing
 * the two xx_write() calls might help. But that requires a lot of testing.
 * Don't do it unless certain there will be a substantial benefit.
 *
 * Args: brickpos, lod, fileoffset, rawdata, brickstatus, align.
 */
void
ZgyInternalBulk::writeOneNormalBrickWithRetry(const WriteNowArgPack& args)
{
  if (blogger(2))
    blogger(2, "writeOneNormalBrickWithRetry(" + args.toString() + ")\n");
  std::int64_t fileoffset = args.fileoffset;
  if (fileoffset) {
    try {
      this->_file->xx_write(args.rawdata.first.get(), fileoffset, args.rawdata.second, UsageHint::Data);
    }
    catch (const OpenZGY::Errors::ZgySegmentIsClosed&) {
      // The update mode doesn't need to be checked again here
      // unless we want to support UpdateMode.NoLeaks which
      // would require raising an exception.
      // TODO-Test: testing will require careful thought.
      if (blogger(3, ""))
        blogger(3, "Abandon brick at " + std::to_string(fileoffset));
      fileoffset = 0; // Write a new brick, abandoning the old one.
    }
  }
  if (!fileoffset) { // New brick or abandoned brick
    fileoffset = this->_file->xx_eof();
    if (blogger(2))
      blogger(2, std::stringstream()
              << "Allocating space at EOF: "
              << std::hex << fileoffset << std::dec << ")\n");
    // Normally fileoffset will be aligned already. The test is needed
    // if updating a file created by the old accessor because that
    // code won't add padding until the first brick is written. When
    // updating ZFP-compressed files or files on the cloud those must
    // have been created by OpenZGY so there is no issue.
    if (args.align && !this->_file->xx_iscloud()) {
      const std::int64_t alignto = this->_metadata->ih().bytesperbrick();
      fileoffset = ((fileoffset + alignto - 1) / alignto) * alignto;
    }
    this->_file->xx_write(args.rawdata.first.get(), fileoffset, args.rawdata.second, UsageHint::Data);
    LookupTable::setBrickFilePosition
      (args.brickpos[0], args.brickpos[1], args.brickpos[2], args.lod,
       LookupTable::LutInfo(args.brickstatus, fileoffset, args.rawdata.second, 0),
       this->_metadata->ih().lodsizes(),
       this->_metadata->ih().brickoffsets(),
       &this->get_metadata_rw()->blup().lup(),
       &this->get_metadata_rw()->blup().lupend());
  }
}

/**
 * \brief padding, compression, convert to void*
 *
 * \details \callgraph \callergraph
 *
 * Process a single brick that is to be written to storage.
 * Input must be a regular brick, not a scalar.
 * In this function:
 *
 *   * Apply padding
 *   * Apply compression
 *   * Apply byte swapping (doesn't work yet)
 *   * Convert DataBuffer to void* + size suitable for the file i/o layer.
 *
 * Thread safety:
 * The function is thread safe in the sense that multiple bricks
 * belong to the same write request may be processed in parallel.
 * No other changes, such as data being written to the file, are
 * allowed while this is going on. Nor can there be a separate
 * write request initiated by the user because (among other reasons)
 * the bricks might overlap.

 * The next step, actually writing the file, needs to be serialized
 * both because it may want to write at EOF and because even if the
 * file offset is explicitly known the lower level write might not be
 * thread safe. Also when serializing the order of bricks should be
 * preserved. Otherwise performance on read might suffer.
 *
 * Args: brickpos, lod, data, compressor, fileoffset.
 */
std::shared_ptr<const WriteNowArgPack>
ZgyInternalBulk::prepareWriteOneNormalBrick(const WriteBrickArgPack& args)
{
  if (blogger(2))
    blogger(2, "prepareWriteOneNormalBrick(" + args.toString() + ")\n");

  if (args.data->isScalar())
    throw OpenZGY::Errors::ZgyInternalError("Wrong brick type in prepareWriteOneNormalBrick");
  if (args.data->size3d() != this->_metadata->ih().bricksize())
    throw OpenZGY::Errors::ZgyInternalError("Wrong brick size in prepareWriteOneNormalBrick");

  std::array<std::int64_t,3> used =
    this->usedPartOfBrick(args.data->size3d(), args.brickpos, args.lod);
  this->setPaddingSamples(args.data,used,missingValue(false),args.compressor);

  // rawdata_t::first is a smart pointer and it needs to be in case it
  // is later replaced with a buffer of compressed data. Besides, it is
  // convenient if we decide to queue the write request.
  rawdata_t rawdata(args.data->voidData(), args.data->allocsize() * args.data->itemsize());

  BrickStatus brickstatus = BrickStatus::Normal;
  if (args.compressor) {
    // TODO-WIP-BrickedAPI: Skip compression if any NaN is present,
    // ZFP doesn't like it. Also add unit test for this.
    blogger(5, "Attempting compression");
    rawdata_t cdata = args.compressor(rawdata, args.data->size3d());
    if (cdata.first) {
      blogger(5, "compression successful");
      rawdata = cdata;
      brickstatus = BrickStatus::Compressed;
    }
  }
  if (brickstatus == BrickStatus::Compressed) {
    // TODO-Low, isn't it Normal data that needs byte swap?
    // TODO-Low, need to implement BYTE SWAPPING.
    // TODO-High, in shortcut mode we might not own the buffer.
    // It might be safer to unconditionally copy the data.
    //data->byteswap();
  }
  // Arguments that our caller needs to be pass to writeOneNormalBrickWithRetry()
  // All normal bricks are flagged as needing alignment if on-prem.
  // Technically I could have passed !compressor so that uncompressed
  // bricks that were attempted compressed but didn't make it will not
  // need to be aligned. But that case is too obscure to warrant extra
  // testing.
  return std::make_shared<const WriteNowArgPack>
    (args.brickpos, args.lod, args.fileoffset, rawdata, brickstatus,
     brickstatus == BrickStatus::Normal);
}

/**
 * \brief Write one constent-value brick to the file header.
 *
 * \details \callgraph \callergraph
 *
 * Apply conversion of the constant value to a lookup table entry
 * Pass on the write request to the function that updates the lut.
 *
 * The size stored in "data" is ignored. The constant applies to
 * the whole brick. Also, args.fileoffset and args.compressor are N/A.
 *
 * Args: brickpos, lod, data, compressor(N/A), fileoffset(N/A).
 */
void
ZgyInternalBulk::writeOneConstantBrick(const WriteBrickArgPack& args)
{
  if (blogger(2))
    blogger(2, "writeOneConstantBrick(" + args.toString() + ")\n");
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  ILookupTableAccess& blup = this->get_metadata_rw()->blup();
  LookupTable::LutInfo info
    (BrickStatus::Constant, /*offset=*/0, /*size=*/0,
     encodeConstant(args.data->scalarAsDouble(), ih.datatype()));
  LookupTable::setBrickFilePosition
    (args.brickpos[0], args.brickpos[1], args.brickpos[2], args.lod, info,
     ih.lodsizes(), ih.brickoffsets(), &blup.lup(), &blup.lupend());
}

/**
 * Set every single brick in every lod level to the same constant
 * value. If force is true, already allocated bricks will be converted
 * back to all-const. THIS WILL LEAK THOSE ALLOCATED BRICKS.
 * Which might not be a big deal for cloud storage. Which would
 * leak the space anyway if the segment is closed.
 *
 * TODO-WIP-BrickedAPI: BUG ****** if force=false, any brick already
 * inflated will not be touched. Need to write such bricks normally.
 * Pick some code from prepareWriteOneBrick() and prepareWriteOneNormalBrick(),
 * then call writeOneNormalBrickWithRetry(). Or dumb down the code so it won't
 * try to recognize a fill-entire-survey if some bricks have already
 * been allocated.
 *
 * It might be possible to instead re-use some of the code for the
 * normal case. Tricky, because here we want to clear all the lowres
 * data as well.
 */
std::int64_t
ZgyInternalBulk::writeAllConstantBrick(double value, bool force)
{
  if (!force)
    throw OpenZGY::Errors::ZgyInternalError
      ("writeAllConstantBrick(value, false) is not implemented yet");

  std::int64_t leaks{0};
  const IInfoHeaderAccess& ih = this->_metadata->ih();

  // TODO: Verify that brick indices, not traces/samples, are expected.
  index3_t survey_size
    {((ih.size()[0] - 1) / ih.bricksize()[0]) + 1,
     ((ih.size()[1] - 1) / ih.bricksize()[1]) + 1,
     ((ih.size()[2] - 1) / ih.bricksize()[2]) + 1};

  ILookupTableAccess& blup = this->get_metadata_rw()->blup();
  const LookupTable::LutInfo info
    (BrickStatus::Constant, /*offset=*/0, /*size=*/0,
     encodeConstant(value, ih.datatype()));

  for (std::int32_t lod = 0; lod < ih.nlods(); ++lod) {
    for (std::int64_t ii = 0; ii < survey_size[0]; ++ii) {
      for (std::int64_t jj = 0; jj < survey_size[1]; ++jj) {
        for (std::int64_t kk = 0; kk < survey_size[2]; ++kk) {
          const LookupTable::LutInfo oldinfo =
            LookupTable::getBrickFilePosition
            (ii, jj, kk, lod,
             ih.lodsizes(), ih.brickoffsets(), blup.lup(), blup.lupend(),
             ih.bytesperbrick());
          const bool leak =
            (oldinfo.status != BrickStatus::Missing &&
             oldinfo.status != BrickStatus::Constant);
          leaks += (leak ? 1 : 0);
          if (force || !leak) {
            LookupTable::setBrickFilePosition
              (ii, jj, kk, lod, info,
               ih.lodsizes(), ih.brickoffsets(), &blup.lup(), &blup.lupend());
          }
        }
      }
    }
    for (int dim = 0; dim < 3; ++dim)
      survey_size[dim] = (survey_size[dim] + 1) / 2;
  }
  if (blogger(2))
    blogger(2, std::stringstream().flush()
            << "writeAllConstantBrick(" << value
            << ", force=" << std::boolalpha << force << ")"
            << " -> " << leaks << (force ? " leaks" : " inflated"));
  return leaks;
}

/**
 * Return True if this block needs to be leaked by pretending
 * that its block offset has not been allocated yet.
 * Raise an exception if the update is disallowed by _update_mode.
 *
 * Note, in "Always" mode I am only supposed to leak the block
 * if the new one is larger. But that is too much work to
 * figure out here. So, treat "Always" as if it were "Pedantic".
 * If I assume compressed blocks are never larger than uncompressed
 * (which is currently enforced) then I could avoid leaking the
 * brick completely. But that logic is fragile and the rule might not
 * be enforced in the future. And the issue might be academic because
 * compressed files should probably be written in Constant mode.
 *
 * Note, there is one other place where bricks might leak:
 * in writeOneNormalBrickWithRetry() when a ZgySegmentIsClosed is caught.
 * The update mode doesn't need to be checked again at that point
 * unless we want to support UpdateMode.NoLeaks.
 *
 * "pos" and "lod" are only used for logging and error reporting.
 */
bool
ZgyInternalBulk::mustLeakOldBrick(
      const std::shared_ptr<DataBuffer>& data,
      const compressor_t& compressor,
      BrickStatus brickstatus,
      std::array<std::int64_t,3>& pos,
      std::int32_t lod) const
{
  bool error = false;
  bool leak = false;
  if (brickstatus != BrickStatus::Missing) {
    // Data has already been written somehow.
    if (this->_update_mode == UpdateMode::Never)
      error = true;
  }
  if (brickstatus == BrickStatus::Normal) {
    // Data has already been written uncompressed.
    if (this->_update_mode == UpdateMode::Never ||
        this->_update_mode == UpdateMode::Constant)
      error = true;
  }
  if (brickstatus == BrickStatus::Compressed ||
      (brickstatus == BrickStatus::Normal && compressor)) {
    // Data already written, was and/or will become compressed.
    //
    // Note: Technically it would be possible to re-use the
    // old brick if it happens to be big enough. DON'T DO THAT.
    // It gets way too complicated. E.g. after the update the
    // lookup table might need to change from uncompressed to
    // compressed without changing the offset. There are
    // probably other problems as well.
    //
    // TODO-Worry: Can I assume that the caller will veto the
    // compression if an uncompressed brick already exists on the
    // file? If so then I don't need to complain about or leak the
    // "will become compressed" case because it won't.
    leak = true;
    if (this->_update_mode != UpdateMode::Pedantic &&
        this->_update_mode != UpdateMode::Always)
      error = true;
  }
  if (error || blogger(3, "")) {
    std::string bstatus;
    switch (brickstatus) {
    case BrickStatus::Missing    : bstatus = "missing";    break;
    case BrickStatus::Constant   : bstatus = "constant";   break;
    case BrickStatus::Normal     : bstatus = "normal";     break;
    case BrickStatus::Compressed : bstatus = "compressed"; break;
    default:
      bstatus = "BrickStatus::" + std::to_string((int)brickstatus);
      break;
    }
    std::string umode;
    switch (this->_update_mode) {
    case UpdateMode::Never      : umode = "never";      break;
    case UpdateMode::Constant   : umode = "constant";   break;
    //case UpdateMode::NoCompress : umode = "NoCompress"; break;
    //case UpdateMode::NoLeaks    : umode = "NoLeaks";    break;
    case UpdateMode::Always     : umode = "always";     break;
    case UpdateMode::Pedantic   : umode = "pedantic";   break;
    default:
      umode = "UpdateMode::" + std::to_string((int)this->_update_mode);
      break;
    }
    std::stringstream ss;
    ss << "Updating a " << (this->_file->xx_iscloud() ? "cloud" : "local")
       << " " << bstatus << " brick with "
       << (compressor ? "compressed" : data->isScalar() ? "constant" : "normal")
       << " data in update mode " << umode;
    ss << (error ? " is illegal" : " is ok");
    ss << " at pos (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")"
       << " lod " << lod;
    if (blogger(1, ""))
      blogger(1, ss.str());
    if (error)
      throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
  return leak;
}

/**
 * Scalar vs regular, update vs. allocate, veto compression.
 *
 * \details \callgraph \callergraph
 *
 * Process a single brick that is to be written to storage.
 * In this function:
 *
 *   * Might convert from a scalar to a regular buffer or vice versa.
 *   * Might decide whether to update an existing brick, and where.
 *   * Might decide to veto the compression.
 *
 * The veto and the potential convert to scalar is why compression
 * should not be done earlier.
 *
 * brickpos is given relative to this lod level. For lod 0 the valid
 * range is the survey size in bricks. For lod 1 it is half that,
 * rounded upwards.
 *
 * The data must already have been scaled to storage values and
 * converted to the correct value type. As well as split to bricks.
 * The data can be either a scalar or a regular buffer.
 *
 * Thread safety:
 * The function is thread safe in the sense that multiple bricks
 * belong to the same write request may be processed in parallel.
 * No other changes, such as data being written to the file, are
 * allowed while this is going on. Nor can there be a separate
 * write request initiated by the user because (among other reasons)
 * the bricks might overlap.
 *
 * Args: brickpos, lod, data, compressor, fileoffset.
 */
std::shared_ptr<const WriteBrickArgPack>
ZgyInternalBulk::prepareWriteOneBrick(const WriteBrickArgPack& args_in)
{
  WriteBrickArgPack args(args_in);
  // brickpos and lod will not change.
  // data, compressor, fileoffset might all change but usually don't.

  // TODO-Low note copy/pasted code from partsNeeded().
  // On a more general note, would it have been better to
  // call partsNeeded() early on, and have this function
  // and a few others operate on lists? This would be more
  // complicated but might allow better parallelization.
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const ILookupTableAccess& blup = this->_metadata->blup();
  LookupTable::LutInfo info = LookupTable::getBrickFilePosition
    (args.brickpos[0], args.brickpos[1], args.brickpos[2], args.lod,
     ih.lodsizes(), ih.brickoffsets(), blup.lup(), blup.lupend(),
     ih.bytesperbrick());

  // We can easily pretend the brick did not exist and needs to be allocated.
  if (mustLeakOldBrick(args.data, args.compressor, info.status, args.brickpos, args.lod))
    info = LookupTable::LutInfo(BrickStatus::Missing, 0, 0, 0);

  std::array<std::int64_t,3> bs = ih.bricksize();
  const bool data_const = args.data->isScalar();
  const bool file_const = (info.status != BrickStatus::Normal &&
                           info.status != BrickStatus::Compressed);

  if (file_const) { // also true if never written yet.
    if (data_const) {
      // Caller asked to store a constant value.
      blogger(2, "Explicit store constant.\n");
    }
    else if (isUsedPartOfBrickAllConstant(args.data, args.brickpos, args.lod)) {
      // Caller did not explicitly ask for a constant value,
      // but all useable samples (i.e. all samples that are
      // inside the survey boundaries) have the same value.
      blogger(2, "Implicit store constant.\n");
      // Need to change to a real scalar. This is a bit roundabout.
      // TODO-Low can I have a DataBuffer::makeSimilar() that retains
      // data type and NDim but can turn it into a scalar or
      // create a deep copy or maybe an uninitialized copy.
      // Or the other way, create an inflated brick from a constant.
      double value = args.data->scalarAsDouble();
      args.data = DataBuffer::makeScalarBuffer3d(value, bs, args.data->datatype());
    }
    else {
      blogger(2, "Store a new brick.\n");
      args.fileoffset = 0;  // Says to allocate a new block on file.
    }
  }
  else { // !file_const
    if (data_const) {
      blogger(1, "Const data expanded before update of brick.\n");
      // The brick has already been allocated. Cannot set it to
      // constant-value because the file storage for this brick
      // would then be leaked. Cannot compress it either.
      // The assumption is that the old brick was uncompressed,
      // because if it had been compressed then mustLeakOldBrick()
      // would have returned true and we wouldn't have gotten here.
      //
      // TODO-Medium: Fragile code:
      //
      //   * mustLeakOldBrick() might change to allow overwriting
      //     a compressed brick with a smaller one. If so, must
      //     *not* veto the compression.
      //
      //   * mustLeakOldBrick() will always, if compressor is set,
      //     make sure we don't get here. Maybe that is wrong.
      //     Currently resetting the compressor below won't be needed.
      //
      //   * A more robust approach may be to veto the compression
      //     if the old brick is uncompressed. And throw an error
      //     if the old brick is compressed and the new one won't be.
      //     Because we don't update the brick status. Ouch. Maybe
      //     that would have been a good idea?
      //
      double value = args.data->scalarAsDouble();
      args.data = DataBuffer::makeNewBuffer3d(bs, args.data->datatype());
      args.data->fill(value);
      args.fileoffset = info.offset_in_file; // Re-use if possible.
      args.compressor = compressor_t();      // No compress on update.
    }
    else {
      blogger(1, "Update a brick.\n");
      args.fileoffset = info.offset_in_file; // Re-use if possible.
      args.compressor = compressor_t();      // No compress on update.
    }
  }
  // Next, prepareWriteOneNormalBrick or writeOneConstantBrick
  // (depending on the type that args.data has now)
  // should be called with the modified argument pack.
  // TODO-Low: Don't copy args so often.
  return std::make_shared<const WriteBrickArgPack>(args);
}

#if 0
/**
 * Add one level of compression noise.
 *
 * This might not be needed if compressed files never allow updating
 * bricks. The statistics and histogram can then always be made
 * from original data without compression noise. Otherwise, one level
 * of compression noise needs to be added before stats & histo is
 * calculated. Because this is what will be seen when it is time to
 * subtract the previous data. Long term, might even consider storing
 * stats and histo for each brick?
 *
 * If it does turn out to be needed, it is possible that the already
 * compressed buffer is available. So only decompress is needed here.
 *
 * If r/m/w processing is needed then the old data to be kept must
 * already have been copied to the new buffer. Preferably the brick
 * should also have been written at this point. So the operation
 * won't need a scratch buffer.
 */
void
ZgyInternalBulk::addCompressionNoise(
     std::shared_ptr<DataBuffer> data,
     compressor_t compressor)
{
  if (compressor) {
    std::int64_t usize = data->allocsize() * data->itemsize();
    rawdata_t rawdata(data->voidData(), usize);
    rawdata_t cdata = compressor(rawdata, data->size3d());
    if (cdata.first) {
      rawdata_t ddata = CompressFactoryImpl::decompress
        (cdata, BrickStatus::Compressed, data->size3d());
      if (ddata.second == usize) {
        memcpy(data->voidData().get(), ddata.first.get(), ddata.second);
      }
      else {
        throw OpenZGY::Errors::ZgyInternalError("Compression roundtrip size.");
      }
    }
    else {
      // Brick not compressible. That is ok; it will be stored that way also.
    }
  }
}
#endif

/**
 * Given two ranges specified as start and size, return the overlap.
 * If there is no overlap, return all zeros. Caller can test for no overlap
 * by looking at the size in just one dimension.
 */
std::pair<index3_t,index3_t>
ZgyInternalBulk::clip(
     const std::pair<index3_t,index3_t>& a,
     const std::pair<index3_t,index3_t>& b)
{
  std::pair<index3_t,index3_t> result;
  for (int dim = 0; dim < 3; ++dim) {
    result.first[dim]  = std::max(a.first[dim], b.first[dim]);
    result.second[dim] = (std::min((a.first[dim] + a.second[dim]),
                                   (b.first[dim] + b.second[dim])) -
                          result.first[dim]);
    if (result.second[dim] <= 0)
      return std::pair<index3_t,index3_t>{{0,0,0},{0,0,0}};
  }
  return result;
}

/**
 * \brief Handle r/m/w and edge padding if needed.
 *
 * \details \callgraph \callergraph
 *
 * Handle r/m/w and edge padding if needed.
 * The returned brick might be old_brick, new_brick, or a new buffer.
 *
 * The input old_brick is assumed to belong to us, and can be used for
 * whatever purpose we want. The new_brick will not be touched.
 *
 * TODO-WIP-BrickedAPI: If immutable_buffer=false, this might save
 * time in the padding-only case where we can drop one buffer copy and
 * in some cases also a buffer allocate.
 *
 * TODO-WIP-BrickedAPI: Edge padding is also in prepareWriteOneNormalBrick().
 * It isn't needed both places.
 *
 * TODO-WIP-BrickedAPI: Allocate from a buffer pool where needed.
 *
 * Thread safety:
 * Safe because the method is static and does not access shared state.
 * Unless, of course, the caller tries to process the same brick more
 * than once. This applies both to shallow and deep copies.
 */
std::shared_ptr<DataBuffer> /*static*/
ZgyInternalBulk::doReadModifyWrite(
     std::shared_ptr<DataBuffer> old_brick,
     std::shared_ptr<DataBuffer> new_brick,
     const std::pair<index3_t,index3_t>& brick_range,
     const std::pair<index3_t,index3_t>& inside_survey,
     const std::pair<index3_t,index3_t>& inside_user,
     double missing_val,
     const compressor_t& compressor,
     bool immutable_buffer)
{
  if (old_brick &&
      (old_brick->size3d() != new_brick->size3d() ||
       old_brick->datatype() != new_brick->datatype()))
    throw OpenZGY::Errors::ZgyInternalError("Incompatible buffers r/m/w");
  if (inside_survey != inside_user) {
    // r/m/w and possibly edge padding is needed.
    // No special handling needed for scalar new_brick; copyFrom()
    // understands this. A scalar old_brick might need to be inflated.
    if (!old_brick)
      throw OpenZGY::Errors::ZgyInternalError("r/m/w for lod > 0.");
    if (old_brick->isScalar()) {
      // Inflate the brick.
      const double value = old_brick->scalarAsDouble();
      old_brick = DataBuffer::makeNewBuffer3d
        (old_brick->size3d(), old_brick->datatype());
      old_brick->fill(value);
    }
    // Copy the samples to be changed from new_brick to old_brick
    old_brick->copyFrom(new_brick.get(),
                        brick_range.first.data(), brick_range.first.data(),
                        inside_user.first.data(), inside_user.second.data());
    if (brick_range != inside_survey) {
      // Set predictable and easily compressible value in padding area.
      setPaddingSamples
        (old_brick, inside_survey.second, missing_val, compressor);
    }
    return old_brick;
  }
  else if (brick_range != inside_survey && !new_brick->isScalar()) {
    // Edge padding but no r/m/w needed.
    if (!old_brick || old_brick->isScalar()) {
      // Allocate scratch buffer since old_brick cannot be used.
      // In case old_brick is missing, use size & vt of the new brick.
      // They should be identical.
      old_brick = DataBuffer::makeNewBuffer3d
        (new_brick->size3d(), new_brick->datatype());
    }
    memcpy(old_brick->voidData().get(), new_brick->voidData().get(),
           old_brick->totalsize() * old_brick->itemsize());
    setPaddingSamples
      (old_brick, inside_survey.second, missing_val, compressor);
    return old_brick;
  }
  else {
    return new_brick;
  }
}

/**
 * \brief Write a list of bricks to storage.
 *
 * \details \callgraph \callergraph
 *
 * Similar to ZgyInternalBulk::_writeAlignedRegion(), but with a list of bricks.
 *
 * Expects a list of bricks, already converted to storage.
 * Sizes refer to the specified lod level. At lod level N the valid
 * range for size is half (rounded up) of that in lod N-1.

 * For each brick to be written, 3 rectangles are provided or implied.
 * All are in trace/sample numbers at the provided lod.
 *
 *  - brick_range  -- what this memory buffer corresponds to on disk.
 *  - survey_range -- Outside this, fill with predictable value.
 *  - user_clip    -- Samples outside retain their old value.
 *
 * Rules for this function:
 *  - brick_range.start is brick aligned.
 *  - brick_range.size is the file's brick size, found in the file's metadata.
 *  - survey_range.start is implicitly {0,0,0}, and size is in metadata.
 *  - survey_range and user_clip are quietly clipped to brick_range.
 *  - There is just one user_clip, applying to every brick.
 *
 * The "predicatble value" in the padding area is the file's defaultvalue
 * (as close to zero as possible) except that if the user writes a scalar
 * to an edge brick then that same scalar will be used for padding.
 *
 * TODO-WIP-BrickedAPI: The check that deflates a buffer filled with a
 * single value and converts that to a scalar buffer should already
 * have been made.
 *
 * The old contents on disk will always need to be read. This is part
 * of the cost of not deferring all low resolution compute until the
 * end. If each brick is only written once (the most common case),
 * this will be cheap because scalar buffers will be returned and no
 * scratch buffers will be needed. Otherwise the missing buffer pool
 * becomes a performance issue. Buffers get allocated and destroyed in
 * each step.
 *
 * TODO-WIP-BrickedAPI: The function needs to calculate changes to
 * statistics and report dirty state to the module that keeps track of
 * when lowres bricks are ready to be written. If user_start/size
 * caused additional clipping, assume that there will be additional
 * writes to the same bricks. This affects when lowres bricks are
 * ready.
 *
 * TODO-WIP-BrickedAPI: There is some code here that facilitates r/m/w
 * when called from the general API. This smells of feature envy, but
 * might help reduce buffer copies.
 *
 * When updating statistics, do the work in a temporary
 * instance and then merge the result while protected by a lock. If
 * the histogram is fixed width (8-bit data and maybe also 16-bit)
 * then it is possible to update the per-file setting directly. But
 * the lock then needs to be held longer.
 *
 * TODO-WIP-BrickedAPI: A minor improvement would be to allocate all
 * the temporary HistogramBuilder outside the parallel loop, then
 * merge them when the parallel region ends. Or even better: keep a
 * pool of builders owned by the open file and don't merge them until
 * finalize. Unfortunately there is the chicken and egg problem with
 * the WeightedAverage algorithm. So, level 2 and above would need to
 * be deferred to finalize. Or at least to a point where a certain
 * percentage of non-constant bricks have been processed.
 *
 * TODO-WIP-BrickedAPI: If there are more bricks than threads,
 * temporaries such as a HistogramBuilder and scratch buffers will be
 * allocated one per brick instead of just one per thread. This can be
 * fixed, but complicates the code. There is another plan that might
 * render this unnecessary: Implement a pool of brick-sized buffers.
 * And another for HistogramBuilder as detailed above. That change is
 * not trivial.
 *
 * TODO-WIP-BrickedAPI: If we can trust immutable_buffer=false then we
 * might be able to avoid a buffer copy. And in some cases a buffer
 * allocation as well. This happens when survey padding is needed and
 * user clipping is not needed. With even more effort there might be
 * more that can be done. E.g. with a reshape implementation that
 * copies data outside a particular region instead of inside. But this
 * is starting to sound ridiculous.
 *
 * TODO-WIP-BrickedAPI: If the application writes unaligned data then
 * the code here might end up doing one more buffer copy that strictly
 * necessary. But applications where performance is critical won't be
 * doing that in the first place.
 *
 * TODO-WIP-BrickedAPI: Unrelated: If the application writes small
 * unaligned chunks where no write request includes full bricks, this
 * will effectively move all lowres processing back to finalize().
 *
 * TODO-WIP-BrickedAPI: Unrelated: If the first few brick columns are
 * dead, the first few lod2 bricks suffer from the chicken and egg
 * problem. BUT if the inputs are all constant then they end up
 * correct anyway. Maybe.
 *
 * TODO-WIP-BrickedAPI: If there are situations where updating a
 * compressed brick is allowed, then statistics and histogram needs
 * to be computed as if the lod0 bricks were read from disk and
 * decompressed. Otherwise it won't be possible to subtract the old
 * statistics correctly. The code can either do an actual read, or
 * add noise by compressing and decompressing in memory.
 *
 * TODO-WIP-BrickedAPI: Allocating buffers or histogram builders
 * inside the parallelFor will be done once per brick instead of once
 * per thread. Mitigate with a double loop. Or implement buffer pools.
 *
 * TODO-WIP-BrickedAPI: Detect, or have caller tell us, whether this
 * write sets the entire survey lod 0 to a const value. See reset_all
 * in the code below. Several shortcuts are desirable:
 *
 * - Skip reading old contents.
 * - Skip subtracting old and adding new statistics.
 * - Skip r/m/w processing in lod 0.
 * - Optionally skip r/m/w processing in lod > 0, might be messy.
 * - Run single threaded.
 * - Do a recursive call to reset all lod levels to the same constant.
 *
 * It might be better to write a new function to handle this case.
 *
 * The new code will need to invoke newTrackedBricksSetAllConst() to
 * reset statistics and histogram, and some other function to reset
 * the dirty bits to become all clean.
 *
 * CAVEAT: Remember that the file might contain already allocated
 * bricks. And some of these might reside in closed segments on the
 * cloud. Had it not been for this, the handling would have been a lot
 * simpler. Is it feasible to truncate the file? Or maybe delete and
 * re-create it? Or should the code forego the short cuts if written
 * bricks exist? Or just leak the already written bricks? Leaking
 * might be better for cloud access. Because that could leak in any
 * case. Or reclaim the bricks bypassing most of the existing code,
 * filling all the allocated bricks with the const value, then assign
 * to the first /n/ bricks? All of these suggestions might violate the
 * principle of least surprise. And add to the testing effort.
 */
void
ZgyInternalBulk::writeAlignedBrickList(
     const std::vector<std::shared_ptr<DataBuffer>>& new_contents,
     const std::vector<index3_t>& position,
     const std::pair<index3_t,index3_t>& user_clip_in,
     int lod, bool is_storage, const compressor_t& compressor,
     bool immutable_buffer)
{
  if (blogger(2))
    blogger(2, std::stringstream().flush()
            << "writeAlignedBrickList " << position.size() << " bricks at lod " << lod);
  if (new_contents.empty())
    return;

  SimpleTimerEx t1(_timers.mt);
  SimpleTimerEx t01(_timers.mt1);

  typedef std::pair<index3_t,index3_t> range_t;
  const std::int64_t lodfactor  = static_cast<std::int64_t>(1) << lod;
  const IInfoHeaderAccess& ih   = this->_metadata->ih();
  const index3_t     bs          = ih.bricksize();
  const index3_t     survey_size = (ih.size() + (lodfactor - 1)) / lodfactor;
  const range_t      survey_range  {{0,0,0}, survey_size};
  const std::int64_t worksize    = new_contents.size();

  const std::pair<index3_t,index3_t>& user_clip =
    user_clip_in.second[0] != 0 ? user_clip_in : survey_range;
  std::vector<index3_t> brickposlist;
  for (const index3_t& pos : position)
    brickposlist.push_back(index3_t{pos[0]/bs[0],pos[1]/bs[1],pos[2]/bs[2]});
  std::vector<int> rmwlist(brickposlist.size(), -1);

  if (new_contents.size() != position.size())
    throw OpenZGY::Errors::ZgyInternalError
      ("writeAlignedBrickList: Inconsistent size of lists.");
  for (const auto& pos : position)
    if ((pos[0]%bs[0]) != 0 || (pos[1]%bs[1]) != 0 || (pos[2]%bs[2]) != 0)
      throw OpenZGY::Errors::ZgyInternalError
        ("writeAlignedBrickList: Misaligned write.");
  for (const auto& buf : new_contents)
    if (buf->size3d() != bs)
      throw OpenZGY::Errors::ZgyInternalError
        ("writeAlignedBrickList: Buffer size is not a single brick.");

  // Possibly not very useful to check this here. If the application
  // wants to reset the survey to a constant value, it ought to use
  // writeconst(), not writebricks() which will require sending a list
  // of all bricks in the survey. writeconst() should have its own
  // check for all constant.
  if (lod == 0) {
    const std::pair<bool, double> is_reset_survey =
      isResettingEntireSurvey(position, survey_size, bs, new_contents);
    if (is_reset_survey.first) {
      writeAllConstantSurvey(is_reset_survey.second);
      return;
    }
  }

  t01.stop();
  SimpleTimerEx t02(_timers.mt2);

  // The old contents of the bricks are always needed when statistics
  // and histogram is collected during write.
  // For lod > 0 there are no statistics to update, and genlod
  // (which is the only user that writes higher lods) always writes
  // full bricks. So, reading old can be skipped there.
  typedef std::vector<std::shared_ptr<DataBuffer>> buffers_t;
  buffers_t old_contents = (lod > 0) ? buffers_t{} : readBricksToNewBuffers
    (position, buffers_t{}, lod, /*as_float=*/false, /*check_constant=*/true);

  // Actual writes are handled outside the parallel region.
  std::vector<std::shared_ptr<const WriteBrickArgPack>> const_queue(worksize);
  std::vector<std::shared_ptr<const WriteNowArgPack>>   normal_queue(worksize);

  // Delta statistics and histogram accumulated in a temporary variable in
  // case prepareWriteOneBrick throws "not allowed to update compressed brick".
  // Make an empty builder with identical histogram range and isfixed.
  // This still needs to be locked when accessed insted the parallel loop.
  // Might as well use the same mutex that locks *this.
  std::shared_ptr<StatisticData> delta_stats{};
  std::shared_ptr<HistogramData> delta_histo{};
  {
    // TODO-WIP-BrickedAPI: This is seriously convoluted!
    // Don't create a new builder here, or if possible, create only a
    // builder and pass that to trackChanges.
    // TODO-WIP-BrickedAPI: Maybe make a copy of the builder as a
    // template for trackChanges to use so it won't need to lock. This
    // also makes the result more predictable, but might also cause
    // the histogram to get expanded many more times.
    std::shared_ptr<ExpandableBuilder> delta{};
    std::lock_guard<std::mutex> lk(_new_mutex);
    delta = std::make_shared<ExpandableBuilder>
      (*_new_modified_histo, StatisticData(), /*copy=*/false);
    delta_stats = std::make_shared<StatisticData>();
    delta_histo = std::make_shared<HistogramData>(delta->gethisto());
  }

  t02.stop();
  SimpleTimerEx t03(_timers.mt3);
  WorkOrderRunner::parallelFor(worksize, [&](std::int64_t ix)
    {
      SimpleTimerEx t3(_timers.ptimer);
      const range_t brick_range   {position[ix], bs};
      const range_t inside_survey = clip(brick_range, survey_range);
      const range_t inside_user   = clip(inside_survey, user_clip);
      const range_t range_inside_brick // shift to brick-relative
        {{inside_user.first[0] - brick_range.first[0],
          inside_user.first[1] - brick_range.first[1],
          inside_user.first[2] - brick_range.first[2]},
         inside_user.second
        };
      std::shared_ptr<DataBuffer> new_brick = new_contents[ix];
      std::shared_ptr<DataBuffer> old_brick;
      if (!old_contents.empty())
        old_brick.swap(old_contents[ix]);
      const std::array<std::int64_t,3> brick_number
        {position[ix][0] / bs[0],
         position[ix][1] / bs[1],
         position[ix][2] / bs[2]};

      bool noop = false;
      if (!new_brick)
        noop = true; // Caller didn't actually want to store this.
      else if (inside_user.second[0] == 0)
        noop = true; // Nothing left after clipping.
      else if (old_brick && new_brick->isScalar() && old_brick->isScalar() &&
               new_brick->scalarAsDouble() == old_brick->scalarAsDouble())
        noop = true; // Regardless of any clipping.
      else if (range_inside_brick.first[0] < 0 ||
               range_inside_brick.first[1] < 0 ||
               range_inside_brick.first[2] < 0 ||
               range_inside_brick.second[0] <= 0 ||
               range_inside_brick.second[1] <= 0 ||
               range_inside_brick.second[2] <= 0 ||
               range_inside_brick.second[0] > bs[0] ||
               range_inside_brick.second[1] > bs[1] ||
               range_inside_brick.second[2] > bs[2])
        noop = true; // internal error?

      if (!noop) {
        if (lod == 0 && this->_new_stathist_good) {
          // Subtract old histogram and statistics.
          // Only consider the actual changed region, i.e. inside_user.
          // TODO-WIP-BrickedAPI: Passing "delta" itself might be faster?
          // If the file is opened for update and was missing statistics
          // and histogram then only add, don't subtract old. The result
          // will be bad. And worse with r/m/w, but still somewhat useful
          // for weighted average computation. Do not flush the result.
          // See ZgyInternalBulk::newTrackedBricksTryEnable() for how
          // to test the above.
          SimpleTimerEx t06(_timers.mt6);
          trackChanges(old_brick, range_inside_brick,
                       delta_stats.get(), delta_histo.get(),
                       /*add=*/false, _new_mutex, _timers);
        }

        // A "real" r/m/w, not just caused by edge padding,
        // affects the dirty state (partial or full write).
        const bool rmw = (inside_survey != inside_user);
        rmwlist[ix] = rmw ? 1 : 0;

        SimpleTimerEx t07(_timers.mt7);
        new_brick = doReadModifyWrite
          (old_brick, new_brick, brick_range, inside_survey, inside_user,
           missingValue(/*as_float=*/false), compressor, immutable_buffer);
        t07.stop();

        // Add new histogram and statistics.
        if (lod == 0) {
          SimpleTimerEx t08(_timers.mt8);
          trackChanges(new_brick, range_inside_brick,
                       delta_stats.get(), delta_histo.get(),
                       /*add=*/true, _new_mutex, _timers);
        }

        SimpleTimerEx t09(_timers.mt9);
        std::shared_ptr<const WriteBrickArgPack> args =
          std::make_shared<const WriteBrickArgPack>
          (brick_number, lod, new_brick, compressor, /*offset=*/0);
        // Change scalar to normal or vice versa. Decide allocation.
        args = prepareWriteOneBrick(*args);
        if (args->data->isScalar()) {
          // Queue task to update the lookup table.
          const_queue[ix] = args;
        }
        else {
          // Padding, compression, etc. and then enqueue.
          normal_queue[ix] = prepareWriteOneNormalBrick(*args);
        }
      }
    });

  t03.stop();

  // The code can now assume that the data will either be written
  // or the file will be considered corrupted. The statistics,
  // histogram, and map of dirty bricks can now all be updated.
  // Locking is a bit paranoid because the loop iss done. And
  // higher-level writes are supposed to be serialized.
  {
    SimpleTimerEx t04(_timers.mt4);
    std::lock_guard<std::mutex> lk(_new_mutex);
    t04.stop();
    SimpleTimerEx t05(_timers.mt5);
    std::int64_t neg1 = negativeBins(*_new_modified_histo);
    // Cannot use just a += because target might need to expand.
#if 1
    ExpandableBuilder::staticAddOrSub
      (_new_modified_stats.get(), _new_modified_histo.get(),
       *delta_stats, *delta_histo, /*add=*/true);
#else
    // Logging wrapper
    addChanges
      (_new_modified_stats.get(), _new_modified_histo.get(),
       *delta_stats, *delta_histo, /*add=*/true,
       index3_t{0,0,0}, user_clip.second, true);
#endif
    for (std::size_t ix = 0, end = brickposlist.size(); ix < end; ++ix)
      _new_modified_bricks->set1Written(brickposlist[ix], lod, rmwlist[ix] > 0);
    std::int64_t neg2 = negativeBins(*_new_modified_histo);
    // Not necessarily a bug; it might be numerical instability.
    if (neg1 != neg2)
      blogger(0, std::stringstream().flush()
              << "WARNING: Negative bin count"
              << " in " << neg2 - neg1 << " samples"
              << " (now " << neg2 << ")");
    t05.stop();
  }

  // Note errorhandling:
  // If there are any errors during actual writing this probably
  // means the entire file is a lost cause. This is true also for
  // ZgyUserError raised in the file layer, because at that layer
  // the "user" is really OpenZGY and not some client code. The only
  // acceptable error is ZgySegmentIsClosed, and that will be caught
  // and handled at lower levels.
  // Might implement retrying of writes at a lower level.
  // In that case we still shouldn't see those errors here.
  // There doesn't seem to be anything to gain fron that change.

  // The writeOneNormalBrickWithRetry() method is not threadsafe and *must* be
  // called sequentially. It is also highly recommended to write
  // the data in the same order as it occurred in the work array.
  // Otherwise the bricks get scrambled and reading will be less
  // efficient. This limitation prevents us from having a separate
  // worker process draining the normal_queue as data becomes
  // available. writeOneConstantBrick() doesn't care about order
  // but might as well follow the same rules. That function is
  // very lightweight.

  t1.stop();
  SimpleTimerEx t2(_timers.st);
  ErrorsWillCorruptFile watchdog(this);
  for (const auto& it : const_queue)
    if (it)
      writeOneConstantBrick(*it);
  for (const auto& it : normal_queue)
    if (it)
      writeOneNormalBrickWithRetry(*it);
  t2.stop();

  // Generate some low resolution bricks if possible.
  // Only call genlod when lod==0, because it will process all layers.
  //
  //  The most useful is (1, -1) in tryToCall() and (2, -1) on finalize.
  //  (1, 1) in tryToCall() might get better WeightedAverage results
  //  at the cost of deferring lods 2 and up until the end.
  //  (0, 0) for tryToProcess and (3, -1) on finalize will do a full.
  if (this->_new_genlod && lod == 0) {
    SimpleTimerEx t4(_timers.gl);
    this->_new_genlod->tryToCall
      (_new_internal_lod_mode.incr, brickposlist, nullptr);
  }

  watchdog.disarm();
}

/**
 * \brief Set every single brick in every lod level to the same constant value.
 *
 * \details
 * Set every single brick in every lod level to the same constant
 * value. If force is true, already allocated bricks will be converted
 * back to all-const. THIS WILL LEAK THOSE ALLOCATED BRICKS.
 *
 * This is a higher level version of \c writeAllConstantBrick.
 * The function is at the same level as writeAlignedBrickList.
 * In addition to writing bulk data, the function also resets the statistics,
 * histogram, and dirty block tracking. It also does some more
 * error handling.
 */
void
ZgyInternalBulk::writeAllConstantSurvey(double value)
{
  if (blogger(1, ""))
    blogger(1, "reset entire survey to " + std::to_string(value));
  // TODO-WIP-BrickedAPI: Is it ok to leak allocated bricks here?
  // See more comments in ZgyInternalBulk::newTrackedBricksTryEnable().
  // This might have helped make a better histogram.
  // TODO-WIP-BrickedAPI: Unittest "api.stats_one", possibly others.
  ErrorsWillCorruptFile watchdog(this);
  (void)writeAllConstantBrick(value, /*force=*/true);
  newTrackedBricksSetAllConst(value);
  _new_modified_bricks->setAllClean();
  watchdog.disarm();
}

/**
 * Finish the work done at the end of writeAlignedBrickList() by
 * processing all the partial bricks. The API layer will request this
 * as part of the finalize step. The API needs to call us (the bulk
 * layer) because, in spite of the _new_genloc instance being created
 * by the API layer, only we keep a reference to it.
 */
void
ZgyInternalBulk::newFinalize(const std::function<bool(std::int64_t,std::int64_t)>& p)
{
  if (this->_new_genlod) {
    SimpleTimerEx tt(_timers.gl);
    this->_new_genlod->tryToCall
      (_new_internal_lod_mode.last, std::vector<index3_t>{}, p);
    // Finalize is often called twice. Explicitly by the application,
    // and unconditionally when the file is closed. If the user asked
    // for LodMode::Rebuild it doesn't make sense to do this more than
    // once. If the application writes more data between finalize and
    // close then the "force" parameter is changed to process dirty
    // data only.
    if (_new_internal_lod_mode.last.level < 0 && _new_internal_lod_mode.last.force == 3)
      _new_internal_lod_mode.last.force = 2;
  }
}

/**
 * \brief Split region into bricks.
 * \details
 *
 * The assumption is that the data is to be written.
 *
 * The data is normally copied into newly allocated brick sized
 * buffers, ready to be written from there. If the input is already
 * a single, aligned brick then the input data buffer is returned
 * without copying.
 *
 * TODO-WIP-BrickedAPI: Worry: CAVEAT: Caller must not assume that it
 * owns the buffer. This is unfortunate. Because owning the buffer
 * means that padding etc. can be written directly into the buffer.
 * Testing is also tricky.
 *
 * \callgraph \callergraph
 */
void
ZgyInternalBulk::splitRegionIntoBricks(
     std::vector<std::shared_ptr<DataBuffer>>* bricks /*out*/,
     std::vector<index3_t>* positions /*out*/,
     const std::shared_ptr<DataBuffer>& data,
     const std::array<std::int64_t,3>& start,
     const index3_t& bs)
{
  if (expedited_write() &&
      (start[0] % bs[0]) == 0 &&
      (start[1] % bs[1]) == 0 &&
      (start[2] % bs[2]) == 0 &&
      data->size3d()[0] == bs[0] &&
      data->size3d()[1] == bs[1] &&
      data->size3d()[2] == bs[2])
  {
    bricks->push_back(data);
    positions->push_back(start);
    return;
  }

  const std::pair<index3_t,index3_t> user_range{start, data->size3d()};
  const index3_t realstart // rounded down
    {(start[0] / bs[0]) * bs[0],
     (start[1] / bs[1]) * bs[1],
     (start[2] / bs[2]) * bs[2]};
  const index3_t end // not rounded up
    {start[0] + data->size3d()[0],
     start[1] + data->size3d()[1],
     start[2] + data->size3d()[2]};

  std::shared_ptr<DataBuffer> scalarbrick = !data->isScalar() ? nullptr :
    DataBuffer::makeScalarBuffer3d(data->scalarAsDouble(),bs,data->datatype());

  std::pair<index3_t,index3_t> pos{{0,0,0}, bs};
  for (pos.first[0]=realstart[0];pos.first[0]<end[0];pos.first[0]+=bs[0]) {
    for (pos.first[1]=realstart[1];pos.first[1]<end[1];pos.first[1]+=bs[1]) {
      for (pos.first[2]=realstart[2];pos.first[2]<end[2];pos.first[2]+=bs[2]) {
        if (scalarbrick) {
          bricks->push_back(scalarbrick);
          positions->push_back(pos.first);
        }
        else {
          const bool partial = clip(pos, user_range) != pos;
          std::shared_ptr<DataBuffer> brick =
            DataBuffer::makeNewBuffer3d(bs, data->datatype());
          if (partial)
            brick->fill(0);
          brick->copyFrom(data.get(), start.data(), pos.first.data(), 0, 0);
          bricks->push_back(brick);
          positions->push_back(pos.first);
        }
      }
    }
  }
}

/**
 * \brief Write an arbitrary region.
 *
 * \details \callgraph \callergraph
 *
 * Apply conversion float -> storage and read/modify/write logic. Keep
 * track of min and max sample range for the file as a whole. Pass on
 * the write request to the next level down.
 *
 * The buffer's value type should be either the same as the file's
 * value type (if is_storage) or float (if !is_storage). The is_storage
 * argument is thus redundant since it could have been derived from the
 * value type. But it is better to have the caller be explicit here.
 *
 * The start position refers to the specified lod
 * level. At lod 0 start + data.size can be up to the survey size. At
 * lod 1 the maximum is just half that, rounded up.
 *
 * Note that when writing to the cloud or writing compressed files it
 * is highly recommended to only write regions where both start and
 * size are aligned to the brick size. Or start is aligned and size is
 * up to the end of the survey. This will skip the read/modify/write
 * step. Read/modify/write may cause wasted space if on the cloud and
 * severe compression artifacts if compressed. For the compression case
 * aligned regions might even be enforced. See UpdateMode.
 *
 * \internal
 * Performance note: The read/modify/write could also have been done
 * one brick at a time. Doing it here means that for large requests
 * a number of bricks which were already full will also be read.
 * On the other hand a single read might help parallelism.
 *
 * The issue gets more important with support for incremental finalize,
 * because this triggers the r/m/w logic on every write, not just those
 * that are misaligned. This is to keep the statistics and histogram
 * information correct.
 *
 * TODO-@@@: A writeconst() on a region less than the entire survey but
 * still larger than available memory might fail if the region is
 * misaligned *or* if the file is open for update.
 *
 * TODO-@@@: Implement a compromise where the write request is split
 * into brick columns before the r/m/w logic and into bricks after.
 * Note:
 *
 *  - Do not split a write() or writeconst() covering the entire survey.
 *    That implies that the file is being overwritten, not updated.
 *    The lower level needs to know.
 *
 *  - It is not critical to do this split for a write() because the caller
 *    must then already have allocated a buffer for the whole region,
 *    So we can probably allocate another one. There might still be a
 *    performance benefit (less data read) or not (less parallelism
 *    and brick consolidation) and what is best depends on whether
 *    we are tracking changes (so every brick must be read anyway) and
 *    whether the file is on the cloud and how its bricks are sorted.
 *
 * TODO-@@@: If the file is open for update and if old contents are
 * only needed to update stats not for r/m/w then try a readconst
 * first. Especially if the new data to be written is a large
 * writeconst. This can save inflating some scalar buffers. Which can
 * also reduce the risk of running out of memory in a few special
 * cases with a large writeconst that is still smaller than the entire
 * survey. Unfortunately this change complicates the r/m/w logic even
 * further. It adds a bad code smell.
 */
void
ZgyInternalBulk::writeRegion(
      const std::shared_ptr<DataBuffer>& data_in,
      const std::array<std::int64_t,3>& start_in,
      int32_t lod, bool is_storage,
      const compressor_t& compressor)
{
  validateUserLod(lod);
  validateUserPosition(start_in, data_in->size3d(), lod);
  std::shared_ptr<DataBuffer> data = data_in;
  std::array<std::int64_t,3> start = start_in;

  if (data->datatype() != (!is_storage ? RawDataType::Float32 :
                           _metadata->ih().datatype()))
    throw OpenZGY::Errors::ZgyUserError("Invalid data type in writeRegion");

  if (blogger(2))
    blogger(2, std::stringstream()
            << "write(start="
            << "(" << start[0]
            << "," << start[1]
            << "," << start[2]
            << "), size="
            << "(" << data->size3d()[0]
            << "," << data->size3d()[1]
            << "," << data->size3d()[2]
            << "), lod=" << lod
            << std::boolalpha << ", is_storage=" << is_storage << ")");

  // TODO-Performance: Combining range() and scaleDataToStorage()
  // might save some time.

  if (!is_storage) {
    SimpleTimerEx t1(_timers.sts);
    data = scaleDataToStorage(data); // @@ "data" can also be changed if r/m/w
  }

  // The application can call writeconst() immediately after creating a file.
  // This can be useful because it establishes the value for "empty" samples.
  // This can also prevent getting some arbitrary (usually zero) value added
  // to the statistical and histogram range even if it doesn't exist in the
  // samples.
  //
  // The code must treat this as a special case. Among other things,
  // the lowres can be trivially computed and marked as ready.
  //
  // It is currently unspecified what happens in the following cases:
  //   - Reset the survey after real data has been written.
  //   - write() with a huge data that happens to be all constant.
  //   - writebricks() listing all bricks in the survey.
  //
  // It is debatable whether recognizing any of the above is useful.
  // Currently a reset is recognized after writing real data, but
  // space is leaked. writebricks() with a list of same value scalar
  // bricks currently works (currently checked for in
  // writeAlignedBrickList) but is hardly something the user would
  // need to do. It *might* be useful if the histogram range needs to
  // be reset. But deleting and re-creating the file should work just
  // as well.
  const bool entire_survey =
    covers(start, data->size3d(), std::array<std::int64_t,3>{0,0,0}, _metadata->ih().size());
  if (entire_survey && data->isScalar()) {
    writeAllConstantSurvey(data->scalarAsDouble());
    return;
  }

  const std::pair<index3_t,index3_t> user_range{start_in, data->size3d()};
  std::vector<std::shared_ptr<DataBuffer>> bricks;
  std::vector<index3_t> positions;
  splitRegionIntoBricks(&bricks, &positions, data, start_in,
                        this->_metadata->ih().bricksize());
  writeAlignedBrickList(bricks, positions, user_range,
                        lod, is_storage, compressor,
                        /*immutable_buffer=*/true);
}

/**
 * Make sure the caller provided a single aligned brick and optionally
 * a matching buffer.
 *
 * If strict=true, all errors will throw. If false, a few errors will
 * quietly return false, instead.
 *
 * Positions outside the survey are not an error, but misaligned bricks are.
 *
 * See also validateUserPosition() and validateUserLod().
 */
bool
ZgyInternalBulk::checkBrickInternal(
    const std::array<int64_t,3>& position,
    const std::shared_ptr<DataBuffer>& data,
    bool as_float, bool strict) const
{
  // Buffer is present.
  if (!data || !data->voidData() || data->isScalar()) {
    if (strict)
      throw OpenZGY::Errors::ZgyUserError
        ("Read or write to missing or scalar buffer is not allowed.");
    else
      return false;
  }
  // Buffer is of the correct value type.
  if (as_float) {
    if (data->datatype() != RawDataType::Float32)
      throw OpenZGY::Errors::ZgyInternalError
        ("Converted samples can only be stored in float buffers.");
  }
  else {
    if (data->datatype() != this->_metadata->ih().datatype())
      throw OpenZGY::Errors::ZgyInternalError
        ("Buffer type does not match storage value type.");
  }
  // Buffer is of the correct size (a single brick).
  const std::array<std::int64_t,3> bs = this->_metadata->ih().bricksize();
  if (data->size3d() != bs)
      throw OpenZGY::Errors::ZgyUserError
        ("The brick API expects brick-sized data buffers.");
  // Buffer position aligned to bricks.
  for (int dim=0; dim<3; ++dim)
    if ((position[dim] % bs[dim]) != 0)
      throw OpenZGY::Errors::ZgyUserError
        ("Brick start must be aligned with brick size.");
  return true;
}

/**
 * Make sure the caller provided the correct number and type
 * of brick-sized data buffers.
 *
 * See also validateUserPosition() and validateUserLod().
 */
void
ZgyInternalBulk::checkBricksInternal(
    const std::vector<std::array<int64_t,3>>& position,
    const std::vector<std::shared_ptr<DataBuffer>>& data,
    bool as_float, bool strict) const
{
  if (position.size() != data.size())
    throw OpenZGY::Errors::ZgyUserError
      ("The brick API expects the same number of positions as buffers");
  for (std::size_t ii = 0; ii < data.size(); ++ii)
    (void)checkBrickInternal(position[ii], data[ii], as_float, strict);
}

/////////////////////////////////////////////////////////////////////////////
// bulk.h --- BRICK API --- public parts ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \brief Internal reads from genlod end up here.
 *
 * \details \callgraph \callergraph
 *
 * As readBricksInternal(), but passing buffers to store
 * the result is now optional. And This function is not
 * required to use them. The result can be a mix of
 * regular and scalar buffers. No result will be null.
 *
 * check_constant=true means the code will also return
 * a scalar buffer for a brick not explicitly flagged
 * as const, but just contains a single sample value.
 * Caveat: this test will also look at padding samples.
 */
std::vector<std::shared_ptr<DataBuffer>>
ZgyInternalBulk::readBricksToNewBuffers(
     const std::vector<std::array<std::int64_t,3>>& position_in,
     const std::vector<std::shared_ptr<DataBuffer>>& data_in,
     int lod, bool as_float, bool check_constant) const
{
  std::vector<std::array<std::int64_t, 3>> position(position_in);
  std::vector<std::shared_ptr<DataBuffer>> data(data_in);
  data.resize(position_in.size());

  // Invariants for this file.
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const std::array<std::int64_t,3> bricksize =
    ih.bricksize();
  const RawDataType dtype =
    (as_float ? RawDataType::Float32 : ih.datatype());
  const std::array<std::int64_t, 3> ssize =
    ih.lodsizes()[lod] * bricksize;

  std::vector<std::pair<bool, double>> cvbricks =
    readConstValueBricks(position, lod, as_float);

  // For any brick flagged as all-const, ignore the provided
  // buffer (if any) and return a scalar buffer instead.
  // In the same loop, allocate any missing data buffers.
  for (std::size_t ii = 0; ii < position.size(); ++ii)
  {
    const bool outside = (position[ii][0] >= ssize[0] ||
      position[ii][1] >= ssize[1] ||
      position[ii][2] >= ssize[2]);
    if (outside) {
      data[ii] = DataBuffer::makeScalarBuffer3d
        (missingValue(as_float), bricksize, dtype);
    }
    else if (cvbricks[ii].first) {
      data[ii] = DataBuffer::makeScalarBuffer3d
        (cvbricks[ii].second, bricksize, dtype);
      position[ii] = std::array<std::int64_t, 3>{ -bricksize[0], 0, 0 };
    }
    else if (data[ii] == nullptr || data[ii]->isScalar()) {
      data[ii] = DataBuffer::makeNewBuffer3d(bricksize, dtype);
    }
  }

  // This will ignore any entries with null or scalar buffers.
  // Also outside-survey since we already handled those.
  readBricksInternal(position, data, lod, as_float);

  // Caller wants a thorough check for all const.
  if (check_constant)
    for (std::size_t ii = 0; ii < position.size(); ++ii)
      if (data[ii] != nullptr && !data[ii]->isScalar())
        if (data[ii]->isAllSame(data[ii]->size3d().data()))
          data[ii] = DataBuffer::makeScalarBuffer3d
            (data[ii]->scalarAsDouble(),
             data[ii]->size3d(),
             data[ii]->datatype());

  return data;
}

/**
 * \brief Calls to IZgyReader::readbricks() end up here.
 *
 * \details \callgraph \callergraph
 *
 * The caller needs to convert the data buffers between simple smart
 * pointers and our internal DataBuffer type.
 *
 * The caller is responsible for providing the correct number and type
 * of brick-sized data buffers.
 */
void
ZgyInternalBulk::readBricksInternal(
    const std::vector<std::array<int64_t,3>>& position,
    const std::vector<std::shared_ptr<DataBuffer>>& data_in,
    int lod, bool as_float) const
{
  checkBricksInternal(position, data_in, as_float, /*strict=*/false);

  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const std::array<std::int64_t,3> bs  = ih.bricksize();
  const std::array<std::int64_t,3> ssize = ih.lodsizes()[lod] * bs;

  // Eligible for a shortcut if start & size is ok.
  // For cloud access it is probably more important to consolidate
  // neighboring reads than to maybe elide a memcopy. So, only try
  // that if the request was for just a single brick.
  const bool maybe_expedited = position.size()==1 || !this->_file->xx_iscloud();

  std::vector<std::shared_ptr<DataBuffer>> data(data_in);
  for (std::size_t ii=0; ii<position.size(); ++ii) {
    const bool outside = (position[ii][0] < 0 || position[ii][0] >= ssize[0] ||
                          position[ii][0] < 0 || position[ii][1] >= ssize[1] ||
                          position[ii][0] < 0 || position[ii][2] >= ssize[2]);
    if (!data[ii] || !data[ii]->voidData() || data[ii]->isScalar())
    {
      if (blogger(3, ""))
        blogger(3, std::stringstream()
          << "    No buffer   " << lod << " pos " << fmt(position[ii]));
    }
    else if (!outside) {
      if (maybe_expedited && expeditedRead
          (position[ii],
           data[ii]->size3d(),
           data[ii]->voidData().get(),
           lod,
           data[ii]->datatype()))
      {
        if (blogger(3, ""))
          blogger(3, std::stringstream()
            << "    Expedited read lod " << lod << " at " << fmt(position[ii]));
        data[ii] = nullptr; // mark as alrady done.
      }
      else
      {
        // Reading now done at the end of the function.
        //this->readToExistingBuffer(data[ii], position[ii], lod, as_float);
        if (blogger(3, ""))
          blogger(3, std::stringstream()
            << "    Read later lod " << lod << " at " << fmt(position[ii])) ;
      }
    }
    else {
      if (blogger(3, ""))
        blogger(3, std::stringstream()
                << "    Not reading " << lod << " pos " << fmt(position[ii]));
      // TODO-WIP-BrickedAPI: If there are all-const buffers among
      // the inputs then perhaps choose one of those. But be careful
      // to keep the behavior consistent. The idea is that if the
      // user initialized the entire survey with a const value
      // then this should also be used for padding, so the lowres
      // bricks can still be stored as all-const instead of being
      // a mix of user's and system't default.
      // Note that in readBricksToNewBuffers() this would (should?)
      // cause the corresponding output buffer to be scalar.
      data[ii]->fill(missingValue(as_float));
      data[ii] = nullptr; // mark as alrady done.
    }
  }
  // Do all remaining real reads in one call.
  readBricksToExistingBuffer(data, position, lod, as_float);
}

/**
 * \brief Calls to IZgyWriter::writebricks() end up here.
 *
 * \details \callgraph \callergraph
 *
 * The caller needs to convert the data buffers between simple
 * smart pointers and our internal DataBuffer type.
 *
 * \param position           Position of each brick, in traces and samples
 *                           at this lod. Not the brick number.
 *
 * \param data               The bulk data to be written, as a list of
 *                           DataBuffer instances. Which is a type not
 *                           exposed in the public API. Both regular and
 *                           scalar (all-const) buffers are supported.
 *                           However, if the request came from writeconst()
 *                           then it is slightly more performant to call
 *                           (TODO-WIP-BrickedAPI: not yet implemented)
 *                           writeConstBricksInternal().
 *
 * This internal version has parameters not exposed in the
 * planned writebricks() in the public API.
 *
 * \param lod                If not zero, low resolution data is being
 *                           written. Only genlod is allowed to do that.
 *
 * \param is_storage         If false, the caller passed a float buffer
 *                           regardless of what the storage type is.
 *                           Which means the samples may need to be scaled,
 *                           This can often be derived from the value type
 *                           of the buffer.
 *
 * \param compressor         The compressor functor for this lod, specified
 *                           when the file was opened.
 *
 * \param immutable_buffer   If false, function is allowed to modify the
 *                           input data buffers in-place. In most cases the
 *                           caller won't have any problems with that.
 *                           E.g. genlod will always pass false.
 *                           If true, an extra buffer copy may be needed.
 *                           In some cases a buffer copy is done just
 *                           because it *might* be needed.
 *
 * \param maybe_more         Assume this brick will be written many more
 *                           times. Typically this means that the general
 *                           API was used, and read/modify/write was needed.
 *                           Generating low resolution bricks for this
 *                           data needs to wait until the file is closed.
 *                           Note, the code might be able to figure this
 *                           outself when it applies "user" clipping.
 */
void
ZgyInternalBulk::writeBricksInternal(
     const std::vector<std::array<int64_t,3>>& position,
     const std::vector<std::shared_ptr<DataBuffer>>& data,
     int lod, bool is_storage, const compressor_t& compressor,
     bool immutable_buffer, bool maybe_more)
{
  checkBricksInternal(position, data, !is_storage, /*strict=*/false);
  writeAlignedBrickList
    (data, position,
     std::pair<index3_t,index3_t>{{0,0,0},{0,0,0}},
     lod, is_storage, compressor, immutable_buffer);
}

/**
 * \brief Get the constant sample value of each provided brick.
 *
 * \details \callgraph \callergraph
 *
 * TODO-WIP-BrickedAPI: Performance issue.
 * Instead of just looping and calling the old code,
 * implement a new readConstantValue() and partsNeeded()
 * that work on list of bricks.
 *
 * TODO-WIP-BrickedAPI: Expose in the public api
 * as readconstbricks().
 */
std::vector<std::pair<bool, double>>
ZgyInternalBulk::readConstValueBricks(
  const std::vector<std::array<std::int64_t, 3>>& position,
  int32_t lod, bool as_float) const
{
  const std::array<std::int64_t, 3> bs = this->_metadata->ih().bricksize();
  std::vector<std::pair<bool, double>> result;
  for (const auto& it : position)
    if (singleBrickOutsideSurvey(it, bs, lod))
      result.push_back(std::make_pair(true, missingValue(as_float)));
    else
      result.push_back(readConstantValue(it, bs, lod, as_float));
  return result;
}

/**
 * \brief This is currently a stub.
 *
 * \details \callgraph \callergraph
 *
 * Equivalent to calling writeBricksInternal() to set every sample
 * in the region to the same value. Called from writeconst in the
 * general API, and maybe internally. Using writeBricksInternal()
 * might lead to running out of memory. As well as being slower.
 *
 * It is particularly important to use this function instad of
 * writeBricksInternal() when the entire survey is being set to
 * the same value. Because in that case the lod and statistics
 * tracking and min/max tracking should all be reset.
 *
 * Value should always be a scalar buffer. The scalar value,
 * valuetype, and size are extracted from the buffer.
 *
 * Implementation note: Except for the most simple cases the code
 * may end up falling back to writeBricksInternal() anyway.
 * E.g. for already allocated bricks or r/m/w processing.
 * But it can do things such as splitting the request into
 * brick-columns if it looks like we might run out of memory.
 */
void
ZgyInternalBulk::writeConstBricksInternal(
     const std::vector<std::array<int64_t,3>>& position,
     const std::shared_ptr<InternalZGY::DataBuffer>& value,
     bool is_storage)
{
  throw std::runtime_error("writeConstBricksInternal not implemented yet");
}

} // namespace
