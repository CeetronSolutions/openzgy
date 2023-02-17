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
 *     easier to implement (and obfuscate the code a bit), _writeOneBrick()
 *     and _writeOneNormalBrick() instead return to the calling
 *     _writeAlignedRegion() with an argument package that should be sent
 *     to the next method in the chain.
 *
 *     The order of calls is:
 *
 *     writeRegion (arbitrary region, user's value_type)
 *
 *         Apply conversion float -> storage and read/modify/write logic. Keep
 *         track of min and max sample range for the file as a whole.
 *
 *     _writeAlignedRegion (arbitrary aligned region, native vt)
 *
 *         Apply splitting into bricks.
 *         Expand bricks at survey edge to full size unless r/m/w already did.
 *         Convert sample position to brick number. I.e. divide by bricksize.
 *
 *     _writeOneBrick (single brick, native vt)
 *
 *         Might convert from a scalar to a regular buffer or vice versa.
 *         Might decide whether to update an existing brick, and where.
 *         Might decide to veto the compression.
 *         The next step depends on whether the buffer ended up scalar
 *         (_writeOneConstantBrick) or not (_writeOneNormalBrick eventually
 *         followed by _writeWithRetry)
 *
 *     _writeOneNormalBrick
 *
 *         Apply padding, compression, byte swapping, and convert from
 *         DataBuffer to void* + size suitable for the file i/o layer.
 *
 *     _writeOneConstantBrick
 *         Apply conversion of the constant value to a lookup table entry,
 *         Pass on the write request to the function that updates the lut.
 *
 *     _writeWithRetry
 *
 *        Allocate space on the bulk file as needed.
 *        Handle retries needed if the segment (on the cloud) was closed.
 *        Update the lookup table if a new block was allocated.
 *
 *     Code flow for reading
 *
 *     readToNewBuffer
 *     readToExistingBuffer | readConstantValue
 *     _partsNeeded
 *     _getBrickFilePosition
 *     _deliverOneBrick
 *     _file.xx_readv
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
#include "mtguard.h"
#include "histogrambuilder.h"
#include "statisticdata.h"
#include "histogramdata.h"
#include "../exception.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <sstream>
#include <omp.h>

namespace InternalZGY {
#if 0
}
#endif

using namespace InternalZGY::ArrayOps;
using namespace InternalZGY::Formatters;

namespace {
  /**
   * For testing only; might be removed. Apps should have no reason to
   * reset this variable because it is unlikely that they have another
   * multi-threaded loop going at the same time.
   */
  static bool enable_compress_mt()
  {
    static int enable = Environment::getNumericEnv("OPENZGY_ENABLE_COMPRESS_MT", 1);
    return enable > 0;
  }

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
   * \brief Add or subtract this buffer's samples in the statistics.
   *
   * \details
   *
   * This is used to keep track of statistics while data is being updated.
   * For the initial write this will not and cannot be used because
   * the width of the histogram is not known yet.
   *
   * Unlike the old accessor the histogram range is not allowed to change.
   * And the statistical min/max range can only grow. So if some spikes
   * are overwritten their values are still visible in both ranges.
   *
   * The old code would try to keep the statistical range somewhat
   * correct by trimming it to the first and last histogram bin with
   * samples. See HistogramBuildet::operator+=() and StatisticData::
   * trimRange(). The problem is that since the histogram cannot grow,
   * that trick might not work. For consistency don't even try. Just
   * document that the range will never shrink as long as the derived
   * information is generated incrementally.
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
   * because of an obscure issue with count == 0 but in that case
   * is the min/max ramge valid and should the final ran
   */
  template <typename T>
  void
  trackChangesT(const std::shared_ptr<const DataBuffer>& data_in,
                const std::array<std::int64_t,3>& valid_size,
                StatisticData *stats, HistogramData *histo,
                bool add)
  {
    if (!data_in || !stats || !histo)
      return;                   // Not tracking changes.
    std::shared_ptr<const DataBufferNd<T,3>> data =
      std::dynamic_pointer_cast<const DataBufferNd<T,3>>(data_in);
    if (!data)
      return;                   // Maybe throw in this case?
    HistogramBuilder hb(histo->getsize(), histo->getmin(), histo->getmax());
    std::int64_t len = valid_size[0] * valid_size[1] * valid_size[2];
    if (data->isScalar()) {
      hb.add(data->data(), data->data() + 1);
      hb *= len;
    }
    else if (data->size3d() == valid_size) {
      hb.add(data->data(), data->data() + len);
    }
    else {
      const std::array<std::int64_t,3> size = data->size3d();
      for (int ii=0; ii<valid_size[0]; ++ii) {
        for (int jj=0; jj<valid_size[1]; ++jj) {
          const T* ptr = data->data() + ii*size[1]*size[2] + jj*size[2];
          hb.add(ptr, ptr + valid_size[2]);
        }
      }
    }
    if (add) {
      *stats += hb.getstats();
      *histo += hb.gethisto();
    }
    else {
      *stats -= hb.getstats();
      *histo -= hb.gethisto();
    }
    //std::cout << "track changes " << (add?"add":"sub")
    //          << " size " << data->size3d() << " valid " << valid_size
    //          << " stats " << hb.getstats().toString()
    //          << " histo " << hb.gethisto().toString(true)
    //          << "\ncurrent value now"
    //          << " stats " << stats->toString()
    //          << " histo " << histo->toString(true)
    //          << std::endl;
  }

  /**
   * \copydoc TrackChangesT
   * This method has the usual boilerplare for not-quite templated code.
   */
  void
  trackChanges(const std::shared_ptr<const DataBuffer>& data_in,
               const std::array<std::int64_t,3>& valid_size,
               StatisticData *stats, HistogramData *histo,
               bool add)
  {
    switch (data_in->datatype()) {
    case RawDataType::SignedInt8:
      trackChangesT<std::int8_t>(data_in, valid_size, stats, histo, add);
      break;
    case RawDataType::SignedInt16:
      trackChangesT<std::int16_t>(data_in, valid_size, stats, histo, add);
      break;
    case RawDataType::Float32:
      trackChangesT<float>(data_in, valid_size, stats, histo, add);
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
 * Argument package used by _writeOneBrick(), _writeOneNormalBrick(),
 * _writeOneConstantBrick(). One size fits all. Which is in general
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
 *   * fileoffset:    Only when known. Never known yet in _writeOneBrick(),
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
 * Argument package used by _writeWithRetry() which is at a lower level
 * than _writeOneBrick() etc. and it got too awkward to use the
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
 * instead call a _get_metadata_rw() method. If ZgyInternalBulk is
 * const-correct then you get a compile time error if trying to call
 * _get_metadata_rw() from inside a method declared as const. If not,
 * _get_metadata_rw() is still preferable because it can raise a proper
 * error message instead of the null pointer exception.
 */
ZgyInternalBulk::ZgyInternalBulk(
    const std::shared_ptr<IFileADT>& file,
    const std::shared_ptr<const ZgyInternalMeta>& metadata,
    const std::shared_ptr<ZgyInternalMeta>& metadata_rw,
    bool compressed_write,
    const LoggerFn& logger)
  : _file(file)
  , _metadata(metadata)
  , _metadata_rw(metadata_rw)
  , _update_mode(compressed_write? UpdateMode::Constant : UpdateMode::Always)
  , _compressed_write(compressed_write)
  , _is_bad(false)
  , _written_sample_min(std::numeric_limits<double>::infinity())
  , _written_sample_max(-std::numeric_limits<double>::infinity())
  , _loggerfn(logger ? logger : LoggerBase::standardCallback(LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"), "openzgy-bulk: ", ""))
  , _ptimer_st(new SummaryPrintingTimerEx("writeAligned[S]"))
  , _ptimer_mt(new SummaryPrintingTimerEx("writeAligned[M]"))
  , _ptimer(new SummaryPrintingTimerEx("writeAligned"))
  , _ststimer(new SummaryPrintingTimerEx("scaleToStorage"))
  , _modified_bricks()
  , _modified_stats()
  , _modified_histo()
{
  // If the file is being opened for update there may already be data in it.
  // The histogram will show the range of written samples to date, even
  // if the histogram itself (and statistics etc.) have not been generated.
  // And if it isn't generated it will need to decide on a value range later.
  // And if so it needs to include the range of data written earlier, not
  // just what is written in this session. So the written range will be
  // initialized from the histogram.
  //
  //  - The range will never shrink, even if the application overwrites
  //    all the blocks with spikes in them. This limitation applies also
  //    when all writes are in the same session. Could be fixed but I
  //    doubt anybody cares much.
  //
  //  - Only handling float cubes here. Currently the histogram has a
  //    fixed range when the data is integral. So it doesn't need to
  //    keep track of the written range. If that changes then the code
  //    here needs to be updated. _written_sample_min/max is in
  //    storage values while the histogram is in real values so there
  //    needs to be an explicit conversion.
  //
  //  - If the histogram already contains samples then this suggests
  //    we will be doing a partial finalize. In that case the histogram
  //    limits cannot currently be changed. _written_sample_min/max
  //    still gets kept up to date but it wouldn't matter if it didn't.
  //
  //  - So, the situation where this code is needed is when the file
  //    is written more than once and the caller decided not to make
  //    use of incremental finalize.
  //
  //  - If the application wrote some zero samples to the file and
  //    nothing else, then re-opened the file to write all positive or
  //    all negative values, then zero should obviously be included in
  //    the histogram range but it might not be. This is technically a
  //    bug but it is really unlikely anybody will notice. In fact,
  //    other code might decide to include zero anyway because empty
  //    bricks exist. The reason for the problem in this function is
  //    that it is not possible to distinguish between a never written
  //    histogram and one that tells us that zeros were written prior
  //    to this histogram being cleared.
  if (metadata_rw) {
    const IHistHeaderAccess& hh = metadata_rw->hh();
    const IInfoHeaderAccess& ih = metadata_rw->ih();
    if (hh.minvalue() <= hh.maxvalue()) {
      if (hh.samplecount() != 0 || hh.minvalue() != 0 || hh.maxvalue() != 0) {
        if (ih.datatype() == RawDataType::Float32) {
          if (_logger(1))
            _logger(1, std::stringstream()
                    << "Include previously seen data range: "
                    << hh.minvalue() << " .. " << hh.maxvalue());
          _written_sample_min = hh.minvalue();
          _written_sample_max = hh.maxvalue();
        }
      }
    }
  }

  // Valid statistics, histogram, and lowres means we might want to do
  // an incremental finalize. So keep track of changes. The method
  // will only turn it on if valid data exists.
  this->trackedBricksTryEnable(true);
}

/**
 * \brief Get hint about all constant region.
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
  _validateUserPosition(start, size, lod);
  const double defaultstorage = this->_metadata->ih().defaultstorage();
  std::vector<LutInfoEx> bricklist = _partsNeeded(start, size, lod);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  double result = nan;
  bool first = true;
  for (const LutInfoEx& brick : bricklist) {
    if (_logger(5))
      _logger(5, std::stringstream()
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
      if (!first && result != defaultstorage)
        return std::make_pair(false,nan);
      result = defaultstorage;
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
  if (_logger(2))
    _logger(2, std::stringstream()
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
  _validateUserPosition(start, result->size3d(), lod);
  RawDataType expected_dtype = as_float ? RawDataType::Float32 : _metadata->ih().datatype();
  // When called from ZgyReader.read(), the check on expected_dtype represents
  // a user error trying to convert to an intergral type while the others
  // should not be possible to trigger.
  if (!result || result->isScalar() || !result->voidData())
    throw OpenZGY::Errors::ZgyInternalError("Reading to missing or wrong buffer type.");
  if (result->datatype() != expected_dtype)
    throw OpenZGY::Errors::ZgyUserError("Requested data type not supported for this file.");

  if (_logger(2))
    _logger(2, std::stringstream()
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

  // Need a default value to use when trying to read a brick that
  // was never written, or to fill in a brick that was only partly
  // written. To avoid non intuitive behavior the same value should
  // be used for both cases, and it should be a value that would
  // also be allowed to store as a regular sample. Use the value
  // that becomes 0.0 after conversion to float. Or as close as
  // possible to that value if the coding range is not zero centric.
  // Always use zero for floating point data. (Not NaN...)
  // Dead traces don't have any special handling apart from the
  // alpha flag. They still contain whatever data was written to them.

  const double defaultstorage = this->_metadata->ih().defaultstorage();
  const double defaultvalue = this->_metadata->ih().defaultvalue();

  // Make a separate pass to gather all the bricks we need to read.
  // The lower levels might fetch some of them in parallel and we might be
  // able to combine bricks to read larger blocks at a time. Both changes
  // can have a dramatic performance impact effect on cloud access.

  std::vector<LutInfoEx> bricklist = _partsNeeded(start, result->size3d(), lod);

  // After all bricks have been processed, the padding past the
  // end if the survey might still not have been touched. Just in
  // case the request did in fact include such samples we will
  // initialize the entire result buffer to the default value.
  // TODO-Performance can I refactor the code so this isn't needed?
  // Or can I make it fill only the actual padding area?
  result->fill(as_float ? defaultvalue : defaultstorage);

  std::vector<ReadRequest> requests;
  for (const LutInfoEx& brick : bricklist) {
    // Functor that accepts a raw void* data + size and copies it into the
    // correct place in result. result, start, and as_float come straight
    // from user code. b.position and b.status vary per brick which is why
    // I copy them to make sure the up to date contents are captured. raw and
    // rawsize (the arguments to the lambda) is what was read from storage.

    const std::array<std::int64_t,3> bpos
      {brick.survey_position.i0,
       brick.survey_position.j0,
       brick.survey_position.k0};
    const BrickStatus bstatus = brick.status;
    auto deliverance = [this,result,start,as_float,bpos,bstatus](ReadRequest::data_t raw, std::int64_t rawsize) {
                         this->_deliverOneBrick
                           (result, start, bpos, raw, rawsize,
                            bstatus, as_float);
                       };
    switch (brick.status) {
    case BrickStatus::Missing:
      if (_logger(2))
        _logger(2, std::stringstream()
                << "  Reading brick at " << brick.survey_position
                << " not found, use " << defaultstorage << "\n");
      deliverance(std::make_shared<double>(defaultstorage), sizeof(double));
      break;

    case BrickStatus::Constant:
      if (_logger(2))
        _logger(2, std::stringstream()
                << "  Reading constant brick at "
                << brick.survey_position << "\n");
      deliverance(std::make_shared<double>(brick.double_constvalue), sizeof(double));
      break;

    case BrickStatus::Normal:
      if (_logger(2))
        _logger(2, std::stringstream()
                << "  Reading brick at " << brick.survey_position
                << " from file offset " << std::hex
                << (std::intptr_t)brick.offset_in_file << std::dec << "\n");
      requests.push_back(ReadRequest{brick.offset_in_file, this->_metadata->ih().bytesperbrick(), deliverance});
      break;

    case BrickStatus::Compressed:
      if (_logger(2))
        _logger(2, std::stringstream()
                << "  Reading compressed brick at " << brick.survey_position
                << " from file offset " << std::hex
                << (std::intptr_t)brick.offset_in_file << std::dec
                << " size " << brick.size_in_file << "\n");
      // TODO-Worry obscure corner case, might need to re-try if we didn't
      // get enough data. I try to prevent this with the rule that
      // no compressed brick may be larger than the uncompressed version.
      requests.push_back(ReadRequest{brick.offset_in_file, brick.size_in_file, deliverance});
      break;

    default:
      throw OpenZGY::Errors::ZgyInternalError("Internal error, bad brick status");
    }
  }
    if (_logger(2))
      _logger(2, std::stringstream()
              << requests.size() << " read requests are queued\n");
    if (!requests.empty())
      this->_file->xx_readv(requests, true, false, true, UsageHint::Data);

  // Note-Performance: If passing true in the second arguent above this
  // could help performance a lot. Especially for reading compressed files
  // where the user sends large requests without multithreading. Also
  // when finalizing compressed files. parallel_ok=true will cause the
  // decompression step to be multi-threaded. Also the conversion to
  // float (if needed) and the copy-out to the applicaton's buffer will
  // be multi-threaded. But there are caveats:
  //
  // * _deliverOneBrick() must be thread safe.
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
  // * See commants in LocalFileLinux::xx_readv() and _deliverOneBrick().
  //   And GenLodImpl::_calculate().
}

/**
 * Read bulk data starting at "start" in index space size "size".
 * Return the result in a newly allocated DataBuffer. The result will
 * either be a scalar (constant-value) buffer or a regular buffer.
 *
 * Pass check_constant=true to check extra hard for all-constant data.
 * A region written with all samples identical, as opposed to a region
 * flagged as constant without taking up space, will also be detected.
 *
 * Start should be in the range (0,0,0) to Size-1. It is valid to pass
 * a size that includes the padding area between the survey and the
 * end of the current brick. But not past that point.
 */
std::shared_ptr<DataBuffer>
ZgyInternalBulk::readToNewBuffer(
    const std::array<std::int64_t,3>& start,
    const std::array<std::int64_t,3>& size,
    int32_t lod, bool as_float, bool check_constant) const
{
  _validateUserPosition(start, size, lod);
  std::shared_ptr<DataBuffer> result;
  RawDataType dtype = (as_float ?
                       RawDataType::Float32 :
                       this->_metadata->ih().datatype());
  std::pair<bool,double> cvalue = this->readConstantValue(start, size, lod, as_float);
  if (cvalue.first) {
    result = DataBuffer::makeScalarBuffer3d(cvalue.second, size, dtype);
  }
  else {
    result = DataBuffer::makeNewBuffer3d(size, dtype);
    this->readToExistingBuffer(result, start, lod, as_float);
    if (check_constant && result->isAllSame(result->size3d().data())) {
      double scalar = result->scalarAsDouble(); // returns the first value
      result = DataBuffer::makeScalarBuffer3d(scalar, size, dtype);
    }
  }
  return result;
}

/**
 * Start or stop keeping track of changed bricks to support
 * incremental finalize. Any existing information is discarded.
 * Tracking requires pre-existing statistics, histogram, and low
 * resolution bricks. If missing then tracking is turned off,
 * regardless of what the caller asked for.
 *
 * This method must be nonvirtual because it is also called from the
 * constructor.
 */
void
ZgyInternalBulk::trackedBricksTryEnable(bool on)
{
  _modified_bricks.clear();
  _modified_stats.reset();
  _modified_histo.reset();
  if (!_metadata_rw) {
    _logger(1, "Will not track changes. Not open for write.");
  }
  else if (!on) {
    _logger(1, "Will not track changes. Explicitly disabled.");
  }
  else if (!_metadata_rw->can_finalize_incremental()) {
    _logger(1, "Will not track changes. Compressed or missing data.");
  }
  else if (readConstantValue(std::array<std::int64_t,3>{0,0,0},
                             _metadata_rw->ih().size(),
                             /*lod*/0, /*as_float*/false).first)
  {
    _logger(1, "Will not track changes. File is empty.");
  }
  else {
    const IHistHeaderAccess& hh = _metadata_rw->hh();
    const IInfoHeaderAccess& ih = _metadata_rw->ih();
    _modified_bricks.resize(ih.brickoffsets().back(), 0);
    _modified_stats.reset(new StatisticData
                          (ih.scnt(), /*inf=*/0, ih.ssum(), ih.sssq(),
                           ih.smin(), ih.smax()));
    _modified_histo.reset(new HistogramData
                          (hh.bins(), (int)hh.bincount(),
                           hh.minvalue(), hh.maxvalue()));
    if (_logger(2))
      _logger(2, std::stringstream()
              << "Will Track changes. Initial stats (real)    "
              << _modified_stats->toString()
              << " histo " << _modified_histo->toString());
    // Statistics saved in the ZGY file uses float data, but internal
    // calculation needs to be in storage because otherwise any clipped
    // values would mess things up.
    const std::array<double,2> factors = ih.storagetofloat();
    _modified_stats->scale(factors[1], factors[0] + factors[1], 0, 1);
    _modified_histo->scale(factors[1], factors[0] + factors[1], 0, 1);
    if (_logger(1))
      _logger(1, std::stringstream()
              << "Will Track changes. Initial stats (storage) "
              << _modified_stats->toString()
              << " histo " << _modified_histo->toString());
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
  _validateUserPosition(start, size, lod); // Could be simplified.

  // TODO-Performance: Could be simplified.
  std::vector<LutInfoEx> bricklist = _partsNeeded(start, size, lod);
  if (bricklist.size() != 1)
    throw OpenZGY::Errors::ZgyInternalError("expeditedRead messed up.");
  const LutInfoEx& brick = bricklist.front();

  switch (brick.status) {
  case BrickStatus::Missing:
    if (_logger(2, ""))
      _logger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of missing brick");
    fillMe(data, size, this->_metadata->ih().defaultstorage(), result_type);
    break;

  case BrickStatus::Constant:
    if (_logger(2, ""))
      _logger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of constant brick");
    fillMe(data, size, brick.double_constvalue, result_type);
    break;

  case BrickStatus::Normal:
    if (brick.size_in_file != size[0] * size[1] * size[2] *
        static_cast<std::int64_t>(RawDataTypeDetails(result_type).size)) {
      throw OpenZGY::Errors::ZgyInternalError("Bad size in expeditedRead.");
    }
    if (_logger(2, ""))
      _logger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of regular brick");
    this->_file->xx_read(data, brick.offset_in_file, brick.size_in_file, UsageHint::Data);
    // TODO-Low: If _deliverOneBrick() clears padding samples, do likewise here.
    break;

  case BrickStatus::Compressed:
    // TODO-Performance: Also support decompression.
    if (_logger(2, ""))
      _logger(2, std::stringstream()
              << "Expedited read " << fmt(start) << " of compressed brick (not implemented)");
    return false;

  default:
    throw OpenZGY::Errors::ZgyInternalError("Internal error, bad brick status");
  }
  return true;
}


/**
 * \brief Mark bricks as dirty.
 *
 * \details
 * Mark all affected bricks as dirty, both those written now and those
 * at a higher lod level that will need to be re-generated.
 *    - lod=0, flag bit 0x01: This fullres brick written by the user.
 *    - lod>0, flag bit 0x01: Brick depends on a fullres that was written.
 *    - lod>0, flag bit 0x02: This lowres brick was written by genlod
 * The second one is redundant as it could be calculated on the fly,
 * and the third one might not be useful outside of debugging.
 *
 * Unlike trackedBricksDirty() the region has already been converted
 * to a list of bricks but the numbers are still in samples not bricks.
 */
void
ZgyInternalBulk::_trackedBricksSetDirty(
     const std::vector<index3_t>& work,
     int32_t lod)
{
  if (_modified_bricks.size() != 0) {
    const IInfoHeaderAccess& ih    = this->_metadata->ih();
    const index3_t           bs    = ih.bricksize();
    const std::int32_t       nlods = ih.nlods();

    int new_dirty{0}, old_dirty{0};
    // number of lods we will build, not the current count.
    for (const index3_t& pos : work) {
      std::int64_t ii = pos[0]/bs[0], jj = pos[1]/bs[1], kk = pos[2]/bs[2];
      for (std::int32_t mylod = lod; mylod < nlods; ++mylod) {
        const std::int64_t bix = LookupTable::getBrickLookupIndex
          (ii, jj, kk, mylod, ih.lodsizes(), ih.brickoffsets());
        if (bix >= (std::int64_t)_modified_bricks.size())
          throw OpenZGY::Errors::ZgyInternalError("overflow _modified_bricks.");
        std::uint8_t flagval = (lod==0) ? 1 : 2;
        if ((_modified_bricks[bix] & flagval) == 0)
          ++new_dirty;
        else
          ++old_dirty;
        _modified_bricks[bix] |= flagval;
        ii /= 2;
        jj /= 2;
        kk /= 2;
      }
    }
    if (_logger(5)) {
      auto clean = std::count_if(_modified_bricks.begin(),
                                 _modified_bricks.end(),
                                 [](std::uint8_t x){return (x&1) == 0;});
      _logger(5, std::stringstream()
              << "dirty: " << new_dirty << " new, " << old_dirty
              << " already soiled, " << clean << "/"
              << _modified_bricks.size() << " still clean");
    }
  }
}

/**
 * \brief For debugging, show dirty bricks.
 *
 * \details
 * The ascii art display shows brick-columns (combines all bricks in Z)
 * while the numbers printed after each LOD are per brick.
 */
void
ZgyInternalBulk::trackedBricksShowDirty(int loglevel) const
{
  if (!_logger(loglevel) || _modified_bricks.size() == 0)
    return;

  if (std::is_sorted(_modified_bricks.begin(), _modified_bricks.end()) &&
      _modified_bricks.front() == _modified_bricks.back()) {
    _logger(loglevel, std::stringstream()
            << "All " << _modified_bricks.size() << "bricks are "
            << (_modified_bricks.front() ? "dirty" : "clean"));
  }
  else {
    _logger(loglevel, "Dumping all tracked bricks:");
    const InternalZGY::IInfoHeaderAccess& ih = this->_metadata->ih();
    for (std::int64_t lod = ih.lodsizes().size() - 1; lod >= 0; --lod) {
      _logger(loglevel, "LOD " + std::to_string(lod));
      int num_dirty{0}, num_total{0};
      for (std::int64_t ii=0; ii<ih.lodsizes()[lod][0]; ++ii) {
        std::stringstream line;
        for (std::int64_t jj=0; jj<ih.lodsizes()[lod][1]; ++jj) {
          char flag = '.';
          for (std::int64_t kk=0; kk<ih.lodsizes()[lod][2]; ++kk) {
            ++num_total;
            std::int64_t index = InternalZGY::LookupTable::getBrickLookupIndex
              (ii, jj, kk, lod, ih.lodsizes(), ih.brickoffsets());
            if (index >= (std::int64_t)_modified_bricks.size()) {
              flag =  '?';    // Actually a fatal error
              break;
            }
            if (_modified_bricks[index] & 0x01) {
              ++num_dirty;
              flag = '*'; // written (if lod 0) or need rebuild (lod != 0)
            }
            else if (flag == '.' && (_modified_bricks[index] & 0x02) != 0) {
              flag = 'o'; // written, but only by genlod
            }
          }
          line << flag;
        }
        _logger(loglevel, line.str());    // After each inline
      }
      _logger(loglevel, std::stringstream()
              << "Total " << num_dirty << "/" << num_total
              << " dirty."); // After each LOD
    }
  }
}

/**
 * \brief Check if any brick in the region has been modified.
 *
 * \details
 * Used by incremental low resolution computation.
 * The input start and size is in sample coordinates.
 */
std::uint8_t
ZgyInternalBulk::trackedBricksDirty(
     const std::array<std::int64_t,3>& start,
     const std::array<std::int64_t,3>& size,
     int32_t lod) const
{
  if (_modified_bricks.size() == 0)
    return 0xFF;                // Or maybe treat this as an error.
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  std::uint8_t result{0};
  const std::array<std::int64_t,3> bs = ih.bricksize();
  const std::array<std::int64_t,3> bstart
    {start[0] / bs[0], start[1] / bs[1], start[2] / bs[2]};
  const std::array<std::int64_t,3> bend
    {(start[0] + size[0] + bs[0] - 1) / bs[0],
     (start[1] + size[1] + bs[1] - 1) / bs[1],
     (start[2] + size[2] + bs[2] - 1) / bs[2]};
  for (std::int64_t ii = bstart[0]; ii < bend[0]; ++ii) {
    for (std::int64_t jj = bstart[1]; jj < bend[1]; ++jj) {
      for (std::int64_t kk = bstart[2]; kk < bend[2]; ++kk) {
        const std::int64_t bix = LookupTable::getBrickLookupIndex
          (ii, jj, kk, lod, ih.lodsizes(), ih.brickoffsets());
        // "cannot happen" because getBrickLookupIndex() should have checked.
        if (bix >= (std::int64_t)_modified_bricks.size())
          throw OpenZGY::Errors::ZgyInternalError("overflow _modified_bricks.");
        result |= _modified_bricks[bix];
      }
    }
  }
  return result;
}

/**
 * Mitigate the problem that ZgyInternalBulk is used both for reading
 * and writing, making it difficult to keep track of thread safety.
 * The function below forces the code to be more explicit read vs write.
 *
 * See ZgyInternalBulk::ZgyInternalBulk for more details.
 */
std::shared_ptr<ZgyInternalMeta>
ZgyInternalBulk::_get_metadata_rw()
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
ZgyInternalBulk::_logger(int priority, const std::string& message) const
{
  return _loggerfn(priority, message);
}

/**
 * Convenience for invoking _loggerfn with a stringstream.
 * Due to a somewhat naughty cast, the function can be caller as:
 *
 *   if(_logger(pr1))
 *    _logger(pri, std::stringstream() << some << data << here);
 *
 * The first line is optional. It just prevents the expression in
 * the second line from being evaluatet if debugging is disabled.
 */
bool
ZgyInternalBulk::_logger(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return _logger(priority, sstream ? sstream->str() : std::string());
}

/**
 * Check that i,j,k,lod are all inside valid bounds. Throw if not.
 * If used to validate an alpha tile then k should be passed as 0.
 * This is similar to LookupTable::_validatePosition() but checks
 * the entire region, not just one brick at a time. And the input
 * is in trace numbers not brick numbers.
 *
 * The exception message is currently rather verbose.
 * I might want to shorten it a bit; skipping the actual position.
 *
 * Note that a file currently being written to but not finalized yet
 * will see nlods() == 1 (which is entirely correct) but to avoid
 * some chicken and egg problems it will still be permitted to read
 * back higher LOD levels. See ZgyInternalBulk::_validateUserPosition.
 * When an incomplete file is opened read only, attempts to read low
 * resolution data will throw an exception.
 *
 * Note on test coverage: For defensive coding these checks are done
 * both here at a high level and in lookuptable. This hurts coverage
 * since it might not be possible to get those other checks to fail.
 * And some errors may be caught even earlier e.g. if trying to make
 * a DataBuffer with negative size.
 */
void
ZgyInternalBulk::_validateUserPosition(
    const std::array<std::int64_t,3>& start,
    const std::array<std::int64_t,3>& size,
    int32_t lod) const
{
  std::string error;
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const std::array<std::int64_t,3> one{1,1,1} ;
  const std::array<std::int64_t,3> bs  = this->_metadata->ih().bricksize();
  const std::array<std::int64_t,3> end = ((start + size + bs - one) / bs) * bs;
  const std::int32_t nlods = static_cast<std::int32_t>(ih.lodsizes().size());
  const bool open_for_write = this->_metadata_rw != nullptr;
  // TODO-Low: Performance: Cache usable_nlods.
  const std::int32_t usable_nlods = open_for_write ? nlods :
    LookupTable::usableBrickLOD
    (this->_metadata->ih().lodsizes(),
     this->_metadata->ih().brickoffsets(),
     this->_metadata->blup().lup(),
     this->_metadata->blup().lupend());
  if (lod < 0 || lod >= usable_nlods) {
    std::stringstream ss;
    ss << "Requested lod " << lod <<
      " is outside the valid range 0 to " << usable_nlods - 1 << " inclusive";
    _logger(1, ss.str() + "\n");
    throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
  const std::array<std::int64_t,3>& ssize = ih.lodsizes()[lod] * bs;
  if (start[0] < 0 || end[0] > ssize[0] || size[0] <= 0 ||
      start[1] < 0 || end[1] > ssize[1] || size[1] <= 0 ||
      start[2] < 0 || end[2] > ssize[2] || size[2] <= 0) {
    std::stringstream ss;
    ss << "Requested region"
       << " from (" << start[0] << ", " << start[1] << ", " << start[2] << ")"
       << " to (" << end[0] << ", " << end[1] << ", " << end[2] << ")"
       << " lod " << lod
       << " is empty or outside the valid range"
       << " (0, 0, 0)"
       << " to (" << ssize[0] << ", " << ssize[1] << ", " << ssize[2] << ")"
       ;
    _logger(1, ss.str() + "\n");
    throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
}

/**
 * Scale integral data from storage according to the coding range.
 * The operation is a no-op if file_dtype is float. The data must
 * be known to be in storage domain, we don't try to guess based on
 * the DataBuffer's valuetype.
 */
std::shared_ptr<DataBuffer>
ZgyInternalBulk::_scaleDataToFloat(
      const std::shared_ptr<DataBuffer>& data) const
{
  // The data buffer type should match file_dtype but I won't bother to check.
  // Also, if result ends up empty this means no conversion is needed.
  std::array<double,2> factors = this->_metadata->ih().storagetofloat();
  std::shared_ptr<DataBuffer> result = data->scaleToFloat(factors);
  return result ? result : data;
}

std::shared_ptr<DataBuffer>
ZgyInternalBulk::_scaleDataToStorage(
      const std::shared_ptr<DataBuffer>& data) const
{
  std::array<double,2> factors = this->_metadata->ih().storagetofloat();
  std::shared_ptr<DataBuffer> result =
    data->scaleToStorage(factors, _metadata->ih().datatype());
  return result ? result : data;
}

namespace {
  template<typename T> double _decodeConstantT(std::uint32_t in)
  {
    if (in == 1) // ZGY v1 could only store constant-zero.
      return 0;
    T value;
    byteswapT(&in);
    memcpy(&value, &in, sizeof(T));
    byteswapT(&value);
    return static_cast<double>(value);
  }

  template<typename T> std::uint32_t _encodeConstantT(double in)
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
 * The same reasoning applies to the inverse _encodeConstantT().
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
ZgyInternalBulk::_decodeConstant(std::uint32_t in, RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::SignedInt8:    return _decodeConstantT<std::int8_t>(in);
  case RawDataType::UnsignedInt8:  return _decodeConstantT<std::uint8_t>(in);
  case RawDataType::SignedInt16:   return _decodeConstantT<std::int16_t>(in);
  case RawDataType::UnsignedInt16: return _decodeConstantT<std::uint16_t>(in);
  case RawDataType::SignedInt32:   return _decodeConstantT<std::int32_t>(in);
  case RawDataType::UnsignedInt32: return _decodeConstantT<std::uint32_t>(in);
  case RawDataType::Float32:       return _decodeConstantT<float>(in);
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Inverse of _decodeConstant(). See that function for details.
 */
/*static*/ std::int32_t
ZgyInternalBulk::_encodeConstant(double in, RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::SignedInt8:    return _encodeConstantT<std::int8_t>(in);
  case RawDataType::UnsignedInt8:  return _encodeConstantT<std::uint8_t>(in);
  case RawDataType::SignedInt16:   return _encodeConstantT<std::int16_t>(in);
  case RawDataType::UnsignedInt16: return _encodeConstantT<std::uint16_t>(in);
  case RawDataType::SignedInt32:   return _encodeConstantT<std::int32_t>(in);
  case RawDataType::UnsignedInt32: return _encodeConstantT<std::uint32_t>(in);
  case RawDataType::Float32:       return _encodeConstantT<float>(in);
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
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
ZgyInternalBulk::_partsNeeded(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const
{
  std::vector<LutInfoEx> bricklist;
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const std::array<std::int64_t,3> bs = this->_metadata->ih().bricksize();
  const ILookupTableAccess& blup = this->_metadata->blup();
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
        IJK pos{iter[0], iter[1], iter[2], bs[0], bs[1], bs[2]};
        LookupTable::LutInfo info = LookupTable::getBrickFilePosition
          (pos.i0/pos.ni, pos.j0/pos.nj, pos.k0/pos.nk, lod,
           /**/ih.lodsizes(), ih.brickoffsets(), blup.lup(), blup.lupend(),
           ih.bytesperbrick());
        bricklist.push_back
          (LutInfoEx(info, pos,
                     _decodeConstant(info.raw_constant, ih.datatype())));
      }
    }
  }
  return bricklist;
}

/**
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
ZgyInternalBulk::_deliverOneBrick(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& startpos,
      const ReadRequest::data_t& raw, std::int64_t rawsize,
      BrickStatus brickstatus, bool as_float) const
{
  if (_logger(5))
    _logger(5, std::stringstream()
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
    data = _scaleDataToFloat(data);

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
ZgyInternalBulk::_usedPartOfBrick(
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
ZgyInternalBulk::_isUsedPartOfBrickAllConstant(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& brickpos,
      int32_t lod) const
{
  const std::array<std::int64_t,3> used =
    this->_usedPartOfBrick(data->size3d(), brickpos, lod);
  return data->isAllSame(used.data());
}

/**
 * Pad unused parts of the data buffer by replicating the last samples,
 * but only up to a multiple of 'modulo' samples. Handles just one
 * dimension, so caller will typically invoke us three times.
 */
/*static*/ void
ZgyInternalBulk::_setPaddingToEdge(
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
    data->copyFrom(data.get(), srcorig.data(), dstorig.data(), cpyorig.data(), cpysize.data());
    ++slice;
  }
}

/**
 * Pad unused parts of the data buffer with a constant.
 * Handles just one dimension, so caller should invoke us three times.
 */
/*static*/ void
ZgyInternalBulk::_setPaddingToConst(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, int dim)
{
  if (used[dim] > 0 && used[dim] < data->size3d()[dim]) {
    std::array<std::int64_t,3> srcorig{0,0,0};
    std::array<std::int64_t,3> dstorig{0,0,0};
    std::array<std::int64_t,3> cpyorig{0,0,0};
    std::array<std::int64_t,3> cpysize=data->size3d();
    cpyorig[dim] = used[dim];
    cpysize[dim] -= used[dim];

  std::shared_ptr<DataBuffer> filler =
    DataBuffer::makeScalarBuffer3d(missingvalue, data->size3d(), data->datatype());
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
ZgyInternalBulk::_setPaddingSamples(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, const compressor_t& compressor)
{
  // Currently not possible because "compressor" is a
  // std::function which means it has just one entry point.
  //if (compressor)
  //  return compressor.pad(data, used);

  for (int dim=0; dim<3; ++dim)
    _setPaddingToConst(data, used, missingvalue, dim);

  for (int dim=0; dim<3; ++dim)
    _setPaddingToEdge(data, used, 4, dim);
}

/**
 * Done in _writeOneNormalBrick() which is called first: Apply padding,
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
ZgyInternalBulk::_writeWithRetry(const WriteNowArgPack& args)
{
  if (_logger(2))
    _logger(2, "_writeWithRetry(" + args.toString() + ")\n");
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
      fileoffset = 0; // Write a new brick, abandoning the old one.
    }
  }
  if (!fileoffset) { // New brick or abandoned brick
    fileoffset = this->_file->xx_eof();
    if (_logger(2))
      _logger(2, std::stringstream()
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
       &this->_get_metadata_rw()->blup().lup(),
       &this->_get_metadata_rw()->blup().lupend());
  }
}

/**
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
ZgyInternalBulk::_writeOneNormalBrick(const WriteBrickArgPack& args)
{
  const double defaultstorage = this->_metadata->ih().defaultstorage();
  if (_logger(2))
    _logger(2, "_writeOneNormalBrick(" + args.toString() + ")\n");

  if (args.data->isScalar())
    throw OpenZGY::Errors::ZgyInternalError("Wrong brick type in _writeOneNormalBrick");
  if (args.data->size3d() != this->_metadata->ih().bricksize())
    throw OpenZGY::Errors::ZgyInternalError("Wrong brick size in _writeOneNormalBrick");

  std::array<std::int64_t,3> used =
    this->_usedPartOfBrick(args.data->size3d(), args.brickpos, args.lod);
  this->_setPaddingSamples(args.data, used, defaultstorage, args.compressor);

  // rawdata_t::first is a smart pointer and it needs to be in case it
  // is later replaced with a buffer of compressed data. Besides, it is
  // convenient if we decide to queue the write request.
  rawdata_t rawdata(args.data->voidData(), args.data->allocsize() * args.data->itemsize());

  BrickStatus brickstatus = BrickStatus::Normal;
  if (args.compressor) {
    _logger(5, "Attempting compression");
    rawdata_t cdata = args.compressor(rawdata, args.data->size3d());
    if (cdata.first) {
      _logger(5, "compression successful");
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
  // Arguments that our caller needs to be pass to _writeWithRetry()
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
 * Apply conversion of the constant value to a lookup table entry
 * Pass on the write request to the function that updates the lut.
 *
 * The size stored in "data" is ignored. The constant applies to
 * the whole brick. Also, args.fileoffset and args.compressor are N/A.
 *
 * Args: brickpos, lod, data, compressor(N/A), fileoffset(N/A).
 */
void
ZgyInternalBulk::_writeOneConstantBrick(const WriteBrickArgPack& args)
{
  if (_logger(2))
    _logger(2, "_writeOneConstantBrick(" + args.toString() + ")\n");
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  ILookupTableAccess& blup = this->_get_metadata_rw()->blup();
  LookupTable::LutInfo info
    (BrickStatus::Constant, /*offset=*/0, /*size=*/0,
     _encodeConstant(args.data->scalarAsDouble(), ih.datatype()));
  LookupTable::setBrickFilePosition
    (args.brickpos[0], args.brickpos[1], args.brickpos[2], args.lod, info,
     ih.lodsizes(), ih.brickoffsets(), &blup.lup(), &blup.lupend());
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
 * in _writeWithRetry() when a ZgySegmentIsClosed is caught.
 * The update mode doesn't need to be checked again at that point
 * unless we want to support UpdateMode.NoLeaks.
 */
bool
ZgyInternalBulk::_mustLeakOldBrick(
      const std::shared_ptr<DataBuffer>& data,
      const compressor_t& compressor,
      BrickStatus brickstatus) const
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
  if (error) {
    // TODO-Low: Display symbolic names for enums.
    std::stringstream ss;
    ss << "Updating a " << (this->_file->xx_iscloud() ? "cloud" : "local")
       << " BrickStatus::" << (int)brickstatus << " brick with "
       << (compressor ? "compressed" : data->isScalar() ? "constant" : "normal")
       << " data is illegal in UpdateMode::" << (int)this->_update_mode;
    throw OpenZGY::Errors::ZgyUserError(ss.str());
  }
  return leak;
}

/**
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
ZgyInternalBulk::_writeOneBrick(const WriteBrickArgPack& args_in)
{
  WriteBrickArgPack args(args_in);
  // brickpos and lod will not change.
  // data, compressor, fileoffset might all change but usually don't.

  // TODO-Low note copy/pasted code from _partsNeeded().
  // On a more general note, would it have been better to
  // call _partsNeeded() early on, and have this function
  // and a few others operate on lists? This would be more
  // complicated but might allow better parallelization.
  const IInfoHeaderAccess& ih = this->_metadata->ih();
  const ILookupTableAccess& blup = this->_metadata->blup();
  LookupTable::LutInfo info = LookupTable::getBrickFilePosition
    (args.brickpos[0], args.brickpos[1], args.brickpos[2], args.lod,
     ih.lodsizes(), ih.brickoffsets(), blup.lup(), blup.lupend(),
     ih.bytesperbrick());

  // We can easily pretend the brick did not exist and needs to be allocated.
  if (_mustLeakOldBrick(args.data, args.compressor, info.status))
    info = LookupTable::LutInfo(BrickStatus::Missing, 0, 0, 0);

  std::array<std::int64_t,3> bs = ih.bricksize();
  const bool data_const = args.data->isScalar();
  const bool file_const = (info.status != BrickStatus::Normal &&
                           info.status != BrickStatus::Compressed);

  if (file_const) { // also true if never written yet.
    if (data_const) {
      // Caller asked to store a constant value.
      _logger(2, "Explicit store constant.\n");
    }
    else if (_isUsedPartOfBrickAllConstant(args.data, args.brickpos, args.lod)) {
      // Caller did not explicitly ask for a constant value,
      // but all useable samples (i.e. all samples that are
      // inside the survey boundaries) have the same value.
      _logger(2, "Implicit store constant.\n");
      // Need to change to a real scalar. This is a bit roundabout.
      // TODO-Low can I have a DataBuffer::makeSimilar() that retains
      // data type and NDim but can turn it into a scalar or
      // create a deep copy or maybe an uninitialized copy.
      // Or the other way, create an inflated brick from a constant.
      double value = args.data->scalarAsDouble();
      args.data = DataBuffer::makeScalarBuffer3d(value, bs, args.data->datatype());
    }
    else {
      _logger(2, "Store a new brick.\n");
      args.fileoffset = 0;  // Says to allocate a new block on file.
    }
  }
  else { // !file_const
    if (data_const) {
      _logger(1, "Const data expanded before update of brick.\n");
      // The brick has already been allocated. Cannot set it to
      // constant-value because the file storage for this brick
      // would then be leaked. Cannot compress it either.
      // The assumption is that the old brick was uncompressed,
      // because if it had been compressed then _mustLeakOldBrick()
      // would have returned true and we wouldn't have gotten here.
      //
      // TODO-Medium: Fragile code:
      //
      //   * _mustLeakOldBrick() might change to allow overwriting
      //     a compressed brick with a smaller one. If so, must
      //     *not* veto the compression.
      //
      //   * _mustLeakOldBrick() will always, if compressor is set,
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
      _logger(1, "Update a brick.\n");
      args.fileoffset = info.offset_in_file; // Re-use if possible.
      args.compressor = compressor_t();      // No compress on update.
    }
  }
  // Next, _writeOneNormalBrick or _writeOneConstantBrick
  // (depending on the type that args.data has now)
  // should be called with the modified argument pack.
  // TODO-Low: Don't copy args so often.
  return std::make_shared<const WriteBrickArgPack>(args);
}

/**
 * Apply splitting into bricks.
 * Expand bricks at survey edge to full size unless r/m/w already did
 * Convert sample position to brick number. I.e. divide by bricksize.
 * Pass on the write request, one brick at a time, to the next level.
 *
 * The data must already be converted to storage. Start must be a
 * multiple of brick size. start + size must either be a multiple of
 * block size or set to the end of the survey. Sizes refer to the
 * specified lod level. At lod level N the valid range both for size
 * and bricksize is half (rounded up) of that in lod N-1.
 *
 * With current usage the data may or may not have been padded to
 * the brick size depending on whether a read/modify/write was used.
 */
void
ZgyInternalBulk::_writeAlignedRegion(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      std::int32_t lod,
      const compressor_t& compressor)
{
  SimpleTimerEx t1(*_ptimer_mt);

  const double   defaultstorage = this->_metadata->ih().defaultstorage();
  const std::int64_t lodfactor  = static_cast<std::int64_t>(1) << lod;
  const IInfoHeaderAccess& ih   = this->_metadata->ih();
  const index3_t    bs          = ih.bricksize();
  const RawDataType file_dtype  = ih.datatype();
  const index3_t    survey_beg    {0,0,0};
  const index3_t    survey_size = (ih.size() + (lodfactor - 1)) / lodfactor;
  const index3_t    count       = data->size3d();
  std::vector<index3_t> work;
  {
    const index3_t beg_brick      = (start / bs) * bs;
    const index3_t end_brick      = ((start + count + bs - (std::int64_t)1) / bs) * bs;
    index3_t it;
    for (it[0] = beg_brick[0]; it[0] < end_brick[0]; it[0] += bs[0])
      for (it[1] = beg_brick[1]; it[1] < end_brick[1]; it[1] += bs[1])
        for (it[2] = beg_brick[2]; it[2] < end_brick[2]; it[2] += bs[2])
          work.push_back(it);
  }
  const std::size_t worksize = work.size();

  if (_logger(2))
    _logger(2, std::stringstream()
            << "_writeAlignedRegion("
            << "start " << fmt(start)
            << ", size " << fmt(count)
            << ", type " << (int)file_dtype
            << ")\n");

  // Mark all affected bricks as dirty, both those written now and those
  // at a higher lod level that will need to be re-generated.
  this->_trackedBricksSetDirty(work, lod);

  // If the application asked to write just a single brick than we can
  // skip the multithreading logic and in most cases skip the copying
  // into the temp buffer. That shortcut is probably not worth the
  // risk. Especially if the buffer belongs to the application and not
  // something reallocated as part of the r/m/w processing. Because we
  // cannot control the lifecycle of the user's data buffer. It is a
  // shared_ptr but it will typically have a bogus destructor. Also it
  // cannot be modified in place for e.g. byteswapping. Allowing this
  // buffer in the lower layers creates all kinds of corner cases. And
  // how often will we hit the shortcut case anyway?
  //
  // One temporary brick buffer is allocated for each brick to be written,
  // consuming the same amount of data as the original request had.
  // Since the requests are queued up the code may also end up allocating
  // memory for compressed data etc. and perhaps not freeing the not
  // compressed buffer immediately.
  //
  // If this becomes problematic then consider:
  //
  //   * If multi-threading is disabled at runtime then the "brick" buffer
  //     can be re-used. _writeOneConstantBrick() or _writeWithRetry(*it)
  //     would be called inside the loop, one brick at a time.
  //
  //   * Instead of processing every single brick and then writing all
  //     at the end, consider splitting huge requests into several
  //     slightly less huge requests and process those in sequence.
  //     (the same applies to the r/m/w logic).
  //
  // If the user-provided data is a scalar then it can almost be used as-is
  // and using the same brick for every write. But the size needs to be
  // adjusted to become just a single brick. _writeOneConstantBrick()
  // ignores the size but _writeOneBrick() etc. might not.
  std::shared_ptr<DataBuffer> constbrick;
  if (data->isScalar()) {
    double value = data->scalarAsDouble();
    constbrick = DataBuffer::makeScalarBuffer3d(value, bs, data->datatype());
  }

  std::vector<std::shared_ptr<const WriteBrickArgPack>> const_queue(worksize);
  std::vector<std::shared_ptr<const WriteNowArgPack>>   normal_queue(worksize);
  int numthreads = std::min(omp_get_max_threads(), std::max((int)worksize, 1));
  MTGuard guard("copy-in", numthreads);
#pragma omp parallel for num_threads(numthreads) if(enable_compress_mt() && worksize > 1)
  for (std::int64_t ix = 0; ix < static_cast<std::int64_t>(worksize); ++ix) {
    SimpleTimerEx t3(*_ptimer);
    guard.run([&](){
      const index3_t surveypos = work[ix]; // user's start i0,j0,k0 rounded down
      const index3_t brickpos = work[ix] / bs; // as above, but in brick coords
      std::shared_ptr<DataBuffer> brick = constbrick;
      if (!brick) {
        brick = DataBuffer::makeNewBuffer3d(bs, data->datatype());
        // TODO-Performance, any way of avoiding the fill()?
        // TODO-Medium: If the existing brick has a constvalue then
        // initialize with that value insteaf of defaultstorage.
        // The reason is that if the application filled the entire
        // survey with a const value then we should respect that.
        brick->fill(defaultstorage);
        brick->copyFrom(data.get(), // source
                        start.data(), surveypos.data(), // in survey coords
                        survey_beg.data(), survey_size.data()); // clip to srv
      }

      std::shared_ptr<const WriteBrickArgPack> args =
        std::make_shared<const WriteBrickArgPack>
        (brickpos, lod, brick, compressor, 0);
      args = _writeOneBrick(*args);
      if (args->data->isScalar()) {
#pragma omp critical // paranoia?
        const_queue[ix] = args;
      }
      else {
        std::shared_ptr<const WriteNowArgPack> now =_writeOneNormalBrick(*args);
#pragma omp critical // paranoia?
        normal_queue[ix] = now;
      }
    });
  } // end parallel for loop
  guard.finished();

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

  // The _writeWithRetry() method is not threadsafe and *must* be
  // called sequentially. It is also highly recommended to write
  // the data in the same order as it occurred in the work array.
  // Otherwise the bricks get scrambled and reading will be less
  // efficient. This limitation prevents us from having a separate
  // worker process draining the normal_queue as data becomes
  // available. _writeOneConstantBrick() doesn't care about order
  // but might as well follow the same rules. That function is
  // very lightweight.

  t1.stop();
  SimpleTimerEx t2(*_ptimer_st);
  ErrorsWillCorruptFile watchdog(this);
  for (const auto& it : const_queue)
    if (it)
      _writeOneConstantBrick(*it);
  for (const auto& it : normal_queue)
    if (it)
      _writeWithRetry(*it);
  watchdog.disarm();
}

/*
 * \brief Write an arbitrary region.
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
  _validateUserPosition(start_in, data_in->size3d(), lod);
  std::shared_ptr<DataBuffer> data = data_in;
  std::array<std::int64_t,3> start = start_in;

  if (data->datatype() != (!is_storage ? RawDataType::Float32 :
                           _metadata->ih().datatype()))
    throw OpenZGY::Errors::ZgyUserError("Invalid data type in writeRegion");

  if (_logger(2))
    _logger(2, std::stringstream()
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

  // TODO-Performance: Combining range() and _scaleDataToStorage()
  // might save some time.

  if (!is_storage) {
    SimpleTimerEx t1(*_ststimer);
    data = _scaleDataToStorage(data);
  }
  // If this write or (more likely) writeconst covers the entire survey,
  // Start over with keeping track of range of samples ever written.
  // Also turn off the logic to keep track of stats and histogram
  // because a full rebuild will be needed anyway.
  // This case is actually quite useful. The application can issue a single
  // writeconst() for the entire survey to avoid having to guess what the
  // default value will be. This can also prevent getting some arbitrary
  // (usually zero) value added to the statistical and histogram range even
  // if it doesn't exist in the samples.
  const bool entire_survey =
    covers(start, data->size3d(), std::array<std::int64_t,3>{0,0,0}, _metadata->ih().size());
  if (entire_survey) {
    _written_sample_min = std::numeric_limits<double>::infinity();
    _written_sample_max = -std::numeric_limits<double>::infinity();
    trackedBricksTryEnable(false);
  }

  const index3_t beg = start;
  const index3_t end = beg + data->size3d();
  const index3_t bs = this->_metadata->ih().bricksize();
  const index3_t survey_beg{0,0,0};
  const std::int64_t lodfactor = static_cast<std::int64_t>(1) << lod;
  const index3_t survey_size = ((this->_metadata->ih().size()) + (lodfactor-1)) / lodfactor;
  const index3_t survey_end = survey_beg + survey_size;
  // Note: std::min(array1,array2) doesn't do what you think it does;
  // it calls a lexiographical compare not an element-wise min().
  index3_t valid_count{
        std::min(end[0], survey_end[0]) - beg[0],
        std::min(end[1], survey_end[1]) - beg[1],
        std::min(end[2], survey_end[2]) - beg[2]
    };

  // The previous contents of the bricks will be needed if a read/modify/write
  // is indicated or if changes to statistics are to be tracked. For r/m/w
  // not all the old bricks might be needed. But for now read all of them.
  bool need_rmw = false;
  for (int d = 0; d < 3; ++d)
    if ((beg[d]%bs[d]) != 0 || ((end[d]%bs[d]) != 0 && end[d] < survey_end[d]))
      need_rmw = true;
  if (_modified_stats && lod==0)
    need_rmw = true;
  if (need_rmw) {
    // Since we are doing a read/modify/write, also read any padding
    // outside the survey. _writeAlignedRegion will see full
    // bricks also at the end of the survey. If no r/m/w needed
    // the latter is not guaranteed. TODO-Worry does this apply also to C++?
    index3_t new_start = (beg / bs) * bs;
    index3_t new_count = (((end + bs - (std::int64_t)1) / bs) * bs) - new_start;
    std::shared_ptr<DataBuffer> new_data =
      DataBuffer::makeNewBuffer3d(new_count, this->_metadata->ih().datatype());
    // Note: std::min(array1,array2) doesn't do what you think it does;
    // it calls a lexiographical compare not an element-wise min().
    const index3_t new_valid{
        std::min(new_start[0] + new_count[0], survey_end[0]) - new_start[0],
        std::min(new_start[1] + new_count[1], survey_end[1]) - new_start[1],
        std::min(new_start[2] + new_count[2], survey_end[2]) - new_start[2]
    };
    if (_logger(1))
      _logger(1, std::stringstream()
              << "_writeRegion r/m/w:"
              << " user" << fmt(start) << " " << fmt(data->size3d())
              << " padded " << fmt(new_start) << " " << fmt(new_count) << "\n");

    readToExistingBuffer(new_data, new_start, lod, /*as_float=*/false);
    if (lod == 0)
      trackChanges(new_data, new_valid, _modified_stats.get(), _modified_histo.get(), /*add=*/false);
    new_data->copyFrom(data.get(), // source
                       start.data(), new_start.data(), // in survey coords
                       survey_beg.data(), survey_size.data()); // clip to survey
    if (lod == 0)
      trackChanges(new_data, new_valid, _modified_stats.get(), _modified_histo.get(), /*add=*/true);
    data = new_data;
    start = new_start;
    valid_count = new_valid;
  }
  else if (entire_survey) {
    _logger(2, "_writeRegion direct: entire survey");
  }
  else {
    if (_logger(2))
      _logger(2, std::stringstream()
              << "_writeRegion direct:"
              << " user" << fmt(start) << " " << fmt(data->size3d()) << "\n");
  }

  // Keep track of valuerange written, to help generate a histogram.
  // This also considers padding samples. Assuming that even though
  // the values are unspecified they won't be outside the value range
  // of the real data plus "defaultvalue".
  // TODO-Worry: Safer to only check the used part of buffer.
  // The issue is the same as the one for trackChanges().
  // Or somehow replace all padding samples with defaultvalue.
  // The range is in storage values. So e.g. for int8 it will be no
  // more than [-128..+127].
  // Note that range() can be expensive. I have seen it take 5% of
  // total time in a copy application. But a single call to range()
  // is often fast enough to not warrant an OpenMP loop. See range()
  // for more detals.


  // TODO-Worry: read/modify/write processing past the survey edge may
  // cause the defaultvalue to be included in the histogram range anyway.
  // See note in _writeAlignedRegion().
  if (lod == 0) {
    if (_modified_stats && _modified_stats->getmin() <= _modified_stats->getmax()) {
      // _modified_stats refers to everything ever written, not just what
      // data->range() would have given, but here it works just as well.
      _written_sample_min = std::min(_written_sample_min, _modified_stats->getmin());
      _written_sample_max = std::max(_written_sample_max, _modified_stats->getmax());
    }
    else {
      std::pair<double, double> minmax = data->range(valid_count.data());
      if (minmax.first <= minmax.second) {
        _written_sample_min = std::min(_written_sample_min, minmax.first);
        _written_sample_max = std::max(_written_sample_max, minmax.second);
      }
    }
  }

  _writeAlignedRegion(data, start, lod, compressor);
}

} // namespace
