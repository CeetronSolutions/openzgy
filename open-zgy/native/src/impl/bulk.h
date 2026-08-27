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

#include "../declspec.h"
#include "enum.h"
#include "file.h" // Only for the ReadRequest::data_t typedef
#include "fancy_timers.h"

#include <memory>
#include <functional>
#include <array>
#include <vector>
#include <tuple>
#include <sstream>
#include <atomic>
#include <mutex>

/**
 * \file bulk.h
 * \brief Bulk data read/write.
 */

namespace Test { class TestBulk; }
namespace InternalZGY {
#if 0
}
#endif

class IFileADT;
class ZgyInternalMeta;
class DataBuffer;
struct IJK;
struct LutInfoEx;
struct WriteBrickArgPack;
struct WriteNowArgPack;
class SummaryPrintingTimerEx;
class StatisticData;
class HistogramData;
class TrackTouched;
class GenLodC;

/**
 * \brief Performance timers for bulk access.
 *
 * \details
 *
 * Collect the timers in a single struct to make them somewhat
 * less messy. Currently, all these timers are associated with
 * a ZgyInternalBulk instance. Free functions and static methods
 * will need to have the BulkAdHocTimers instance passed in.
 *
 * The struct can be aggregated into class ZgyInternalBulk
 * instead of using a smart pointer. This saves a few cycles.
 * Which might be significant when timers are not enabled.
 */
struct OPENZGY_TEST_API BulkAdHocTimers
{
  SummaryPrintingTimerEx st;
  SummaryPrintingTimerEx mt;
  SummaryPrintingTimerEx mt1;
  SummaryPrintingTimerEx mt2;
  SummaryPrintingTimerEx mt3;
  SummaryPrintingTimerEx mt4;
  SummaryPrintingTimerEx mt5;
  SummaryPrintingTimerEx mt6;
  SummaryPrintingTimerEx mt7;
  SummaryPrintingTimerEx mt8;
  SummaryPrintingTimerEx mt9;
  SummaryPrintingTimerEx mt10;
  SummaryPrintingTimerEx mt11;
  SummaryPrintingTimerEx mt12;
  SummaryPrintingTimerEx mt13;
  SummaryPrintingTimerEx mt14;
  SummaryPrintingTimerEx mt15;
  SummaryPrintingTimerEx mt16;
  SummaryPrintingTimerEx gl;
  SummaryPrintingTimerEx ptimer;
  SummaryPrintingTimerEx sts;

  BulkAdHocTimers(const BulkAdHocTimers&)  = delete;
  BulkAdHocTimers(const BulkAdHocTimers&&) = delete;
  BulkAdHocTimers& operator=(const BulkAdHocTimers&)  = delete;
  BulkAdHocTimers& operator=(const BulkAdHocTimers&&) = delete;

  BulkAdHocTimers()
  : st("writeAligned[S]")
  , mt("writeAligned[M]")
  , mt1("writeAligned[M1]") // setup
  , mt2("writeAligned[M2]") // read old contents
  , mt3("writeAligned[M3]") // outside mt loop, see "writeAligned" for inside
  , mt4("writeAligned[M4]") // lock wait
  , mt5("writeAligned[M5]") // add/sub stats/histo
  , mt6("writeAligned[M6]") // inside mt loop, only subtract old stats/histo
  , mt7("writeAligned[M7]") // inside mt loop, only r/m/w processing
  , mt8("writeAligned[M8]") // inside mt loop, only add new stats/histo
  , mt9("writeAligned[M9]") // inside mt loop, prepare write one brick
  , mt10("trkChangesT[M10]") // setup. Note, M6+M8 ~= M11+M12+M13
  , mt11("trkChangesT[M11]") // subtract,add scalar value
  , mt12("trkChangesT[M12]") // subtract,add full brick
  , mt13("trkChangesT[M13]") // subtract,add edge brick
  , mt14("trkChangesT[M14]") // lock wait
  , mt15("trkChangesT[M15]") // add/sub stats/histo holding lock
  , mt16("trkChangesT[M16]") // add/sub stats/histo not holding lock
  , gl("writeAligned[GL]")
  , ptimer("writeAligned[in]")
  , sts("scaleToStorage")
  {
  }
};

/**
 * Read or write bulk data. The meta data needs to have been read
 * already. The user-callable API will forward its read requests here.
 *
 * Thread safety:
 * \li most data members are only used for write.
 *     Those don't need to be thread safe.
 * \li this->_file points to an IFileADT which is already thread safe
 *     where needed.
 * \li this->_metadata is problematic because it is too easy to
 *     accidentally invoke a method that changes data even if the
 *     file is open for read only. Mitigated by separate const
 *     and mutable pointers to metadata. See ZgyInternalBulk::ZgyInternalBulk.
 *
 * SINGLEPASS: Data members and methods removed from ZgyInternalBulk.
 * The list might not be complete. The names might linger on in comments.
 *
 *   _written_sample_min _written_sample_max
 *   _modified_bricks _modified_stats _modified_histo
 *
 *   _writeAlignedRegion() _calculate() _read() _write()
 *   readToNewBuffer() isCloud() valueRangeWritten() trackedChanges()
 *   trackedBricks() trackedBricksTryEnable() trackedBricksDirty()
 *   trackedBricksShowDirty() _trackedBricksSetDirty()
 *   addCompressionNoise()
 */
class OPENZGY_TEST_API ZgyInternalBulk
{
  friend class Test::TestBulk;
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

private:
  class ErrorsWillCorruptFile;
  std::shared_ptr<IFileADT>         _file;
  std::shared_ptr<const ZgyInternalMeta> _metadata;
  std::shared_ptr<ZgyInternalMeta> _metadata_rw;
  std::shared_ptr<ZgyInternalMeta> get_metadata_rw();
  UpdateMode _update_mode;
  bool       _compressed_write; // If true: do not align bricks.
  std::atomic<bool> _is_bad;    // If true: instance not usable.
  LoggerFn   _loggerfn;

  // The mutex is used inside trackChanges to protect new_xxx;
  // Other access to new_xxx should be single threaded and
  // don't need to set the lock.
  std::mutex                        _new_mutex;
  // The _new_xxx members are initialized in newTrackedBricksTryEnable()
  // which *should* be called after a file has been opened for write,
  // or opened for update, or if a call to writeconst reset the entire file.
  //
  //  - In ZgyInternalBulk constructor.
  //  - In writeRegion() when resetting the entire file
  //  - NOT in finalize.
  //
  // Updates to stats and histo is done in trackChanges called from
  // ZgyInternalBulk::writeAlignedBrickList(). Note that trackChanges()
  // is also called from ZgyInternalBulk::writeRegion(), but there it
  // is used to update the old members.
  //
  std::shared_ptr<TrackTouched>     _new_modified_bricks;
  std::shared_ptr<StatisticData>    _new_modified_stats;
  std::shared_ptr<HistogramData>    _new_modified_histo;
  bool                              _new_stathist_good; // not updating bad file

  // Long lived GenLod instance, to be invoked every time we think we might
  // be able to generate more lowres bricks. Still undecided whether this
  // belongs here, or whether ZgyInternalBulk::writeAlignedBrickList()
  // should instead create a short lived instance. One benefit of having
  // it here is that the instance keeps state that is useful because
  // it isn't readily accessed from ZgyInternalBulk, only from the API layer.
  std::shared_ptr<GenLodC> _new_genlod;

  // Remembering the user's choice of LOD generation. Sort of belongs
  // here, but could also have been inside the genlod instance.
  InternalLodMode _new_internal_lod_mode;

  // Remembering the user's requested histogram range. Currently only honored
  // for float files, and only when initially created or reset to all const.
  std::array<float,2> _force_histo_range;

  // Performance measurements. Note that this isn't a pointer.
  // Order matters. The timers will be printed before other
  // data members are destructed.
  BulkAdHocTimers _timers;

public:
  ZgyInternalBulk(
      const std::shared_ptr<IFileADT>& file,
      const std::shared_ptr<const ZgyInternalMeta>& metadata,
      const std::shared_ptr<ZgyInternalMeta>& metadata_rw,
      bool compressed_write,
      const std::array<float,2>& force_histo_range,
      const LoggerFn& logger = LoggerFn());

  std::pair<bool,double> readConstantValue(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod, bool as_float) const;

  void readToExistingBuffer(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      int32_t lod, bool as_float) const;

  void readBricksToExistingBuffer(
     const std::vector<std::shared_ptr<DataBuffer>>& result,
     const std::vector<std::array<std::int64_t,3>>& start,
     int32_t lod, bool as_float) const;

  bool expeditedRead(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      void* data, int lod, RawDataType result_type) const;

  void writeRegion(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      int32_t lod, bool is_storage,
      const compressor_t& compressor);

  std::vector<std::shared_ptr<DataBuffer>>
  readBricksToNewBuffers(
      const std::vector<std::array<std::int64_t, 3>>& position_in,
      const std::vector<std::shared_ptr<DataBuffer>>& data_in,
      int lod, bool as_float, bool check_constant) const;

  void readBricksInternal(
       const std::vector<std::array<std::int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod, bool as_float) const;

  void writeBricksInternal(
       const std::vector<std::array<std::int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod, bool is_storage, const compressor_t& compressor,
       bool immutable_buffer, bool maybe_more);

  void writeConstBricksInternal(
       const std::vector<std::array<int64_t,3>>& position,
       const std::shared_ptr<InternalZGY::DataBuffer>& value,
       bool is_storage);

public: // actually internal

  bool errorflag() const        { return _is_bad.load(); }
  void set_errorflag(bool flag) { _is_bad.store(flag); }
  LoggerFn set_logger(const LoggerFn& logger) {
    LoggerFn old = _loggerfn;
    _loggerfn = logger;
    return old;
  }
  bool isStatHistGood() const   { return _new_stathist_good; }
  bool wasLowResRequested() const { return _new_internal_lod_mode.last.force >= 2 && _new_internal_lod_mode.last.level < 0; }

private:
  std::vector<std::pair<bool, double>>
  readConstValueBricks(
      const std::vector<std::array<std::int64_t, 3>>& position,
      int32_t lod, bool as_float) const;
  static bool coversEntireSurvey(
       const std::vector<std::array<std::int64_t,3>>& position,
       const std::array<std::int64_t,3>& surveysize,
       const std::array<std::int64_t,3>& bricksize);
  static std::pair<bool, double> isSameScalarValue(
       const std::vector<std::shared_ptr<DataBuffer>>& data);
  static std::pair<bool, double> isResettingEntireSurvey(
       const std::vector<std::array<std::int64_t,3>>& position,
       const std::array<std::int64_t,3>& surveysize,
       const std::array<std::int64_t,3>& bricksize,
       const std::vector<std::shared_ptr<DataBuffer>>& data);
  void newTrackedBricksSetAllConst(double value);
  void newTrackedBricksInitFromFile();
  void newTrackedBricksTryEnable();

public:
  std::tuple<std::shared_ptr<StatisticData>, std::shared_ptr<HistogramData>>
  newTrackedChanges() const {
    return std::make_tuple(_new_modified_stats, _new_modified_histo);
  }
  std::shared_ptr<TrackTouched> newTrackedBricks() const {
    return _new_modified_bricks;
  }
  void newSetGenLodInstance(
       const std::shared_ptr<GenLodC>& genlod,
       const InternalLodMode& internal_lod_mode)
  {
    _new_genlod = genlod;
    _new_internal_lod_mode = internal_lod_mode;
  }
  void newFinalize(const std::function<bool(std::int64_t,std::int64_t)>& p);

private:

  bool blogger(int priority, const std::string& ss = std::string()) const;

  bool blogger(int priority, const std::ios& ss) const;

  void validateUserLod(int32_t lod) const;

  void validateUserPosition(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const;

  bool singleBrickOutsideSurvey(
    const std::array<std::int64_t,3>& start,
    const std::array<std::int64_t,3>& size,
    int32_t lod) const;

  std::shared_ptr<DataBuffer> scaleDataToFloat(
      const std::shared_ptr<DataBuffer>& data) const;

  std::shared_ptr<DataBuffer> scaleDataToStorage(
      const std::shared_ptr<DataBuffer>& data) const;

  static double decodeConstant(std::uint32_t in, RawDataType dtype);

  static std::int32_t encodeConstant(double in, RawDataType dtype);

  double missingValue(bool as_float) const;

  LutInfoEx makeLutInfoEx(
      const std::array<std::int64_t,3>& start,
      int32_t lod) const;

  std::vector<LutInfoEx> partsNeeded(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const;

  void deliverOneBrick(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& startpos,
      const ReadRequest::data_t& raw, std::int64_t rawsize,
      BrickStatus brickstatus, bool as_float) const;

  ReadRequest::delivery_t createDeliverance(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      bool as_float,
      const LutInfoEx& brick) const;

  // --- WRITE SUPPORT ---

  std::array<std::int64_t,3> usedPartOfBrick(
      const std::array<std::int64_t,3>& size,
      const std::array<std::int64_t,3>& brickpos,
      std::int32_t lod) const;

  bool isUsedPartOfBrickAllConstant(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& brickpos,
      int32_t lod) const;

  static void setPaddingToEdge(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      int modulo, int dim);

  static void setPaddingToConst(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, int dim);

  static void setPaddingSamples(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, const compressor_t& compressor);

  void writeOneNormalBrickWithRetry(const WriteNowArgPack& args);

  std::shared_ptr<const WriteNowArgPack>
  prepareWriteOneNormalBrick(const WriteBrickArgPack& args);

  void writeOneConstantBrick(const WriteBrickArgPack& args);

  std::int64_t writeAllConstantBrick(double value, bool force);

  bool mustLeakOldBrick(
      const std::shared_ptr<DataBuffer>& data,
      const compressor_t& compressor,
      BrickStatus brickstatus,
      std::array<std::int64_t,3>& pos,
      std::int32_t lod) const;

  std::shared_ptr<const WriteBrickArgPack>
  prepareWriteOneBrick(const WriteBrickArgPack& args);

#if 0
  static void
  addCompressionNoise(
       std::shared_ptr<DataBuffer> data,
       compressor_t compressor);
#endif

  static std::pair<index3_t,index3_t>
  clip(
       const std::pair<index3_t,index3_t>& a,
       const std::pair<index3_t,index3_t>& b);

  static std::shared_ptr<DataBuffer>
  doReadModifyWrite(
       std::shared_ptr<DataBuffer> old_brick,
       std::shared_ptr<DataBuffer> new_brick,
       const std::pair<index3_t,index3_t>& brick_range,
       const std::pair<index3_t,index3_t>& inside_survey,
       const std::pair<index3_t,index3_t>& inside_user,
       double defaultstorage,
       const compressor_t& compressor,
       bool immutable_buffer);

  void
  writeAlignedBrickList(
       const std::vector<std::shared_ptr<DataBuffer>>& new_contents,
       const std::vector<index3_t>& position,
       const std::pair<index3_t,index3_t>& user_clip,
       int lod, bool is_storage, const compressor_t& compressor,
       bool immutable_buffer);

  void
  writeAllConstantSurvey(double value);

  static void
  splitRegionIntoBricks(
       std::vector<std::shared_ptr<DataBuffer>>* bricks /*out*/,
       std::vector<index3_t>* positions /*out*/,
       const std::shared_ptr<DataBuffer>& data,
       const std::array<std::int64_t,3>& start,
       const index3_t& bs);

  bool
  checkBrickInternal(
       const std::array<int64_t,3>& position,
       const std::shared_ptr<DataBuffer>& data,
       bool as_float, bool strict) const;

  void
  checkBricksInternal(
       const std::vector<std::array<int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       bool as_float, bool strict) const;
};

} // namespace
