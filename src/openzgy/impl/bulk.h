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

#include "declspec.h"
#include "enum.h"
#include "file.h" // Only for the ReadRequest::data_t typedef

#include <memory>
#include <functional>
#include <array>
#include <vector>
#include <tuple>
#include <sstream>
#include <atomic>

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
  std::shared_ptr<ZgyInternalMeta> _get_metadata_rw();
  UpdateMode _update_mode;
  bool       _compressed_write; // If true: do not align bricks.
  std::atomic<bool> _is_bad;    // If true: instance not usable.
  double     _written_sample_min;
  double     _written_sample_max;
  LoggerFn   _loggerfn;
  std::shared_ptr<SummaryPrintingTimerEx> _ptimer_st;
  std::shared_ptr<SummaryPrintingTimerEx> _ptimer_mt;
  std::shared_ptr<SummaryPrintingTimerEx> _ststimer;
  std::vector<std::uint8_t>         _modified_bricks;
  std::shared_ptr<StatisticData>    _modified_stats;
  std::shared_ptr<HistogramData>    _modified_histo;

public:
  ZgyInternalBulk(
      const std::shared_ptr<IFileADT>& file,
      const std::shared_ptr<const ZgyInternalMeta>& metadata,
      const std::shared_ptr<ZgyInternalMeta>& metadata_rw,
      bool compressed_write,
      const LoggerFn& logger = LoggerFn());

  std::pair<bool,double> readConstantValue(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod, bool as_float) const;

  void readToExistingBuffer(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      int32_t lod, bool as_float) const;

  std::shared_ptr<DataBuffer> readToNewBuffer(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod, bool as_float, bool check_constant) const;

  bool expeditedRead(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      void* data, int lod, RawDataType result_type) const;

  void writeRegion(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      int32_t lod, bool is_storage,
      const compressor_t& compressor);

public: // actually internal

  bool errorflag() const        { return _is_bad.load(); }
  void set_errorflag(bool flag) { _is_bad.store(flag); }
  LoggerFn set_logger(const LoggerFn& logger) {
    LoggerFn old = _loggerfn;
    _loggerfn = logger;
    return old;
  }
  std::array<double,2> valueRangeWritten() const {
    return std::array<double,2>{_written_sample_min, _written_sample_max};
  }
  std::tuple<std::shared_ptr<StatisticData>, std::shared_ptr<HistogramData>> trackedChanges() const {
    return std::make_tuple(_modified_stats, _modified_histo);
  }
  const std::vector<std::uint8_t>& trackedBricks() const {
    return _modified_bricks;
  }
  void trackedBricksTryEnable(bool on);
  std::uint8_t trackedBricksDirty(
       const std::array<std::int64_t,3>& start,
       const std::array<std::int64_t,3>& size,
       int32_t lod) const;
  void trackedBricksShowDirty(int loglevel) const;

private:

  void _trackedBricksSetDirty(
     const std::vector<index3_t>& work,
     int32_t lod);

  bool _logger(int priority, const std::string& ss = std::string()) const;

  bool _logger(int priority, const std::ios& ss) const;

  void _validateUserPosition(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const;

  std::shared_ptr<DataBuffer> _scaleDataToFloat(
      const std::shared_ptr<DataBuffer>& data) const;

  std::shared_ptr<DataBuffer> _scaleDataToStorage(
      const std::shared_ptr<DataBuffer>& data) const;

  static double _decodeConstant(std::uint32_t in, RawDataType dtype);

  static std::int32_t _encodeConstant(double in, RawDataType dtype);

  std::vector<LutInfoEx> _partsNeeded(
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& size,
      int32_t lod) const;

  void _deliverOneBrick(
      const std::shared_ptr<DataBuffer>& result,
      const std::array<std::int64_t,3>& start,
      const std::array<std::int64_t,3>& startpos,
      const ReadRequest::data_t& raw, std::int64_t rawsize,
      BrickStatus brickstatus, bool as_float) const;

  // --- WRITE SUPPORT ---

  std::array<std::int64_t,3> _usedPartOfBrick(
      const std::array<std::int64_t,3>& size,
      const std::array<std::int64_t,3>& brickpos,
      std::int32_t lod) const;

  bool _isUsedPartOfBrickAllConstant(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& brickpos,
      int32_t lod) const;

  static void _setPaddingToEdge(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      int modulo, int dim);

  static void _setPaddingToConst(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, int dim);

  static void _setPaddingSamples(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& used,
      double missingvalue, const compressor_t& compressor);

  void _writeWithRetry(const WriteNowArgPack& args);

  std::shared_ptr<const WriteNowArgPack>
  _writeOneNormalBrick(const WriteBrickArgPack& args);

  void _writeOneConstantBrick(const WriteBrickArgPack& args);

  bool _mustLeakOldBrick(
      const std::shared_ptr<DataBuffer>& data,
      const compressor_t& compressor,
      BrickStatus brickstatus) const;

  std::shared_ptr<const WriteBrickArgPack>
  _writeOneBrick(const WriteBrickArgPack& args);

  void _writeAlignedRegion(
      const std::shared_ptr<DataBuffer>& data,
      const std::array<std::int64_t,3>& start,
      std::int32_t lod,
      const compressor_t& compressor);
};

} // namespace
