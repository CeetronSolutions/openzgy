// Copyright 2017-2020, Schlumberger
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

#include "enum.h"
#include "lodalgo.h"
#include "histogramdata.h"
#include "histogrambuilder.h"
#include "statisticdata.h"
#include "databuffer.h"
#include "../exception.h"
#include "fancy_timers.h"

#include <memory>
#include <cstdint>
#include <array>
#include <vector>
#include <functional>
#include <iostream>

/**
 * \file genlod.h
 * \brief Generate low resolution bricks.
 */

namespace InternalZGY {
#if 0
}
#endif

class ZgyInternalMeta;
class ZgyInternalBulk;

/**
 * Abstract class for generating low resolution bricks, histogram,
 * and statistics. At this level only define virtual methods for I/O.
 * The implementation can be used as-is when mocking the class.
 * The optional nlods parameter is only used as a consistency check.
 *
 * Thread safety: Instances are not thread safe. Nor need they be.
 * An instance of this class represents a running low-res computation
 * invoked from finalize() and it makes no sense to have more than one
 * of those. This applies to derived classes as well.
 */
class GenLodBase
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

protected:
  std::int32_t _nlods;
  std::int64_t _total;
  std::int64_t _done;
  index3_t _surveysize;
  index3_t _bricksize;
  RawDataType _dtype;
  std::array<double,2> _histogram_range;
  std::vector<LodAlgorithm> _decimation;
  std::shared_ptr<HistogramData> _wa_histogram;
  double _wa_defaultstorage;
  bool _incremental;
  bool _skip_histogram;
  std::function<bool(std::int64_t,std::int64_t)> _progress;
  LoggerFn _loggerfn;
  SummaryPrintingTimerEx _timer_histogram;

public:
  GenLodBase(const index3_t& size,
             const index3_t& bricksize,
             RawDataType dtype,
             const std::array<double,2>& histogram_range,
             std::int32_t nlods_in,
             const std::vector<LodAlgorithm>& decimation,
             const std::shared_ptr<HistogramData>& histogram,
             double defaultvalue,
             bool incremental,
             bool skip_histogram,
             const std::function<bool(std::int64_t,std::int64_t)>& progress,
             LoggerFn logger);
protected:
  virtual std::int64_t  _willneed() const;
  virtual bool _isclean(
      std::int32_t lod, const index3_t& pos, const index3_t& size) const;
  virtual std::shared_ptr<DataBuffer> _read(
      std::int32_t lod, const index3_t& pos, const index3_t& size) const;
  virtual void _write(
      std::int32_t lod, const index3_t& pos,
      const std::shared_ptr<const DataBuffer>& data) const;
  index3_t _clipsizetosurvey(std::int32_t lod, const index3_t& pos, const index3_t& size) const;
  void _report(const DataBuffer* data) const;
  void _reporttotal(std::int64_t total) { _total = total; }
  std::string _prefix(std::int32_t lod) const;
  static std::string _format_result(const std::shared_ptr<DataBuffer>& data);
  bool _logger(int priority, const std::string& message) const;
  bool _logger(int priority, const std::ios& ss) const;
};

/**
 * Abstract class for generating low resolution bricks, histogram,
 * and statistics. The inherited methods for I/O are still stubs.
 *
 * Thread safety: Instances are not thread safe. See GenLodBase.
 */
class GenLodImpl: public GenLodBase
{
  std::shared_ptr<StatisticData> _stats;
  std::shared_ptr<HistogramData> _histo;

public:
  GenLodImpl(const index3_t& size,
             const index3_t& bricksize,
             RawDataType dtype,
             const std::array<double,2>& histogram_range,
             std::int32_t nlods_in,
             const std::vector<LodAlgorithm>& decimation,
             const std::shared_ptr<HistogramData>& histogram,
             double defaultvalue,
             bool incremental,
             bool skip_histogram,
             const std::function<bool(std::int64_t,std::int64_t)>& progress,
             LoggerFn logger);
  std::tuple<std::shared_ptr<StatisticData>, std::shared_ptr<HistogramData>>
  call();
protected:
  template <typename T>
  void _accumulateT(const std::shared_ptr<const DataBuffer>& data_in);
  void _accumulate(const std::shared_ptr<const DataBuffer>& data);
  std::shared_ptr<DataBuffer>
  _calculate(const index3_t& readpos, const index3_t& readsize, std::int32_t readlod);
  std::shared_ptr<DataBuffer>
  _decimate(const std::shared_ptr<const DataBuffer>& data, std::int64_t lod);
  std::shared_ptr<DataBuffer>
  _paste1(const std::shared_ptr<DataBuffer>& result,
          const std::shared_ptr<const DataBuffer>& more,
          std::int64_t ioff, std::int64_t joff);
  std::shared_ptr<const DataBuffer>
  _paste4(const std::shared_ptr<const DataBuffer>& d00,
          const std::shared_ptr<const DataBuffer>& d01,
          const std::shared_ptr<const DataBuffer>& d10,
          const std::shared_ptr<const DataBuffer>& d11);
  static std::array<double,2>
  suggestHistogramRange(
    const std::array<double,2>& writtenrange,
    RawDataType dtype);
};

/**
 * Implement plan "C". See documentation.
 *
 * Thread safety: Instances are not thread safe. See GenLodBase.
 */
class GenLodC : public GenLodImpl
{
  std::shared_ptr<ZgyInternalBulk> _accessor;
  compressor_t _lodcompressor;
public:
  GenLodC(const std::shared_ptr<ZgyInternalBulk>& accessor,
          const std::shared_ptr<ZgyInternalMeta>& meta,
          const compressor_t& lodcompressor,
          const std::vector<LodAlgorithm>& decimation,
          bool incremental,
          bool skip_histogram,
          const std::function<bool(std::int64_t,std::int64_t)>& progress,
          LoggerFn logger);
protected:
  std::int64_t  _willneed() const override;
  bool _isclean(
      std::int32_t lod, const index3_t& pos, const index3_t& size) const override;
  std::shared_ptr<DataBuffer> _read(
      std::int32_t lod, const index3_t& pos, const index3_t& size) const override;
  void _write(
      std::int32_t lod, const index3_t& pos,
      const std::shared_ptr<const DataBuffer>& data) const override;
};

} // namespace
