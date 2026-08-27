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
class TrackTouched;

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
 *
 * SINGLEPASS: Data members and methods removed from GenLodBase and derived.
 * The list might not be complete. The names might linger on in comments.
 *
 *   _stats _histo _incremental _skip_histogram
 *
 *   call() callNewVersion() callSingleBrick() suggestHistogramRange()
 *   _read() _write() _willneed() _isclean() _use_plan_a()
 *   _clipsizetosurvey() _accumulateT() _accumulate()
 *   _accumulateAndStore() accumulateOneInputBrick() _calculate()
 *   _decimate() _paste1() _paste4()
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
  std::vector<LodAlgorithm> _decimation;
  std::shared_ptr<HistogramData> _wa_histogram;
  double _wa_defaultstorage;
  std::function<bool(std::int64_t,std::int64_t)> _progress;
  LoggerFn _loggerfn;
  SummaryPrintingTimerEx _timer_histogram;

  virtual ~GenLodBase();
  // Technically copyable, but copy/assign doesn't make much sense.
  GenLodBase(const GenLodBase&) = delete;
  GenLodBase(const GenLodBase&&) = delete;
  GenLodBase& operator=(const GenLodBase&) = delete;
  GenLodBase& operator=(const GenLodBase&&) = delete;

public:
  GenLodBase(const index3_t& size,
             const index3_t& bricksize,
             RawDataType dtype,
             std::int32_t nlods_in,
             const std::vector<LodAlgorithm>& decimation,
             const std::shared_ptr<HistogramData>& histogram,
             double defaultvalue,
             const std::function<bool(std::int64_t,std::int64_t)>& progress,
             LoggerFn logger);
protected:
  virtual std::shared_ptr<TrackTouched> _getNewTouched() const;
  virtual std::vector<std::shared_ptr<DataBuffer>>
  _readbricks(
       const std::vector<std::array<int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod) const;
  virtual void
  _writebricks(
       const std::vector<std::array<int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod) const;

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
protected:
  virtual ~GenLodImpl();

public:
  GenLodImpl(const index3_t& size,
             const index3_t& bricksize,
             RawDataType dtype,
             std::int32_t nlods_in,
             const std::vector<LodAlgorithm>& decimation,
             const std::shared_ptr<HistogramData>& wa_histogram,
             double defaultvalue,
             const std::function<bool(std::int64_t,std::int64_t)>& progress,
             LoggerFn logger);
  void tryToCall(
       const LevelAndForce& lodmode,
       const std::vector<index3_t>& brickroi,
       const std::function<bool(std::int64_t,std::int64_t)>& progress);
protected:
  static bool areAllEightSameConst(
       const std::vector<std::shared_ptr<DataBuffer>>& indata,
       std::int64_t ipos);
  void processFourInputColumnsOneLod(
       const int out_lod,
       const std::int64_t pos0,
       const std::int64_t pos1,
       const std::int64_t zcount,
       const LodAlgorithm algorithm,
       const std::vector<std::shared_ptr<DataBuffer>>& all_buffers);
  void decimateOneInputBrick(
       const std::vector<std::shared_ptr<DataBuffer>>& outbuffer,
       const std::vector<std::shared_ptr<DataBuffer>>& indata,
       std::int64_t part,
       LodAlgorithm algorithm);
};

/**
 * Implement plan "C". See documentation.
 *
 * Thread safety: Instances are not thread safe. See GenLodBase.
 *
 * \internal
 *
 * A full build of low resolution bricks, histogram, and statistics is
 * just a single call to GenLodC::call(). A lot of parameters are
 * required. It is impractical to ripple them down to private methods
 * that need them. Instead, a short lived instance is used to hold the
 * information. Instead of the application calling
 * GenLodC().call(many_arguments), it calls GenLodC(many_arguments),call().
 *
 * GenLod instances hold references to:
 *   - Mutable StatisticData and HistogramData to be returned.
 *   - Mutable done and total counts.
 *   - HistogramData used by the WeightedAverage decimation.
 *   - Progress and logger callback.
 *   - Bulk pointer and compressor callback.
 *   - Several other POD data.
 *
 * genlod does not keep an explicit reference to the metatdata passed to
 * the constructor; it is only used during construction. There is an
 * indirect reference from the bulk instance.
 *
 * For incremental builds, the bulk instance will own a long lived genlod
 * instance and will do multiple calls to the tryToCall() method. This
 * bends the original design rather out of shape. E.g. there will be a
 * circular reference between the genlod instance and the bulk instance.
 *
 * It is technically possible to use multiple short-lived genlod
 * instances instead of one long lived one. The mutable stats and
 * histo can be owned by the caller. The mutable done and total are
 * only for progress reports, which might not make sense during an
 * incremental build.
 *
 * TODO-WIP-BrickedAPI: Better fix for cyclical reference problem.
 *
 * A bulk instance now has a reference to a long lived genlod instance,
 * allowing it to write lod bricks incrementaly. This genlod instance
 * has a reference back to the bulk instance to query for work and to
 * write lowres data to file. There is no infinite recursion, because
 * the bulk->genlod call only happens for lod 0, and the genlod->bulk
 * only happens for lod>0. But the cyclical reference prevent the two
 * instances from ever being destroyed. And, much worse, prevents the
 * file from being properly flushed if the application is missing
 * an explicit close.
 *
 * - [current] Make the genlod->bulk reference a weak_ptr. This is the
 *   simplest fix. It is not exactly pretty.
 *
 * - Refactor tryToCall() to have the bulk accessor instance passed as
 *   a parameter instead of keeping a reference. This requires either
 *   a massive kludge (store the pointer in function entry, release it
 *   on exit) or a large restructuring of the code.
 *
 * - Use a short lived genlod instance also on each incremental write.
 *   This is closer to the original design. But there is a small cost
 *   in creating and destroying the genlod instance. Also, only the
 *   API layer can create instances. So, add a functor that the api
 *   layer can give the bulk layer to request an instance. And
 *   suddenly, it isn't that elegant any longer.
 *   Not only that, but the functor could easily end up with its
 *   own reference to the api layer - which would be even worse.
 *
 * Related: Could the ZgyWriter::_accessor be a unique_ptr instead of
 * shared_ptr? This helps ensure that references aren't leaked.
 * Unfortunately the ZombieCheck etc. would no longer work. And many
 * function signatures need to change from shared_ptr to plain
 * references or pointers.
 */
class GenLodC : public GenLodImpl
{
  std::weak_ptr<ZgyInternalBulk> _accessor;
  compressor_t _lodcompressor;
public:
  GenLodC(const std::shared_ptr<ZgyInternalBulk>& accessor,
          const std::shared_ptr<ZgyInternalMeta>& meta,
          const compressor_t& lodcompressor,
          const std::vector<LodAlgorithm>& decimation,
          const std::shared_ptr<InternalZGY::HistogramData>& wa_histogram,
          LoggerFn logger);
  virtual ~GenLodC();

protected:
  std::shared_ptr<TrackTouched> _getNewTouched() const override;

  std::vector<std::shared_ptr<DataBuffer>>
  _readbricks(
       const std::vector<std::array<int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod) const override;

  void
  _writebricks(
       const std::vector<std::array<int64_t,3>>& position,
       const std::vector<std::shared_ptr<DataBuffer>>& data,
       int lod) const override;

};

} // namespace
