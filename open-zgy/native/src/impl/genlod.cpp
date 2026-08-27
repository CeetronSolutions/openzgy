// Copyright 2017-2023, Schlumberger
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

#include "edgebrick.h"

#include "genlod.h"
#include "enum.h"
#include "types.h"
#include "lodalgo.h"
#include "histogramdata.h"
#include "histogrambuilder.h"
#include "statisticdata.h"
#include "databuffer.h"
#include "arrayops.h"
#include "meta.h"
#include "bulk.h"
#include "../exception.h"
#include "wintools.h"
#include "workorder.h"
#include "environment.h"
#include "tracktouched.h"

#include <memory>
#include <cstdint>
#include <array>
#include <vector>
#include <functional>
#include <sstream>
#include <tuple>

namespace InternalZGY {
#if 0
}
#endif

using namespace ArrayOps;
using namespace Formatters;

/**
 * When timers are active, might want more verbose logging,
 * but not so much that it warrants setting OPENZGY_VERBOSE.
 */
static bool openzgy_timers()
{
  static int result = Environment::getNumericEnv("OPENZGY_TIMERS", 0);
  return result > 0;
}

//# TODO-Low: Several tweaks are possible but might not offer
//# much value. See a detailed discussion in doc/lowres.html.

//# Testing notes:
//#   test.black.checkLodContents(), checkStatistics(), checkHistogram()
//#   is part of a reasonably complete end to end test checking this module.
//#   A new unit test testFinalizeProgress() has been written to verify the
//#   progress mechanism. Both that the done and total end up precisely
//#   identical at the end and that the client can abort the calculation.

/**
 * Abstract class for generating low resolution bricks, histogram,
 * and statistics. At this level only define virtual methods for I/O.
 * The implementation can be used as-is when mocking the class.
 * The optional nlods parameter is only used as a consistency check.
 *
 * Note that the WeightedAverage algorithm requires a histogram to work.
 * If no histogram was provided then the current contents of the
 * accumulated histogram will be used. This is unfortunate and might
 * cause brick artifacts. Especially in the first few bricks that are
 * generated. With a non-recursive algorithm (plan B) and with only
 * lod2 and above uses weighted average then this is unproblematic.
 * Because in that case we will be done with the histogram once we
 * need it. TODO-Low consider doing an initial pass with a statistical
 * sampling of the lod0 data, only for use with weighted average.
 * There will be a minor issue with some values appearing to have zero
 * frequency, but this should not cause any trouble. (assume "1").
 *
 * Note that the WeightedAverage and AverageNon0 algorithms expect a
 * defaultstorage to use when all inputs are inf/nan or (for AverageNon0)
 * zero. Only relevant for integral types, to ensure that the default
 * is whatever will produce the value closest to 0 after conversion.
 * And integral data can neither be inf nor nan, so this is a pretty
 * academic issue. For AverageNon0 that algorithm is currently not
 * used. So it isn't even clear what the desired behavior should be.
 */
GenLodBase::GenLodBase(
    const std::array<std::int64_t,3>& size,
    const std::array<std::int64_t,3>& bricksize,
    RawDataType dtype,
    std::int32_t nlods_in,
    const std::vector<LodAlgorithm>& decimation,
    const std::shared_ptr<HistogramData>& histogram,
    double defaultstorage,
    const std::function<bool(std::int64_t,std::int64_t)>& progress,
    LoggerFn logger)
  : _nlods(nlods_in)
  , _total(0)
  , _done(0)
  , _surveysize(size)
  , _bricksize(bricksize)
  , _dtype(dtype)
  , _decimation(decimation.empty() ? std::vector<LodAlgorithm>{
      LodAlgorithm::LowPass, LodAlgorithm::WeightedAverage} : decimation)
  , _wa_histogram(histogram)
  , _wa_defaultstorage(defaultstorage)
  , _progress(progress)
  , _loggerfn(logger)
  , _timer_histogram("GenLod.hist")
{
  // Do a sanity check on the supplied nlods.
  std::int32_t nlods = 1; // The loop stops before counting final level
  std::array<std::int64_t,3> bs = bricksize;
  std::array<std::int64_t,3> sz = size;
  while (sz[0]>bs[0] || sz[1] > bs[1] || sz[2] > bs[2]) {
    nlods += 1;
    sz = (sz + std::int64_t(1)) / std::int64_t(2);
  }
  if (nlods_in > 0 && nlods_in != nlods)
    throw OpenZGY::Errors::ZgyInternalError("GenLod error in nlods computation");
}

GenLodBase::~GenLodBase()
{
}

/**
 * \brief Allow querying for which bricks are ready for low resolution.
 */
std::shared_ptr<TrackTouched>
GenLodBase::_getNewTouched() const
{
  return nullptr;
}

/**
 * This is a stub that can be overridden later. If not overridden,
 * and code attempts to call it, then an exception is raised,
 *
 * Read one or more bricks from the ZGY file (plans A, B and C).
 * The caller can but doesn't need to pass a vector of buffers
 * to be used. And even when it does, the function isn't obliged
 * to use them. E.g. because some bricks were all-const and should
 * be returned as a scalar buffer.
 */

std::vector<std::shared_ptr<DataBuffer>>
GenLodBase::_readbricks(
     const std::vector<std::array<int64_t,3>>& position,
     const std::vector<std::shared_ptr<DataBuffer>>& buffers,
     int lod) const
{
  throw OpenZGY::Errors::ZgyInternalError("GenLodBase::_readbricks() missing");
}

void
GenLodBase::_writebricks(
     const std::vector<std::array<int64_t,3>>& position,
     const std::vector<std::shared_ptr<DataBuffer>>& data,
     int lod) const
{
  throw OpenZGY::Errors::ZgyInternalError("GenLodBase::_writebricks() missing");
}

/**
 * For debugging and logging only.
 */
std::string
GenLodBase::_prefix(std::int32_t lod) const
{
  return std::string(2 * (this->_nlods - 1 - lod), ' ');
}

/**
 * For debugging and logging only.
 */
std::string
GenLodBase::_format_result(const std::shared_ptr<DataBuffer>& data)
{
  return data ? data->toString() : "(null)";
}

/**
 * Convenience for invoking _loggerfn with a simple message.
 * This isn't much different from invoking the callback directly.
 * But it makes debugging slightly simpler to have an easy place
 * to set a breakpoint. It also adds more symmetry with respect
 * to the stringstream version, which does add value.
 */
bool
GenLodBase::_logger(int priority, const std::string& message) const
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
GenLodBase::_logger(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return _logger(priority, sstream ? sstream->str() : std::string());
}

/**
 * Abstract class for generating low resolution bricks, histogram,
 * and statistics. The inherited methods for I/O are still stubs.
 * See doc/lowres.html for details. This class implements plan C or D
 * which is good for compressed data and acceptable for uncompressed.
 * The ordering of low resolution bricks in the file will not be optimal.
 * For optimal ordering but working only for uncompressed data consider
 * implementing plan B in addition to the plan C already implemented.
 * The implementation can be used as-is in a unit test with mocked I/O.
 */
GenLodImpl::GenLodImpl(
    const std::array<std::int64_t,3>& size,
    const std::array<std::int64_t,3>& bricksize,
    RawDataType dtype,
    std::int32_t nlods_in,
    const std::vector<LodAlgorithm>& decimation,
    const std::shared_ptr<HistogramData>& wa_histogram,
    double defaultstorage,
    const std::function<bool(std::int64_t,std::int64_t)>& progress,
    LoggerFn logger)
  : GenLodBase(size, bricksize, dtype, nlods_in, decimation,
               wa_histogram, defaultstorage, progress, logger)
{
}

GenLodImpl::~GenLodImpl()
{
}

/**
 * The caller has a multiple of 8 input data bricks to be decimated,
 * with each group of 8 ordered dim2 fastest, dim0 slowest. I.e. same
 * ordering of bricks as ordering of samples inside each brick.
 * And one output brick per group.
 *
 * Return the offset into the decimated output brick, in samples,
 * corresponding to the first sample in inputs[index].
 */
static std::int64_t
offsetInOutput(std::int64_t index, const index3_t& bs)
{
  //const std::int64_t vbrick = index / 8;
  const std::int64_t subi = (index / 4) % 2;
  const std::int64_t subj = (index / 2) % 2;
  const std::int64_t subk = index % 2;
  return (subk * bs[2] / 2 +
          subj * bs[2] * bs[1] / 2 +
          subi * bs[2] * bs[1] * bs[0] / 2);
}

/**
 * Decimate a single input cube, storing the result in 1/8 of an output
 * cube. The index into outdata will always be input_buffer_num/8.
 * index mod 8 is used to choose the correct octant in the output,
 * so it cannot be chosen freely.
 */
void
GenLodImpl::decimateOneInputBrick(
     const std::vector<std::shared_ptr<DataBuffer>>& outbuffer,
     const std::vector<std::shared_ptr<DataBuffer>>& indata,
     std::int64_t input_buffer_num,
     LodAlgorithm algorithm)
{
  if (input_buffer_num < 0 ||
      (std::size_t)input_buffer_num >= indata.size() ||
      (std::size_t)input_buffer_num / 8 >= outbuffer.size())
  {
    throw OpenZGY::Errors::ZgyInternalError("Buffer overrun prevented");
  }
  if (!outbuffer[input_buffer_num / 8])
    throw OpenZGY::Errors::ZgyInternalError("Output decimation buffer is null");
  if (outbuffer[input_buffer_num / 8]->isScalar())
    throw OpenZGY::Errors::ZgyInternalError("Output decimation buffer is scalar");
  const index3_t bs = _bricksize;
  const std::int64_t vbrick = input_buffer_num / 8;
  const std::int64_t output_offset = offsetInOutput(input_buffer_num, bs);
  // TODO-Medium missing _wa_defaultstorage.
  if (_logger(3, ""))
    _logger(3, std::stringstream()
      << "    createLod buffer " << input_buffer_num
      << " offset " << output_offset);
  createLod(outbuffer[vbrick], indata[input_buffer_num], algorithm,
    output_offset,
    _wa_histogram->getbins(),
    _wa_histogram->getsize(),
    _wa_histogram->getmin(),
    _wa_histogram->getmax());
}

/**
 * See also ZgyInternalBulk::readConstantValue().
 * Much of the logic is the same. And not entirely trivial.
 */
bool /*static*/
GenLodImpl::areAllEightSameConst(
     const std::vector<std::shared_ptr<DataBuffer>>& indata,
     std::int64_t ipos)
{
  if (indata.size() < (std::size_t)ipos + 8)
    throw OpenZGY::Errors::ZgyInternalError("Buffer overrun prevented");
  for (std::int64_t ix = ipos; ix < ipos + 8; ++ix)
    if (!indata[ix]->isScalar())
      return false;
  const double value = indata[ipos]->scalarAsDouble();
  if (std::isnan(value)) {
    for (std::int64_t ix = ipos + 1; ix < ipos + 8; ++ix)
      if (!std::isnan(indata[ix]->scalarAsDouble()))
        return false;
  }
  else {
    for (std::int64_t ix = ipos + 1; ix < ipos + 8; ++ix)
      if (indata[ix]->scalarAsDouble() != value)
        return false;
  }
  return true;
}

namespace {
/**
 * Helper for invoking the user's callback function, if any.
 * Makes sure the rules for progress reporting are adhered to:
 *
 *   - need > 0 even when there isn't actually any work to be done.
 *   - done == 0 as early as possible to tell the app when "need" is.
 *   - done < need in every call but the last one, which has done == need.
 *
 * If created with throw_on_cancel = false, calls to start(), step(),
 * finish() all need to check the result and stop processing data.
 */
class SendProgress
{
public:
  typedef std::function<bool(std::int64_t,std::int64_t)> progress_t;

private:
  progress_t   progress_;
  bool         throw_;
  std::int64_t done_;
  std::int64_t need_;
  bool         finished_;

public:
  SendProgress(progress_t progress, bool throw_on_cancel = true)
    : progress_(progress)
    , throw_(throw_on_cancel)
    , done_(0)
    , need_(0)
    , finished_(false)
  {
  }

  bool start(const TrackTouched::tasklist_t& tasklist)
  {
    done_ = need_ = 0;
    for(const TrackTouched::task_t& task : tasklist)
      need_ += TrackTouched::countIO(task);
    if (need_ < 1)
      need_ = 1;
    return report(0, need_);
  }

  bool step(const TrackTouched::task_t& task)
  {
    done_ += TrackTouched::countIO(task);
    return report(std::min(done_, need_-1), need_);
  }

  bool finish()
  {
    // The next two are errors, but no big deal to silently correct them.
    //if (done_ != need_ && need_ != 1) {
    //  logger(0, std::stringstream().flush()
    //         << "Internal error, progress report is wrong: "
    //         << done << "/" << need << " in last call.");
    //}
    //if (finished_)
    //  logger(0, "ERROR: SendProgress::finish() called more than once.");
    if (!finished_) {
      finished_ = true;
      return report(need_, need_);
    }
    return true;
  }

private:
  bool report(std::int64_t done, std::int64_t need)
  {
    bool ok = (!progress_ || progress_(done, need_));
    if (!ok && throw_)
      throw OpenZGY::Errors::ZgyAborted
        ("Computation of low resolution data was aborted");
    return ok;
  }
};
}

/*
 * Notes from callNewVersion, TODO-WIP-BrickedAPI: move to some other function.
 * EXPERIMENTAL!
 *
 * TODO-WIP-BrickedAPI: Possible shortcut.
 * When processing 8 bricks in, 1 brick out, check whether each of the
 * inputs are scalar and all have the same value. If so, just copy the
 * first input brick to the output brick. Making the output scalar as
 * well. And do the short cut handling for histogram and statistics.
 * This short cut might be important for incremental LOD generation.
 *
 * TODO-WIP-BrickedAPI: Possible shortcut. Need unit test in any case.
 * If one (but not all) input is scalar, fill the corresponding 1/8 of
 * the output with the scalar value and don't call any decimation
 * algorithm. This short cut will probably be required once the brick
 * API is fully implemented and can return scalars. This short cut
 * might belong in createLod() and not here. Also do a similar quick
 * handling of histogram and statistics.
 *
 * TODO-WIP-BrickedAPI: Possible shortcut.
 * After producing one regular output brick, check to see if it can be
 * replaced with a scalar. Not really needed in plan A, but plan C
 * might benefit a lot. In plan A, the write and read-back will take
 * care of it. The test is cheap if it fails, because isAllSame()
 * should return very quickly. And if it succeeds, it just proved
 * itself to be useful.
 *
 * TODO-WIP-BrickedAPI: If the user writes a constant value to the
 * entire survey, then all bricks, even those which are normally
 * padded, will be scalars. The problem is that bricks completely
 * outside the survey will be scalar with the default value which
 * is normally zero, while inside bricks will have the user's
 * scalar. This would lead to regular bricks in low res data even
 * though there are none in lod0.
 *
 * TODO-WIP-BrickedAPI: For incremental lod generation, a significant
 * improvement is possible in plan A. If 1 to 7 of a group of 8 bricks
 * is flagged as clean, then do a read/modify/write by reading the old
 * values into the output buffer, and then skippig the green input
 * bricks. Note that this gets a bit complicated, because the reads
 * from the output lod level ought to be grouped together into a single
 * read. Tip: Keep the vectors the same size as today. Use a pos of
 * (0, 0, -bs[0]) as a placeholder for "do not read". This works
 * better if read can return actual scalar buffers, which should be
 * fairly straight forward in the outside survey case. Not implemented
 * as of this comment.
 */

/**
 * Generate low resolution bricks only when all inputs are available.
 *
 * Thread safety: not thread safe.
 * The function is intended to be called from writeAlignedBrickList()
 * or equivalent, and those functions are also single threaded.
 * TODO-WIP-BrickedAPI: Locking anyway, just in case?
 */
void
GenLodImpl::tryToCall(
     const LevelAndForce& lodmode,
     const std::vector<index3_t>& brickroi,
     const std::function<bool(std::int64_t,std::int64_t)>& progress)
{
  const int loglevel = openzgy_timers() && lodmode.force >= 2 ? 0 : 3;
  const index3_t bs = _bricksize;
  const std::int64_t in_brickcolumn_lod1 =
    (_surveysize[2] + (2 * bs[2] - 1)) / (2 * bs[2]);
  const std::int64_t all_buffer_count = 9 * in_brickcolumn_lod1;

  // It is safer to ask the accessor for the tracker in each tryToCall()
  // to make sure the pointer doesn't go stale. Instead of holding on
  // to the pointer. Currently, only the GenLodC leaf type holds a
  // pointer to the accessor. So another virtual function is needed.
  // Sigh. The spaghetti factor of this code is getting too high.

  std::shared_ptr<TrackTouched> touched = _getNewTouched();
  if (touched) {
    TrackTouched::tasklist_t tasklist =
      touched->getWorkAndClear(lodmode, brickroi);
    SendProgress prog(progress);
    prog.start(tasklist);
    if (_logger(loglevel, "")) {
      if (tasklist.empty())
        _logger(loglevel, std::stringstream().flush()
                << "There are no chunks ready for lod computation."
                << " maxlevel = " << lodmode.level
                << " force = " << std::boolalpha << lodmode.force
                << (progress ? " (report progress)" : ""));
      else
        _logger(loglevel, std::stringstream().flush()
                << "There are " << tasklist.size()
                << " chunks ready for lod computation."
                << " maxlevel = " << lodmode.level
                << " force = " << std::boolalpha << lodmode.force
                << (progress ? " (report progress)" : ""));
    }
    if (!tasklist.empty()) {
      // TODO-WIP-BrickedAPI: Better re-use of buffers.
      std::vector<std::shared_ptr<DataBuffer>> all_buffers;
      for (int ii = 0; ii < all_buffer_count; ++ii)
        all_buffers.push_back(DataBuffer::makeNewBuffer3d(_bricksize, _dtype));
      for(const TrackTouched::task_t& task : tasklist) {
        if (_logger(3, ""))
          _logger(3, std::stringstream().flush()
                  << "  process lod " << task.lod << " brick"
                  << " (" << task.pos[0]
                  << ", " << task.pos[1]
                  << ", " << task.pos[2]
                  << ")");
        const LodAlgorithm algorithm =
          _decimation[std::min(task.lod-1, (int)_decimation.size()-1)];
        // "lod" is the output level; input is "lod-1"
        // positions should be relative to "lod".
        processFourInputColumnsOneLod
          (task.lod, task.pos[0], task.pos[1], task.zsize, algorithm, all_buffers);
        prog.step(task);
      }
    }
    prog.finish();
  }
  else {
    if (_logger(loglevel, ""))
      _logger(loglevel, std::stringstream().flush()
              << "tryToCall(): No TrackTouched instance");
  }
  // TODO-WIP-BrickedAPI: Test this.
  if (lodmode.force >= 2)
    clearLodTimers(true);
}

/**
 * Compute and write out low resolution bricks for level "lod" and coordinates
 * i.e. brick numbers in "lod" (pos0, pos1, 0) to (pos0, pos1, zcount-1).
 */
void
GenLodImpl::processFourInputColumnsOneLod(
     const int out_lod,
     const std::int64_t pos0,
     const std::int64_t pos1,
     const std::int64_t zcount,
     const LodAlgorithm algorithm,
     const std::vector<std::shared_ptr<DataBuffer>>& all_buffers)
{
  const index3_t bs = _bricksize;
  if ((std::int64_t)all_buffers.size() < 9 * zcount)
    throw OpenZGY::Errors::ZgyInternalError("Buffer overrun prevented");
  std::vector<std::shared_ptr<DataBuffer>> level_outbuffers
    (all_buffers.begin(), all_buffers.begin() + zcount);
  std::vector<std::shared_ptr<DataBuffer>> level_inbuffers
    (all_buffers.begin() + zcount, all_buffers.begin() + 9 * zcount);
  std::vector<index3_t> ipos;
  std::vector<index3_t> opos;
  // Handle 4 input brick columns, 1 half height output column.
  for (std::int64_t pos2 = 0; pos2 < zcount; ++pos2) {
    opos.push_back(index3_t{pos0*bs[0],pos1*bs[1],pos2*bs[2]});
    // Handle 8 input bricks, one output brick.
    for (std::int64_t subi = 0; subi < 2; ++subi) {
      for (std::int64_t subj = 0; subj < 2; ++subj) {
        for (std::int64_t subk = 0; subk < 2; ++subk) {
          ipos.push_back(index3_t
            {(2*pos0+subi)*bs[0], (2*pos1+subj)*bs[1], (2*pos2+subk)*bs[2]});
        }
      }
    }
  }
  // Handle 4 input brick columns, 1 half height output column.
  std::vector<std::shared_ptr<DataBuffer>> indata =
    _readbricks(ipos, level_inbuffers, out_lod - 1);
  // Handle the case where all 8 inputs of the lowres brick have
  // the same constant value.
  // TODO-WIP-BrickedAPI: Handle the case where some bricks are
  // completely outside the survey, and disagree about sample values.
  std::vector<bool> same;
  for (std::int64_t ix = 0; ix < zcount; ++ix) {
    same.push_back(areAllEightSameConst(indata, ix * 8));
  }
  WorkOrderRunner::parallelFor(8 * zcount,
    [this,&level_outbuffers,&indata,algorithm,&same]
    (std::int64_t part)
    {
      // Handle 1 input brick, saving the result to 1/8th of an output brick.
      // Or, handle 8 identical scalar bricks by making the output scalar.
      if (same[part/8]) {
        if ((part % 8) == 0) {
          level_outbuffers[part/8] = indata[part];
        }
      }
      else {
        decimateOneInputBrick(level_outbuffers, indata, part, algorithm);
      }
    });
  // Write 1 half height output column made from 4 input columns.
  _writebricks(opos, level_outbuffers, out_lod);
}

/**
 * Generate and store low resolution bricks, histogram, and statistics.
 * See doc/lowres.html for details. I/O is done via ZgyInternalBulk.
 * Use this class as part as finalize().
 */
GenLodC::GenLodC(
     const std::shared_ptr<ZgyInternalBulk>& accessor,
     const std::shared_ptr<ZgyInternalMeta>& meta,
     const compressor_t& lodcompressor,
     const std::vector<LodAlgorithm>& decimation,
     const std::shared_ptr<InternalZGY::HistogramData>& wa_histogram,
     LoggerFn logger)
  : GenLodImpl(meta->ih().size(),
               meta->ih().bricksize(),
               meta->ih().datatype(),
               meta->ih().nlods(),
               decimation,
               wa_histogram,
               meta->ih().defaultstorage(),
               /*progress=*/nullptr,
               logger)
  , _accessor(accessor)
  , _lodcompressor(lodcompressor)
{
  _logger(1, "GenLodC has been created for piecewise write.");
}

GenLodC::~GenLodC()
{
}

/**
 * \brief Allow querying for which bricks are ready for low resolution.
 */
std::shared_ptr<TrackTouched>
GenLodC::_getNewTouched() const
{
  auto accessor = this->_accessor.lock();
  if (accessor)
    return accessor->newTrackedBricks();
  else
    return std::shared_ptr<TrackTouched>{};
}

std::vector<std::shared_ptr<DataBuffer>>
GenLodC::_readbricks(
     const std::vector<std::array<int64_t,3>>& position,
     const std::vector<std::shared_ptr<DataBuffer>>& data_in,
     int lod) const
{
  std::vector<std::shared_ptr<DataBuffer>> data(data_in);
  auto accessor = this->_accessor.lock();
  if (!accessor)
    throw OpenZGY::Errors::ZgyInternalError("GenLod readbricks: no bulk");
  // TODO-WIP-BrickedAPI: The old code passed check_constant=true.
  // Should this be done here as well? Note that this might only
  // matter in obscure cases. A similar dilemma is seen when
  // re-opening a file. The check for an all-constant file is
  // just doing the quick check. See more explanation in
  // ZgyInternalBulk::newTrackedBricksTryEnable()

  if (!data.empty()) {
    data = accessor->readBricksToNewBuffers(position, data, lod, /*as_float=*/false, /*check_constant=*/false);
  }
  return data;
}

void
GenLodC::_writebricks(
     const std::vector<std::array<int64_t,3>>& position,
     const std::vector<std::shared_ptr<DataBuffer>>& data,
     int lod) const
{
  auto accessor = this->_accessor.lock();
  if (!accessor)
    return;
  if (!data.empty()) {
    accessor->writeBricksInternal
      (position, data, lod,
       /*is_storage=*/true, _lodcompressor,
       /*immutable_buffer=*/false,
       /*maybe_more=*/false);
  }
}

} // namespace
