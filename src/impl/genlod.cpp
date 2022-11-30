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
#include "exception.h"

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
    const std::array<double,2>& histogram_range,
    std::int32_t nlods_in,
    const std::vector<LodAlgorithm>& decimation,
    const std::shared_ptr<HistogramData>& histogram,
    double defaultstorage,
    bool incremental,
    bool skip_histogram,
    const std::function<bool(std::int64_t,std::int64_t)>& progress,
    LoggerFn logger)
  : _nlods(nlods_in)
  , _total(0)
  , _done(0)
  , _surveysize(size)
  , _bricksize(bricksize)
  , _dtype(dtype)
  , _histogram_range(histogram_range) // default (-1,1)
  , _decimation(decimation.empty() ? std::vector<LodAlgorithm>{
      LodAlgorithm::LowPass, LodAlgorithm::WeightedAverage} : decimation)
  , _wa_histogram(histogram) // Might point to self._histo later.
  , _wa_defaultstorage(defaultstorage)
  , _incremental(incremental)
  , _skip_histogram(skip_histogram)
  , _progress(progress)
  , _loggerfn(logger)
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

/**
 * Invoke the user's progress callback if any.
 * Keep track of how many bricks we have processed. Both reads and
 * writes increment the same counter. For plan C the reads will cover
 * all blocks in lod 0 and the writes will cover all blocks in lod > 0.
 * For plan D all blocks are written which means the computation of
 * _total done in __init__ might need to change.
 */
void
GenLodBase::_report(const DataBuffer* data) const
{
  if (data) {
    std::array<std::int64_t,3> bricks =
      (data->size3d() + this->_bricksize - std::int64_t(1)) / this->_bricksize;
    const_cast<std::int64_t&>(this->_done) += (bricks[0] * bricks[1] * bricks[2]);
  }
  if (this->_progress && !this->_progress(_done, _total))
    throw OpenZGY::Errors::ZgyAborted("Computation of low resolution data was aborted");
}

/**
 * Return the number of bricks that will be accessed (either read or written)
 * in order to generate low resolution bricks. If the version in the base
 * class is not overridden the assumption is that no lowres will be stored.
 *
 * Note that since this is a virtual method it must not be invoked from
 * the base class constructor.
 */
std::int64_t
GenLodBase::_willneed() const
{
  return 0;
}

/**
 * \brief Return true if the region exists on file and has not changed.
 * \details
 * The method may be overridden to implement imcremental rebuild where
 * some low resolution data might also be available on the file.
 */
bool
GenLodBase::_isclean(
    std::int32_t lod,
    const std::array<std::int64_t,3>& /*pos*/,
    const std::array<std::int64_t,3>& /*size*/) const
{
  return lod==0;
}

/**
 * This is a stub that must be redefined except for low level unit tests.
 * Read a block from the ZGY file (plans B and C) or the application
 * (plan D). The lod parameter will always be 0 for plans C and D.
 * Returns a ScalarBuffer if all constant, else a 3D numpy array.
 */
std::shared_ptr<DataBuffer>
GenLodBase::_read(
    std::int32_t lod,
    const std::array<std::int64_t,3>& pos,
    const std::array<std::int64_t,3>& size) const
{
  std::shared_ptr<DataBuffer> result = DataBuffer::makeScalarBuffer3d
    (_wa_defaultstorage, size, _dtype);
  _report(result.get());
  return result;
}

/**
 * This is a stub that must be redefined except for low level unit tests.
 * Write a block to the ZGY file. Silently ignore writes of data that
 * is known to have been read directly from the file. For plans B and C
 * this means ignoring all writes to lod 0.
 */
void
GenLodBase::_write(
    std::int32_t lod,
    const std::array<std::int64_t,3>& pos,
    const std::shared_ptr<const DataBuffer>& data) const
{
  this->_report(data.get());
}

/**
 * Given a starting position and a size, reduce the size if needed so
 * it fits inside the survey. The result may be empty. Note that there
 * is no need for a correspinding method to adjust the position because
 * padding is always added at the end of the survey.
 */
index3_t
GenLodBase::_clipsizetosurvey(std::int32_t lod, const index3_t& pos, const index3_t& size) const
{
  const std::int64_t lodfactor = std::int64_t(1) << lod;
  const index3_t surveysize =
    (this->_surveysize + (lodfactor - 1)) / lodfactor;
  // Don't use std::min(index3_t,index3_t). It doesn't do what you think.
  return index3_t{
    std::max(std::int64_t(0), std::min(pos[0]+size[0], surveysize[0]) - pos[0]),
    std::max(std::int64_t(0), std::min(pos[1]+size[1], surveysize[1]) - pos[1]),
    std::max(std::int64_t(0), std::min(pos[2]+size[2], surveysize[2]) - pos[2])
  };
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
    const std::array<double,2>& histogram_range,
    std::int32_t nlods_in,
    const std::vector<LodAlgorithm>& decimation,
    const std::shared_ptr<HistogramData>& histogram,
    double defaultstorage,
    bool incremental,
    bool skip_histogram,
    const std::function<bool(std::int64_t,std::int64_t)>& progress,
    LoggerFn logger)
  : GenLodBase(size, bricksize, dtype, histogram_range, nlods_in, decimation,
               histogram, defaultstorage, incremental, skip_histogram, progress, logger)
  , _stats(std::make_shared<StatisticData>())
  , _histo(std::make_shared<HistogramData>(256, histogram_range[0], histogram_range[1]))
{
  // See base constructor. Better than no histogram at all.
  // The WeightedAverage algorithm will now rely on the statistics
  // calculated so far. The exact result thus depends on the order
  // the data bricks is written in.
  if (!this->_wa_histogram || this->_wa_histogram->getcount() == 0)
    this->_wa_histogram = this->_histo;

  // LODs with WeightedAverage need the histogram. It is less disruptive
  // to just ignore the BuildNoHistogram argument because that would
  // be the safe option.
  if (_skip_histogram) {
    for (LodAlgorithm alg : decimation) {
      if (alg == LodAlgorithm::WeightedAverage) {
        _logger(0, "WARNING: NoHistogram ignored because of WeightedAverage");
        _skip_histogram = false;
        break;
      }
    }
  }

#if 0
  // Might want to disable this until BuildNoHistogram can be tested better.
  // Leaving the instances behind will result in incorrect data but fewer
  // corner cases to test. The bad data should be discarded by the caller.
  if (_skip_histogram) {
    _stats.reset();
    _histo.reset();
    _wa_histogram = std::make_shared<HistogramData>(256, -1, 1);
  }
#endif
}

/**
 * Generate and store statistics, histogram, and all low resolution
 * bricks. Works for plans C and D. If we also need an implementation
 * of plan B then this method wold need to iterate over all bricks
 * and lods, and _calculate would not make any recursive calls.
 *
 * TODO-Performance: If the bulk layer is made thread safe for writing
 * it is possible to do parallel processing at this high level.
 * E.g. do this for just one level: Split into 4 sub-tasks that
 * each execute in one thread. Special handling will be needed
 * for the simgle highest-level brick. All 4 threads need to be
 * joined before that one can be done. With one level here and with
 * 4 threads used in _calculate() this means we will be reading with
 * 16 threads. But there are serious caveats:
 *
 * \li Using 4 threads here means there is just 1/4th of the memory
 *     available to each thread, which means requests may become smaller
 *     which means we might not see any benefit at all.
 *
 * \li The added complexity both here and in making the bulk write
 *     thread safe might not be worth the trouble.
 */
std::tuple<std::shared_ptr<StatisticData>, std::shared_ptr<HistogramData>>
GenLodImpl::call()
{
  if (_logger(4, ""))
    _logger(4, std::stringstream()
            << "@ GenLod is running."
            << " Histogram range " << fmt(_histogram_range));
  // Keep this in sync with GenLodC::_willneed(). If incremental is true
  // it is not trivial to use anything other than one brick. which gets
  // changed to one brick-column elsewhere.
  index3_t chunksize = this->_bricksize * std::int64_t(_incremental ? 1 : 2);
  this->_reporttotal(_willneed());
  this->_report(nullptr);
  this->_calculate(index3_t{0,0,0}, chunksize, this->_nlods-1);
  return std::make_tuple(this->_stats, this->_histo);
}

template <typename T>
void
GenLodImpl::_accumulateT(const std::shared_ptr<const DataBuffer>& data_in)
{
  std::shared_ptr<const DataBufferNd<T,3>> data =
    std::dynamic_pointer_cast<const DataBufferNd<T,3>>(data_in);
  if (!data)
    return;
  if (!this->_stats || !this->_histo)
    return; // Probably because _skip_histogram was set.
  // Note that we never allow the histogram to grow dynamically;
  // that was a problematic feature of the old accessor.
  HistogramBuilder hb(256, _histogram_range[0], _histogram_range[1]);
  std::int64_t len = data->size3d()[0] * data->size3d()[1] * data->size3d()[2];
  if (this->_skip_histogram) {
    // Keep track of very coarse statistics and histogram, by using just
    // the first sample of each brick and pretending this is const value.
    // This gets rid of 99.99% of the overhead but leaves enough to avoid
    // some corner cases. NOTE that the bogus statistics should be removed
    // as soon as genlod finishes. See GenLodImpl::GenLodImpl() for whether
    // we can ever get here.
    hb.add(data->data(), data->data() + 1);
    hb *= len;
  }
  else if (data->isScalar()) {
    hb.add(data->data(), data->data() + 1);
    hb *= len;
  }
  else {
    hb.add(data->data(), data->data() + len);
  }
  *this->_stats += hb.getstats();
  *this->_histo += hb.gethisto();
}

/**
 * Keep a running tally of statistics and histogram.
 *
 * When doing an incremental build this data will be discarded
 * so it might as well not be collected it in the first place.
 */
void
GenLodImpl::_accumulate(const std::shared_ptr<const DataBuffer>& data)
{
  if (!data)
    return;
  if (_incremental)
    return;
  switch (data->datatype()) {
  case RawDataType::SignedInt8:    _accumulateT<std::int8_t>(data); break;
  case RawDataType::UnsignedInt8:  _accumulateT<std::uint8_t>(data); break;
  case RawDataType::SignedInt16:   _accumulateT<std::int16_t>(data); break;
  case RawDataType::UnsignedInt16: _accumulateT<std::uint16_t>(data); break;
  case RawDataType::SignedInt32:   _accumulateT<std::int32_t>(data); break;
  case RawDataType::UnsignedInt32: _accumulateT<std::uint32_t>(data); break;
  case RawDataType::Float32:       _accumulateT<float>(data); break;
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Read or recursively compute a chunk (usually 1 or 4 brick-columns)
 * of data at the specified (readpos, readsize, readlod) and return a
 * memory buffer 1/8th of the input size holding a decimated version.
 * Also store the chunk that was obtained (i.e. not the decimated data)
 * back to file, if it wasn't read from file in the first place.
 * If the decimated data is to be saved then the caller must do that.
 *
 * When asked to read (and implicitly write) the highest LOD level then
 * nothing is returned because no more decimation is to be done.
 *
 * The vertical start and size is ignored. Full vertical traces are read.
 * The provided readpos and readsize must both be brick aligned.
 * The actual size of the returned buffer will be clipped to the survey
 * size and might even end up empty. The returned buffer has no padding.
 *
 * When doing an incremental build, if calculate is called with a clean
 * region it will not need to do a recursive call and it will not need
 * to write anything. It still needs to return an in-memory buffer of
 * decimated data but it can get that by reading 1/8th of a chunk of
 * decimated data from the lod+1 level. Because when the chunk is clean
 * then the corresponding part of the lowres brick is still valid.
 * Quality might suffer in incremental mode if the output is compressed.
 *
 * TODO-@@@: Performance: Incremental build: Reading more data than needed.
 * Assuming the chunk size is one brick-column which means that any chunk C
 * not on the survey edge will have 4 brick-columns as input. If two or
 * three of those inputs are marked as clean then there will eventually be
 * two or three reads from different parts of C in order to get hold of
 * the still valid decimated data it contains. Each read will read just
 * one quarter of a brick column. So the low level I/O routine will need
 * to read the full brick column and discard 3/4 or the data. The net
 * effect is that the I/O layer may be reading the same data 2 or 3 times.
 * Note that the issue might not be noticeable in practice because it
 * is more likely for a chunk to have either zero dirty inputs, i.e. no
 * calculation needed, or all dirty inputs.
 * If a fix really is needed, consider:
 *
 *   - Rewrite the two lines handling recursive calls to something a lot
 *     more complicated. Implementing some kind of r/m/w logic.
 *   - Set chunk size to 4 brick columns. Cure may be worse than the disease.
 *   - Some kind of caching.
 *
 * If doing an incremental build, readsize is best set to one brick-column
 * because this determines the granularity of the "dirty" test and that
 * needs to be as small as possible. But see above; this might mean
 * the code will in some cases read data smaller than one brick. TODO-@@@
 *
 * TODO-Performance: If doing a full build it might we a good idea to allow
 * the application to configure how much memory the computation is allowed
 * to use. Increase readsize accordingly. Larger bricks might help the bulk
 * layer and the multi-threaded decimation routines become more efficient.
 * The gain might not be noticeable though.
 *
 * Full resolution data (lod 0) will be read from file (plan C) or the
 * application (plan D). Low resolution is computed by a recursive call
 * to this function (plans C and D) or by reading the file (plan B).
 * Note that currently only plan C is fully implemented.
 *
 * For plans B and C a single call needs to be made to read the brick
 * (there is by definition just one) at the highest level of detail.
 * This will end up computing all possible low resolution bricks and
 * storing them. For plan B the caller must iterate.
 *
 * The function is also responsible for collecting statisics and
 * histogram data when doing a full build. For incremental builds
 * it becomes the callers repsonsibility to track changes as data
 * is written.
 *
 * Note that some of the decimation algorithms need the histogram of the
 * entire file. Ideally the histogram of the entire file should be
 * available before decimation starts but that is impractical. At least
 * make sure the histogram update is done early enough and the decimation
 * late enough that the chunk of data being decimated has already been
 * added to the histogram. TODO-Low: Read e.g. 5% of all bricks up front
 * to get a better approximation of the histogram and then use that
 * result for the enite lowres generation.
 */
std::shared_ptr<DataBuffer>
GenLodImpl::_calculate(const index3_t& readpos_in, const index3_t& readsize_in, std::int32_t readlod)
{
  const std::int64_t lodfactor = std::int64_t(1) << readlod;
  const std::array<std::int64_t,3> surveysize =
    (this->_surveysize + (lodfactor - 1)) / lodfactor;
  const std::array<std::int64_t,3> readpos{readpos_in[0], readpos_in[1], 0};
  if (readpos[0] >= surveysize[0] || readpos[1] >= surveysize[1])
    return nullptr;
  const std::array<std::int64_t,3> chunksize
    {readsize_in[0], readsize_in[1], surveysize[2]}; // Read full traces.
  const std::array<std::int64_t,3> readsize
    {std::min(chunksize[0], (surveysize[0] - readpos[0])),
     std::min(chunksize[1], (surveysize[1] - readpos[1])),
     chunksize[2]};
  // Size of the in-memory buffer of decimated data returned to the caller.
  const std::array<std::int64_t,3> returnsize = (readsize + std::int64_t(1)) / std::int64_t(2);
  const std::array<std::int64_t,3> returnpos = readpos / std::int64_t(2);

  if (_logger(4, ""))
  _logger(4, std::stringstream()
          << "@" << _prefix(readlod)
          << "calculate(lod=" << readlod
          << ", pos=" << fmt(readpos)
          << ", size=" << fmt(readsize) << ")");

  std::shared_ptr<const DataBuffer> data;
  bool wasread = false;
  if (this->_isclean(readlod, readpos, readsize)) {
    std::shared_ptr<DataBuffer> result = (readlod == this->_nlods-1) ? nullptr :
      this->_read(readlod+1, returnpos, returnsize);
    if (_logger(4, "")) {
      _logger(4, std::stringstream()
              << "@" << _prefix(readlod)
              << "calculate returns(lod="
              << readlod+1 << ", pos=" << fmt(returnpos)
              << ", size=" << fmt(returnsize)
              << ", data=" << result->toString() << " (SHORTCUT))");
    }
    return result;
  }
  else if(readlod == 0) {
    // Fullres bricks are always read, not calculated.
    data = this->_read(readlod, readpos, readsize);
    wasread = true;
    this->_accumulate(data);
  }
  else {
    std::array<std::int64_t,3> offsets[4] =
      {{           0,            0,   0},
       {           0, chunksize[1],   0},
       {chunksize[0],            0,   0},
       {chunksize[0], chunksize[1],   0}};
    std::shared_ptr<DataBuffer> hires[4]{nullptr, nullptr, nullptr, nullptr};

    // Compute the requested result by recursively reading 4 chunks
    // at lod-1 and gluing the result together.
    //
    // This loop should normally not be parallelized or consolidated into a
    // single call reading 4x the amount of data. The serial loop is what
    // prevents the recursion from trying to read the entire file into memory.
    //
    // TODO-Performance: A possible exception in algorithm "C" is when
    // readlod==1. Caveat: For incremental builds we might not need all 4
    // sub-parts. So in that case we should loop normally.
    // Caveat for multi threading: nested loops, smaller blocks.
    // Caveat for replacing with a single call: more special cases to test.
    // In particular handling of crops to survey size.

    for (int ii=0; ii<4; ++ii)
      hires[ii] = this->_calculate(readpos*std::int64_t(2) + offsets[ii],
                                   chunksize, readlod-1);
    data = this->_paste4(hires[0], hires[1], hires[2], hires[3]);
  }

  if (!wasread)
    this->_write(readlod, readpos, data);

  std::shared_ptr<DataBuffer> result;
  if (readlod == this->_nlods - 1) {
    result = nullptr; // Caller will discard it anyway.
    if (this->_done != this->_total) {
#if 1 // If trusting _willneed()
      throw OpenZGY::Errors::ZgyInternalError
        ("GenLodImpl: Expected " + std::to_string(this->_total) +
         " reads and writes but saw " + std::to_string(this->_done) + ".");
#else
      _logger(0, std::stringstream()
              << "Warning: GenLodImpl: Expected " << this->_total
              << " reads and writes but saw " << this->_done << ".\n");
      this->_done = this->_total;
      this->_report(nullptr);
#endif
    }
  }
  else if (data->isScalar() || data->isAllSame(nullptr)) {
    // The test for isAllSame() is just a performance tweak.
    // It normally returns false and in that case it normally
    // returns very quickly. So there should be no downside.
    // data is shrink-to-fit i.e. has no padding at the survey
    // edge so isAllSame() doesn't need the "used" parameter.
    result = DataBuffer::makeScalarBuffer3d
      (data->scalarAsDouble(), returnsize, data->datatype());
  }
  else {
    // TODO-Performance: Enable parallelization (in genlod.cpp) of this call.
    // Note that for readlod==0 we might already be running in one of four
    // OpenMP threads. It is not clear whether to force OpenMP to also
    // parallelize this low level loop. // omp_set_nested(1) or better, use
    // omp_set_max_active_levels(...).
    result = this->_decimate(data, readlod+1);
    if (result->size3d() != returnsize) {
      _logger(0, std::stringstream()
              << "Decimation returned " << fmt(result->size3d())
              << ", expected " << fmt(returnsize)
              << ", start " << fmt(readpos)
              << ", size=" << fmt(readsize) << "\n");
      throw OpenZGY::Errors::ZgyInternalError("GenLodImpl: Wrong returnsize.");
    }
  }
  if (_logger(4, ""))
    _logger(4, std::stringstream()
            << "@" << _prefix(readlod)
            << "calculate returns(lod="
            << readlod+1 << ", pos=" << fmt(readpos) << " / 2"
            << ", size=" << fmt(returnsize)
            << ", data=" << data->toString() << ")");
  return result;
}

/**
 * Return a decimated version of the input buffer with half the size
 * (rounded up) in each dimension. In total the result will be ~1/8
 * the size of the input.
 *
 * Lod refers to the level being generated. Must be >= 1.
 */
std::shared_ptr<DataBuffer>
GenLodImpl::_decimate(const std::shared_ptr<const DataBuffer>& data, std::int64_t lod)
{
  if (!data)
    return nullptr;
  if (lod < 1)
    throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_decimate must have lod >= 1");

  const std::array<std::int64_t,3> outsize =
    (data->size3d() + std::int64_t(1)) / std::int64_t(2);
  if (data->isScalar()) {
    double value = data->scalarAsDouble();
    return DataBuffer::makeScalarBuffer3d(value, outsize, data->datatype());
  }
  else {
    std::shared_ptr<DataBuffer> result =
      DataBuffer::makeNewBuffer3d(outsize, data->datatype());
    // Constructor guaranteed !_decimation.empty() but might not be big enough.
    const LodAlgorithm algorithm =
      (lod <= (int)_decimation.size()) ? _decimation[lod-1] : _decimation.back();

    // The 4 last arguments are only used for vvArrayBasic::lod_WeightedAverage.
    // Note that the value range for the histogram should at this
    // point refer to storage and not converted values. For integral
    // types the value range will typically be the natural range of
    // int8_t or int16_t.
    // See BrickedAccessor.cpp in old code
    // See createGenericLevelOfDetail() regarding shrink-to-fit survey,
    // c-contiguous buffers, and odd survey size.
    // TODO-Medium missing _wa_defaultstorage.
    createLod(result, data, algorithm,
              _wa_histogram->getbins(),
              _wa_histogram->getsize(),
              _wa_histogram->getmin(),
              _wa_histogram->getmax());
    return result;
  }
}

/**
 * See _paste4() for details.
 */
std::shared_ptr<DataBuffer>
GenLodImpl::_paste1(
    const std::shared_ptr<DataBuffer>& result,
    const std::shared_ptr<const DataBuffer>& more,
    std::int64_t ioff, std::int64_t joff)
{
  if (more) {
    const index3_t dstorig{0,0,0};
    const index3_t srcorig{ioff, joff, 0};
    result->copyFrom(more.get(), srcorig.data(), dstorig.data(), 0, 0);
  }
  return result;
}

/**
 * Combine 4 buffers into one. Input buffers may be None (do not
 * paste) or ScalarBuffer (paste a constant value). If all not-None
 * buffers are just scalars then the return from this function
 * will also be a scalar. d01 adds more data in the J direction,
 * so it starts at i=0, j>0 in the target. Similarly d10 adds
 * more in the J direction. And d11 in the diagonal.
 *
 * Performance note: It is in theory possible to avoid some buffer
 * copies, this one in particular, by passing our 4 or 8 buffers
 * to the decimation algorithm instead of a combined buffer.
 * The algorithms remain just as efficient. BUT if sizes can vary
 * or bricks can be missing or number of bricks differs in level 1
 * because we read directly from the file then things can get really
 * complicated really fast.
 */
std::shared_ptr<const DataBuffer>
GenLodImpl::_paste4(
    const std::shared_ptr<const DataBuffer>&d00,
    const std::shared_ptr<const DataBuffer>&d01,
    const std::shared_ptr<const DataBuffer>&d10,
    const std::shared_ptr<const DataBuffer>&d11)
{
  if (!d01 && !d10 && !d11)
    return d00; // Nothing to paste. Also works for empty or scalar.

  // Empty d00 + non-empty others is not good.
  if (!d00)
    throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_paste4() assert#1.");
  if (d01 && d01->size3d()[0] != d00->size3d()[0])
    throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_paste4() assert#2.");
  if (d10 && d10->size3d()[1] != d00->size3d()[1])
    throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_paste4() assert#3.");
  if (d01 && d10) {
    // The "diagonal" brick must exist with the right size.
    if (!d11 ||
        d11->size3d()[1] != d01->size3d()[1] ||
        d11->size3d()[0] != d10->size3d()[0])
      throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_paste4() assert#4.");
  }
  else {
    // The "diagonal" brick should not be needed in this case.
    if (d11)
      throw OpenZGY::Errors::ZgyInternalError("GenLodImpl::_paste4() assert#5.");
  }

  const std::int64_t ni = d00->size3d()[0] + (d10 ? d10->size3d()[0] : 0);
  const std::int64_t nj = d00->size3d()[1] + (d01 ? d01->size3d()[1] : 0);
  const std::int64_t nk = d00->size3d()[2];

  bool all_same = true;
  const std::array<std::shared_ptr<const DataBuffer>,4> all{d00, d01, d10, d11};
  for (const auto& e : all)
    if (all_same && e)
      if (!e->isScalar() || e->scalarAsDouble() != d00->scalarAsDouble())
        all_same = false;

  std::shared_ptr<DataBuffer> result;
  if (all_same) {
    double value = d00->scalarAsDouble();
    result = DataBuffer::makeScalarBuffer3d(value,
                                          index3_t{ni,nj,nk}, d00->datatype());
  }
  else {
    // TODO-Test: needs careful testing with weird sizes.
    result = DataBuffer::makeNewBuffer3d(index3_t{ni,nj,nk}, d00->datatype());
    _paste1(result, d00, 0, 0);
    _paste1(result, d01,                0, d00->size3d()[1]);
    _paste1(result, d10, d00->size3d()[0],                0);
    _paste1(result, d11, d00->size3d()[0], d00->size3d()[1]);
  }
  return result;
}

/**
 * Choose the histogram range to use.
 *
 * [1] For cubes stored as integral data use the entire range that can be
 * represented. Not just the possibly smaller range samples written.
 * For int8/uint8 and a 256-bin histogram this is a no-brainer because
 * having less than one possible sample value in each bin inevitably
 * leads to some empty bins even for a completely smooth distribution
 * of the input. For int16/uint16 the strategy is still workable but
 * it had been better to use a compromise: Use a narrower range but
 * narrowed using an integer factor. Or sidestep the issue by
 * internally using a 64k histogram that will get trimmed down to 256
 * entries later.
 *
 * For cubes stored as float or compressed data it gets more complicated.
 *
 * The range of all sample values seen until now, i.e. everything
 * written, is passed in. The possibility to overwrite data or
 * (future) append to existing file makes this not accurate but still
 * probably good enough.
 *
 * [2] In the normal case the file contains a range of different
 *     samples. Map the center of the first bin of the histogram to
 *     the lowest written value and the center of the last bin of the
 *     histogram to the highest written value. If the written values
 *     include zero then adjust the range slightly to make the range
 *     zero-centric. TODO-Low implement zero-centric.
 *
 *     Note that if the distribution of sample values is perfectly
 *     uniform this results in the first and last bin will probably
 *     have lower counts then expected. This is partly because the
 *     initial mapping is to the center of the bin instead of the
 *     outer edge, and partly because of the zero-centric adjustment.
 *     This is a very minor issue compared to how bad the histogram
 *     might look if it isn't zero centric.
 *
 *     One motivation for using bin centers instead of outer edges is
 *     that it helps (but does not guarantee) smooth histograms for a
 *     float cube that contains discrete values.
 *
 * [3] If writtenrange is invalid this means that no finite samples
 *     have been written. Choose an arbitrary range in that case
 *     instead of ending up with a lot of obscure corner cases. The
 *     range needs to include defaultvalue (which for float data is
 *     zero) as there might be unwritten bricks. The range (-1,+1) is
 *     probably as good as any. Or adjust it slightly to make it zero
 *     centric. Or use (-128,+127) giving a bin width of 1 if the
 *     histogram has 256 values.
 *
 * [4] The same applies when only zeros have been written. Use (-1,+1)
 *     or (-128,+127) All samples and up in the center bin which is
 *     128.
 *
 * [5] For robustness also use that value if some internal error
 *     occurred.
 *
 * [6] If only a single non-zero value has been written then choose a
 *     range with defaultvalue (i.e. zero) at one end and the written
 *     value at the other. Only bin 0 and 255 will have non-zero
 *     counts. Which bin contains defaultvalue and which contains the
 *     written value depends on the sign of the written value.
 *
 *     An alternative to the above would be to include 2*writtenvalue
 *     in the range, which would always map writtenvalue to bin 128
 *     just like the all-zero case does. This is arguably more
 *     consistent. But to avoid bias, writtenvalue should map to the
 *     center of a bin. Similar to the zero-centric propery which maps
 *     zero to the center of a bin. So the range should run to
 *     (255.0/128.0)*writtenvalue. Or (255.0/127.0) if writtenvalue is
 *     negative. Or something else if the size of the histogram is not
 *     256. Sigh. The previous choice sounds a lot better now, right?
 *
 * TODO-Medium: Make sure the Python version and the C++ version use
 * the same algorithm. In the Python implementation this logic is in
 * HistogramData.suggestHistogramRange().
 *
 * TODO-Low: If data is being appended the code will still re-compute
 * the entire histogram. To do this, it uses the _written_sample_min/max
 * kept track of while writing lod0 data. The range should be the union
 * of the data range written previously, found in the old histogram, and
 * the samples written this time around. Problem
 * is, the first write might not have bothered to finalize and thus
 * did not save this information. I probably need to "finalize" with
 * an empty histogram. Then include the existing histogram range in
 * the new one. Bear in mind that the stored range will be in user and
 * not storage values.
 *
 * TODO-Low should the suggested range for float data include default
 * value? Which for floats is always zero? Technically the code should
 * keep track of whether all bricks have been explicitly written and
 * if not include zero. It doesn't. This "bug" is probably in the
 * "extreme nitpicking" category. Application code can remove the
 * problem by explicitly choosing the defaultvalue to use, then
 * initialize the entire file with that value.
 *
 * TODO-Test: Unit tests for all these corner cases both C++ and
 * Python for files containing Nothing, 0, -42, +42. I can probably do
 * all those tests manually with zgycopy / zgycopyc but an automated
 * test is better.
 */
std::array<double,2>
GenLodImpl::suggestHistogramRange(
    const std::array<double,2>& writtenrange,
    RawDataType dtype)
{
  const std::array<double,2> bogus{-128,+127};
  const RawDataTypeDetails details(dtype);
  if (details.is_integer)
    return std::array<double,2>{details.lowest, details.highest}; // [1]
  else if (!std::isfinite(writtenrange[0]) ||
           !std::isfinite(writtenrange[1]) ||
           writtenrange[0] > writtenrange[1])
    return bogus; // [3], [5], nothing written or error.
  else if (writtenrange[0] < writtenrange[1])
    return writtenrange; // [2], normal case
  else if (writtenrange[0] > 0)
    return std::array<double,2>{0, writtenrange[0]}; // [6], single value
  else if (writtenrange[0] < 0)
    return std::array<double,2>{writtenrange[0], 0}; // [6], single value
  else
    return bogus; // [4], all zero
}

/**
 * Generate and store low resolution bricks, histogram, and statistics.
 * See doc/lowres.html for details. I/O is done via ZgyInternalBulk.
 * Use this class as part as finalize().
 * The implementation uses plan C, which means the full resolution data
 * will be read from the ZGY file. To implement plan D, make a derived
 * class that redefines _read() to query the client for the required full
 * resolution data. _read() must then also call _write() to store the
 * data it just received.
 */
GenLodC::GenLodC(
    const std::shared_ptr<ZgyInternalBulk>& accessor,
    const std::shared_ptr<ZgyInternalMeta>& meta,
    const compressor_t& lodcompressor,
    const std::vector<LodAlgorithm>& decimation,
    bool incremental,
    bool skip_histogram,
    const std::function<bool(std::int64_t,std::int64_t)>& progress,
    LoggerFn logger)
  : GenLodImpl(meta->ih().size(),
               meta->ih().bricksize(),
               meta->ih().datatype(),
               suggestHistogramRange(accessor->valueRangeWritten(),
                                     meta->ih().datatype()),
               meta->ih().nlods(),
               decimation,
               nullptr, // Will use computed histogram so far.
               meta->ih().defaultstorage(),
               incremental,
               skip_histogram,
               progress,
               logger)
  , _accessor(accessor)
  , _lodcompressor(lodcompressor)
{
  // Doing a full build even when incremental changes have been tracked
  // is ok; that is up to the application to choose. The inverse won't
  // work and should have been tested for earlier.
  if (incremental && accessor->trackedBricks().size() == 0)
    throw OpenZGY::Errors::ZgyInternalError("GenLodC inconsistent parameters.");
  if (_logger(4, ""))
    _logger(4, std::stringstream()
            << "@ GenLod is created."
            << " Written range   " << fmt(accessor->valueRangeWritten()));
}

/**
 * Return the number of bricks that will be accessed (either read or written)
 * in order to generate low resolution bricks. With plan C doing a full rebuild
 * all the lod0 bricks will be read exactly once and all the lowres bricks
 * will be written exactly once.
 *
 * Incremental builds are trickier. Fullres bricks do not directly affect
 * the read+write count. Any lowres brick marked as clean means that it
 * 8 input bricks will no longer be read (subtract 8) and that this
 * brick changes from potentially being written to potentially being read
 * (no change to the total). The brick might end up neither read nor
 * written if all its 7 siblings also are clean. But that adjustment is
 * done by the lowres brick above. With one exception: If the single
 * brick in the last level is clean, implying that the entire data set
 * is clean, this saves not just the input but also the brick itself
 * because there is never any reason to read it. Only write.
 *
 * To complicate the above, fewer than 8 bricks might be involved if
 * close to the survey border.
 *
 * To complicate even further, when _calculate chooses how much data
 * to process then it will read 4 brick-columns at "readlod".
 * I suspect one brick-column would have worked just as well.
 */
std::int64_t
GenLodC::_willneed() const
{
  std::int64_t total = 0; // Total number of bricks in all levels.
  {
    std::array<std::int64_t,3> bs = this->_bricksize;
    std::array<std::int64_t,3> sz = this->_surveysize;
    while (sz[0] > bs[0] || sz[1] > bs[1] || sz[2] > bs[2]) {
      std::array<std::int64_t,3> bricks = (sz + bs - std::int64_t(1)) / bs;
      total += (bricks[0] * bricks[1] * bricks[2]);
      sz = (sz + std::int64_t(1)) / std::int64_t(2);
    }
    ++total; // Loop stopped short of the last level, by definition one brick.
  }

  if (!_incremental)
    return total;

  // The incremental case is trickier.
  // The code below might have worked for full rebuilds as well,
  // but it currently assumes chunksize is one brick-column.
  const std::int64_t full_total = total;
  total = 0; // Start over. full_total is just for logging.

  const auto countbricks = [this](const index3_t& size) {
      return
        ((size[0] + this->_bricksize[0] - 1) / this->_bricksize[0]) *
        ((size[1] + this->_bricksize[1] - 1) / this->_bricksize[1]) *
        ((size[2] + this->_bricksize[2] - 1) / this->_bricksize[2]);
  };

  // How much data to be processed in each chunk, in samples?
  // Needs to match chunksize in GenLodImpl::call(), and will
  // currently ONLY works as shown, chunksize = one brick-column.
  // Other chunk sizes will need non trivial changes below.
  std::array<std::int64_t,3> chunksize
    {this->_bricksize[0],
     this->_bricksize[1],
     this->_surveysize[2]};
  // How much decimated data does this correspond to, in samples?
  const index3_t cshalf
    {(chunksize[0]+1)/2,
     (chunksize[1]+1)/2,
     (chunksize[2]+1)/2};
  // Number of bricks to read or write to get these samples?
  // Note that cshalf might need less than a brick, it still needs to read all.
  for (std::int32_t lod = 0; lod < _nlods; ++lod) {
    const std::int64_t lodfactor = std::int64_t(1) << lod;
    const index3_t sz = (this->_surveysize + (lodfactor - 1)) / lodfactor;
    for (std::int64_t ii = 0; ii < sz[0]; ii += chunksize[0]) {
      for (std::int64_t jj = 0; jj < sz[1]; jj += chunksize[1]) {
        for (std::int64_t kk = 0; kk < sz[2]; kk += chunksize[2]) {
          const index3_t pos{ii,jj,kk};
          const index3_t poshalf{ii/2, jj/2, kk/2}; // in lod+1
          const index3_t csclip =_clipsizetosurvey(lod, pos, chunksize);
          const index3_t cshalfclip =
            _clipsizetosurvey(lod+1, poshalf, cshalf);
          if (!_isclean(lod, pos, csclip)) {
            total += countbricks(csclip);
          }
          else if (lod < (_nlods-1) && !_isclean(lod+1, poshalf, cshalfclip)) {
            total += countbricks(cshalfclip);
            // READ THIS if you want to change chunksize.
            // The region checked using _isclean() should have been the
            // chunksize- sized region in lod+1 that is built from, in
            // part, the chunk being handled here. Because the code needs
            // to know whether that lod+1 region needs to be pocessed.
            // With chunksize exactly one brick-column the test on
            // poshalf, cshalfclip is equivalent since the dirty bricks
            // are tracked per brick so the code will be checking what
            // it needs to. If chunksize is larger then need to figure
            // out the start and size of that lod+1 brick. Remember
            // to take into account rounding up to align with bricksize,
            // clipping to survey size, rounding start down to chunksize,
            // and possibly something else. And possibly in a different
            // order. Possibly if academic interest only because the
            // chunk size of one brick-column might well be optimal
            // for incremental builds.
          }
        }
      }
    }
  }

  if (_logger(4, "")) {
    _logger(4, std::stringstream()
            << "@ Incremental lowres build needs " << total
            << " of " << full_total << " bricks");
  }

  return total;
}

/**
 * \brief Return true if the region exists on file and has not changed.
 *
 * \details
 *
 * If not doing an incremental build than all data is considered dirty.
 * Dirty LOD0 data needs to be read from file (plan C) or app (plan D),
 * Dirty LOD>1 data will need to be re-computed by recursively calculating it.
 *
 * An exception is thrown if any part of the region is outside the survey.
 * This can easily be changed to just ignore the parts that are outside
 * and to reurn a special status (so the return can no longer be a bool)
 * if completely outside. But with current usage the caller handles that
 * check anyway.
 */
bool
GenLodC::_isclean(
    std::int32_t lod,
    const std::array<std::int64_t,3>& pos,
    const std::array<std::int64_t,3>& size) const
{
  if (!_incremental) {
    return false;    // Not tracking dirty bits, or caller wants full anyway.
  }
  else if ((this->_accessor->trackedBricksDirty(pos, size, lod) & 0x01) != 0) {
    return false;    // Something dirty somewhere in the input.
  }
  else {
    return true;     // Hurra! We can read from file instead of recomputing,
  }
}

/**
 * See base class for details.
 */
std::shared_ptr<DataBuffer>
GenLodC::_read(
    std::int32_t lod,
    const index3_t& pos,
    const index3_t& size) const
{
  std::shared_ptr<DataBuffer> result = this->_accessor->readToNewBuffer
    (pos, size, lod, /*as_float=*/false, /*extra_checking*/true);
  if (_logger(4, ""))
    _logger(4, std::stringstream()
            << "@" << _prefix(lod)
            << "read(lod=" << lod << ", pos=" << fmt(pos)
            << ", " << result->toString() << ")");
  _report(result.get());
  return result;
}

/**
 * See base class for details.
 */
void
GenLodC::_write(
    std::int32_t lod, const index3_t& pos,
    const std::shared_ptr<const DataBuffer>& data) const
{
  if (lod > 0) {
    if (_logger(4, ""))
      _logger(4, std::stringstream()
              << "@" << _prefix(lod)
              << "write(lod=" << lod << ", pos=" << fmt(pos)
              << ", " << data->toString() << ")");
    // TODO-Low more const-correctness of the main API.
    auto unconst_data = std::const_pointer_cast<DataBuffer>(data);
    _accessor->writeRegion(unconst_data, pos, lod, /*is_storage=*/true, _lodcompressor);
    _report(data.get());
  }
}

} // namespace
