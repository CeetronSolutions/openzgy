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

#include "api.h"
#include "impl/fancy_timers.h"

namespace OpenZGY { namespace Impl {
#if 0
}}
#endif

/**
 * Wrapper for IZgyWriter that fully synchronizes all access not known to already be thread safe. This is done as a convenience to application developers. IZgyWriter is documented as being single threaded, but by introducing this wrapper the application may use multiple threads and let the library worry about locking.
 *
 * The main downside is the overhead introduced by locking std::mutex for most calls. All applications will see this cost, even those that play nice and use only one thread. If the overhead becomes too high then add a no_lock argument to IZgyWriter::open(). Otherwise just leave it in place in the name of defensive coding. A.k.a. don't expect programmers to read the manual.
 *
 * Only the application should call this class explicitly. Even when one method is implemented in terms of other api methods this will be done one layer below. So those methods should also be synchronized.
 *
 * A single non-recursive mutex is used. If this causes a deadlock them most likely somebody managed to call this class from a lower level. If this happens, please fix the bug and do not simply make the mutex recursive.
 *
 * Users can also cause a deadlock if doing something silly such as invoking api methods from inside a progress callback. If this happens, please fix the user.
 */

class ZgySafeWriter : public IZgyWriter
{
 private:
  std::shared_ptr<IZgyWriter> writer_;
  mutable std::mutex mutex_;
  mutable std::shared_ptr<InternalZGY::SummaryPrintingTimerEx> mtimer_;
  mutable std::shared_ptr<InternalZGY::SummaryPrintingTimerEx> ctimer_;
  mutable std::shared_ptr<InternalZGY::SummaryPrintingTimerEx> ftimer_;
  // The lifetime of this class, not useful unless busy all the time.
  mutable std::shared_ptr<InternalZGY::SummaryPrintingTimerEx> atimer_;
  mutable std::shared_ptr<InternalZGY::SimpleTimerEx> simple_atimer_;

public:
  typedef IZgyWriter::int8_t       int8_t;
  typedef IZgyWriter::int16_t      int16_t;
  typedef IZgyWriter::int32_t      int32_t;
  typedef IZgyWriter::int64_t      int64_t;
  typedef IZgyWriter::float32_t    float32_t;
  typedef IZgyWriter::float64_t    float64_t;
  typedef IZgyWriter::size3i_t     size3i_t;
  typedef IZgyWriter::corners_t    corners_t;
  typedef IZgyWriter::rawdata_t    rawdata_t;
  typedef IZgyWriter::compressor_t compressor_t;

public:
  ZgySafeWriter(const std::shared_ptr<IZgyWriter>& writer)
    : writer_(writer), mutex_()
  {
    mtimer_.reset(new InternalZGY::SummaryPrintingTimerEx("ZgyWriter.mutex"));
    ctimer_.reset(new InternalZGY::SummaryPrintingTimerEx("ZgyWriter.close"));
    ftimer_.reset(new InternalZGY::SummaryPrintingTimerEx("ZgyWriter.final"));
    atimer_.reset(new InternalZGY::SummaryPrintingTimerEx("ZgyWriter.inuse"));
    simple_atimer_.reset(new InternalZGY::SimpleTimerEx(*atimer_));
  }
  ZgySafeWriter(const ZgySafeWriter&) = delete;
  ZgySafeWriter& operator=(const ZgySafeWriter&) = delete;

 public:
  // FUNCTIONS FROM IZgyMeta();
  size3i_t size() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->size();
  }

  SampleDataType datatype() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->datatype();
  }

  std::array<float32_t,2> datarange() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->datarange();
  }

  std::array<float32_t,2> raw_datarange() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->raw_datarange();
  }

  UnitDimension zunitdim() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->zunitdim();
  }

  UnitDimension hunitdim() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->hunitdim();
  }

  std::string zunitname() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->zunitname();
  }

  std::string hunitname() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->hunitname();
  }

  float64_t zunitfactor() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->zunitfactor();
  }

  float64_t hunitfactor() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->hunitfactor();
  }

  float32_t zstart() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->zstart();
  }

  float32_t zinc() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->zinc();
  }

  std::array<float32_t,2> annotstart() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->annotstart();
  }

  std::array<float32_t,2> annotinc() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->annotinc();
  }

  const corners_t corners() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->corners();
  }

  const corners_t indexcorners() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->indexcorners();
  }

  const corners_t annotcorners() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->annotcorners();
  }

  size3i_t bricksize() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->bricksize();
  }

  std::vector<size3i_t> brickcount() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->brickcount();
  }

  int32_t nlods() const override
  {
    // Not allowed to change after file is created.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->nlods();
  }

  // Currently not needed by any client.
  //virtual std::string dataid() const override
  //{
  //  // Not allowed to change after file is opened.
  //  // std::lock_guard<std::mutex> lk(mutex_);
  //  return writer_->dataid();
  //}

  std::string verid()  const override
  {
    // Not allowed to change after file is opened.
    // std::lock_guard<std::mutex> lk(mutex_);
    return writer_->verid();
  }

  // Currently not needed by any client.
  //virtual std::string previd() const override
  //{
  //  // Not allowed to change after file is opened.
  //  // std::lock_guard<std::mutex> lk(mutex_);
  //  return writer_->previd();
  //}

  void dump(std::ostream& os) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->dump(os);
  }

  SampleStatistics statistics() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->statistics();
  }

  SampleHistogram histogram() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->histogram();
  }

  std::shared_ptr<const FileStatistics> filestats() const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->filestats();
  }

  // FUNCTIONS FROM IZgyTools
  void transform(const corners_t& from, const corners_t& to, std::vector<std::array<float64_t,2>>& data) const override
  {
    // Semantically static and is thread safe.
    // Virtual only to satisfy the interface.
    //std::lock_guard<std::mutex> lk(mutex_);
    writer_->transform(from, to, data);
  }

  std::array<float64_t,2> transform1(const corners_t& from, const corners_t& to, const std::array<float64_t,2>& data) const override
  {
    // Semantically static and is thread safe.
    // Virtual only to satisfy the interface.
    //std::lock_guard<std::mutex> lk(mutex_);
    return writer_->transform1(from, to, data);
  }

  std::array<float64_t,2> annotToIndex(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->annotToIndex(data);
  }

  std::array<float64_t,2> annotToWorld(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->annotToWorld(data);
  }

  std::array<float64_t,2> indexToAnnot(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->indexToAnnot(data);
  }

  std::array<float64_t,2> indexToWorld(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->indexToWorld(data);
  }

  std::array<float64_t,2> worldToAnnot(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->worldToAnnot(data);
  }

  std::array<float64_t,2> worldToIndex(const std::array<float64_t,2>& data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->worldToIndex(data);
  }

  // FUNCTIONS FROM IZgyWriter

  void read(const size3i_t& start, const size3i_t& size, float* data) const override
  {
    // TODO-Low: Use std::shared_mutex available in C++17
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->read(start, size, data);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->read(start, size, data);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->read(start, size, data);
  }

  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float = true) const override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    return writer_->readconst(start, size, as_float);
  }

  void write(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->write(start, size, data);
  }

  void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->write(start, size, data);
  }

  void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->write(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->writeconst(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->writeconst(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->writeconst(start, size, data);
  }

  void finalize(
       const std::vector<DecimationType>& decimation,
       const std::function<bool(std::int64_t,std::int64_t)>& progress,
       FinalizeAction action,
       bool force) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    InternalZGY::SimpleTimerEx ff(*ftimer_);
    writer_->finalize(decimation, progress, action, force);
  }

  void close_incomplete() override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    InternalZGY::SimpleTimerEx cc(*ctimer_);
    writer_->close_incomplete();
  }

  void close() override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    InternalZGY::SimpleTimerEx cc(*ctimer_);
    writer_->close();
  }

  bool errorflag() const override
  {
    // No locking, implemented using atomics.
    //std::lock_guard<std::mutex> lk(mutex_);
    return writer_->errorflag();
  }

  void set_errorflag(bool flag) override
  {
    InternalZGY::SimpleTimerEx mm(*mtimer_);
    std::lock_guard<std::mutex> lk(mutex_);
    mm.stop();
    writer_->set_errorflag(flag);
  }
};

}} // namespace
