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

#pragma once

#include "../api.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

// Note: Might move this inside the main OpenZGY API if others need it.

namespace OpenZGY { namespace Tools {
#if 0
}}
#endif

/**
 * Base class for ZgyReaderRelay and ZgyWriterRelay.
 * Implements IZgyMeta and IZgyTools, forwarding all to the relay.
 */
class ZgyReadWriteRelay : public virtual IZgyTools
{
private:
  std::shared_ptr<IZgyTools> relay_;

protected:
  IZgyTools& relay() { return *relay_; }
  const IZgyTools& relay() const { return *relay_; }

public:
  typedef IZgyReader::int8_t       int8_t;
  typedef IZgyReader::int16_t      int16_t;
  typedef IZgyReader::int32_t      int32_t;
  typedef IZgyReader::int64_t      int64_t;
  typedef IZgyReader::float32_t    float32_t;
  typedef IZgyReader::float64_t    float64_t;
  typedef IZgyReader::size3i_t     size3i_t;
  typedef IZgyReader::corners_t    corners_t;
  typedef IZgyReader::rawdata_t    rawdata_t;
  typedef IZgyReader::compressor_t compressor_t;

public:
  explicit ZgyReadWriteRelay(const std::shared_ptr<IZgyTools>& relay)
    : relay_(relay)
  {
    if (!relay_)
      throw std::runtime_error("ZgyReadWriteRelay: empty relay");
  }

  virtual ~ZgyReadWriteRelay()
  {
  }

  ZgyReadWriteRelay(const ZgyReadWriteRelay&) = delete;
  ZgyReadWriteRelay(const ZgyReadWriteRelay&&) = delete;
  ZgyReadWriteRelay& operator=(const ZgyReadWriteRelay&) = delete;
  ZgyReadWriteRelay& operator=(const ZgyReadWriteRelay&&) = delete;

 public:

  // FUNCTIONS FROM IZgyMeta();

  size3i_t size() const override
  {
    return relay().size();
  }

  SampleDataType datatype() const override
  {
    return relay().datatype();
  }

  std::array<float32_t,2> datarange() const override
  {
    return relay().datarange();
  }

  std::array<float32_t,2> raw_datarange() const override
  {
    return relay().raw_datarange();
  }

  UnitDimension zunitdim() const override
  {
    return relay().zunitdim();
  }

  UnitDimension hunitdim() const override
  {
    return relay().hunitdim();
  }

  std::string zunitname() const override
  {
    return relay().zunitname();
  }

  std::string hunitname() const override
  {
    return relay().hunitname();
  }

  float64_t zunitfactor() const override
  {
    return relay().zunitfactor();
  }

  float64_t hunitfactor() const override
  {
    return relay().hunitfactor();
  }

  float32_t zstart() const override
  {
    return relay().zstart();
  }

  float32_t zinc() const override
  {
    return relay().zinc();
  }

  std::array<float32_t,2> annotstart() const override
  {
    return relay().annotstart();
  }

  std::array<float32_t,2> annotinc() const override
  {
    return relay().annotinc();
  }

  const corners_t corners() const override
  {
    return relay().corners();
  }

  const corners_t indexcorners() const override
  {
    return relay().indexcorners();
  }

  const corners_t annotcorners() const override
  {
    return relay().annotcorners();
  }

  size3i_t bricksize() const override
  {
    return relay().bricksize();
  }

  std::vector<size3i_t> brickcount() const override
  {
    return relay().brickcount();
  }

  int32_t nlods() const override
  {
    return relay().nlods();
  }

  // Currently not needed by any client.
  //std::string dataid() const override
  //{
  //  return relay().dataid();
  //}

  std::string verid()  const override
  {
    return relay().verid();
  }

  // Currently not needed by any client.
  //std::string previd() const override
  //{
  //  return relay().previd();
  //}

  void dump(std::ostream& os) const override
  {
    relay().dump(os);
  }

  SampleStatistics statistics() const override
  {
    return relay().statistics();
  }

  SampleHistogram histogram() const override
  {
    return relay().histogram();
  }

  std::shared_ptr<const FileStatistics> filestats() const override
  {
    return relay().filestats();
  }

  // FUNCTIONS FROM IZgyTools

  void transform(const corners_t& from, const corners_t& to, std::vector<std::array<float64_t,2>>& data) const override
  {
    relay().transform(from, to, data);
  }

  std::array<float64_t,2> transform1(const corners_t& from, const corners_t& to, const std::array<float64_t,2>& point) const override
  {
    return relay().transform1(from, to, point);
  }

  // The converers need to be re-implemented and not forwarded,
  // because they needto pick up the correct corner points.

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
 * Wrapper for IZgyReader that forwards all calls to the underlying
 * wrapper. It needs to be subclassed to be of any use. A derived
 * class can override one or more methods in the API, without seeing a
 * lot of clutter from methods it doesn't need to change.
 *
 * If a derived class ends up overriding most of the methods, then
 * this helper should probably not be used.
 */

class ZgyReaderRelay : public ZgyReadWriteRelay, virtual public IZgyReader
{
public:
  explicit ZgyReaderRelay(const std::shared_ptr<IZgyReader>& relay)
    : ZgyReadWriteRelay(relay)
  {
  }

  virtual ~ZgyReaderRelay()
  {
  }
  ZgyReaderRelay(const ZgyReaderRelay&) = delete;
  ZgyReaderRelay(const ZgyReaderRelay&&) = delete;
  ZgyReaderRelay& operator=(const ZgyReaderRelay&) = delete;
  ZgyReaderRelay& operator=(const ZgyReaderRelay&&) = delete;

protected:

  IZgyReader& relay() {
    return dynamic_cast<IZgyReader&>(ZgyReadWriteRelay::relay());
  }

  const IZgyReader& relay() const {
    return dynamic_cast<const IZgyReader&>(ZgyReadWriteRelay::relay());
  }

public:
  // FUNCTIONS FROM IZgyReader

  void read(const size3i_t& start, const size3i_t& size, float* data, int lod) const override
  {
    relay().read(start, size, data, lod);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod) const override
  {
    relay().read(start, size, data, lod);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int8_t* data, int lod) const override
  {
    relay().read(start, size, data, lod);
  }

  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, int lod, bool as_float = true) const override
  {
    return relay().readconst(start, size, lod, as_float);
  }

  void close() override
  {
    relay().close();
  }
};

/**
 * Wrapper for IZgyWriter that forwards all calls to the underlying
 * wrapper. It needs to be subclassed to be of any use. A derived
 * class can override one or more methods in the API, without seeing a
 * lot of clutter from methods it doesn't need to change.
 *
 * If a derived class ends up overriding most of the methods, then
 * this helper should probably not be used.
 */

class ZgyWriterRelay : public ZgyReadWriteRelay, virtual public IZgyWriter
{
public:
  explicit ZgyWriterRelay(const std::shared_ptr<IZgyWriter>& relay)
    : ZgyReadWriteRelay(relay)
  {
  }

  virtual ~ZgyWriterRelay()
  {
  }

  ZgyWriterRelay(const ZgyWriterRelay&) = delete;
  ZgyWriterRelay(const ZgyWriterRelay&&) = delete;
  ZgyWriterRelay& operator=(const ZgyWriterRelay&) = delete;
  ZgyWriterRelay& operator=(const ZgyWriterRelay&&) = delete;

protected:

  IZgyWriter& relay() {
    return dynamic_cast<IZgyWriter&>(ZgyReadWriteRelay::relay());
  }

  const IZgyWriter& relay() const {
    return dynamic_cast<const IZgyWriter&>(ZgyReadWriteRelay::relay());
  }

public:
  // FUNCTIONS FROM IZgyWriter

  void read(const size3i_t& start, const size3i_t& size, float* data) const override
  {
    relay().read(start, size, data);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const override
  {
    relay().read(start, size, data);
  }

  void read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const override
  {
    relay().read(start, size, data);
  }

  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float = true) const override
  {
    return relay().readconst(start, size, as_float);
  }

  void write(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    relay().write(start, size, data);
  }

  void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data) override
  {
    relay().write(start, size, data);
  }

  void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    relay().write(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const float* data) override
  {
    relay().writeconst(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) override
  {
    relay().writeconst(start, size, data);
  }

  void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data) override
  {
    relay().writeconst(start, size, data);
  }

  void finalize(
       const std::vector<DecimationType>& decimation,
       const std::function<bool(std::int64_t,std::int64_t)>& progress,
       FinalizeAction action,
       bool force) override
  {
    relay().finalize(decimation, progress, action, force);
  }

  void close_incomplete() override
  {
    relay().close_incomplete();
  }

  void close() override
  {
    relay().close();
  }

  bool errorflag() const override
  {
    return relay().errorflag();
  }

  void set_errorflag(bool flag) override
  {
    relay().set_errorflag(flag);
  }
};

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
