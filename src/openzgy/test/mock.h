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

#include "../api.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace Test {
#if 0
}
#endif

using namespace OpenZGY;

class ZgyMetaMock : virtual public IZgyMeta, virtual public IZgyTools
{
  size3i_t size_;
  SampleDataType datatype_;
public:
  ~ZgyMetaMock();
  ZgyMetaMock(const size3i_t& size, SampleDataType dt)
    : size_(size)
    , datatype_(dt)
  {
  }
  // Functions from IZgyMeta
  virtual size3i_t       size()           const {return size_;}
  virtual SampleDataType datatype()       const {return datatype_;}
  virtual std::array<float32_t,2> datarange() const {return std::array<float32_t,2>{-1,1};}
  virtual std::array<float32_t,2> raw_datarange() const {return std::array<float32_t,2>{-1,1};}
  virtual UnitDimension  zunitdim()       const {return UnitDimension::unknown;}
  virtual UnitDimension  hunitdim()       const {return UnitDimension::unknown;}
  virtual std::string    zunitname()      const {return std::string();}
  virtual std::string    hunitname()      const {return std::string();}
  virtual float64_t      zunitfactor()    const {return 1;}
  virtual float64_t      hunitfactor()    const {return 1;}
  virtual float32_t      zstart()         const {return 0;}
  virtual float32_t      zinc()           const {return 4;}
  virtual std::array<float32_t,2> annotstart() const {return std::array<float32_t,2>{1,1};}
  virtual std::array<float32_t,2> annotinc()   const {return std::array<float32_t,2>{1,1};}
  virtual const corners_t corners()      const {static corners_t result{{{0,0},{0,0},{0,0},{0,0}}}; return result;}
  virtual const corners_t indexcorners() const {static corners_t result{{{0,0},{0,0},{0,0},{0,0}}}; return result;}
  virtual const corners_t annotcorners() const {static corners_t result{{{0,0},{0,0},{0,0},{0,0}}}; return result;}
  virtual size3i_t       bricksize()      const {return size3i_t{64,64,64};}
  virtual std::vector<size3i_t> brickcount() const {throw std::runtime_error("brickcount() has not been mocked");}
  virtual int32_t        nlods()          const {throw std::runtime_error("nlods() has not been mocked");}
  virtual std::string    dataid()         const {return "0000000a-000b-000c-000d-333333333333";}
  virtual std::string    verid()          const {return "0000000a-000b-000c-000d-222222222222";}
  virtual std::string    previd()         const {return "0000000a-000b-000c-000d-111111111111";}
  virtual void           meta()           const {throw std::runtime_error("meta() has not been mocked");}
  virtual int32_t        numthreads()     const {throw std::runtime_error("numthreads() has not been mocked");}
  virtual void           set_numthreads(int32_t) {throw std::runtime_error("set_numthreads() has not been mocked");}
  virtual void           dump(std::ostream&) const {throw std::runtime_error("dump() has not been mocked");}
  virtual SampleStatistics statistics()   const {throw std::runtime_error("statistics() has not been mocked");}
  virtual SampleHistogram  histogram()    const {throw std::runtime_error("histogram() has not been mocked");}
  virtual std::shared_ptr<const FileStatistics> filestats()    const {throw std::runtime_error("filestats() has not been mocked");}
  // Functions from IZgyTools
  virtual void transform(const corners_t& from, const corners_t& to, std::vector<std::array<float64_t,2>>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> transform1(const corners_t& from, const corners_t& to, const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> annotToIndex(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> annotToWorld(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> indexToAnnot(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> indexToWorld(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> worldToAnnot(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
  virtual std::array<float64_t,2> worldToIndex(const std::array<float64_t,2>&) const {throw std::runtime_error("coord handling has not been mocked");}
};

class ZgyReaderMock : public ZgyMetaMock, virtual public IZgyReader
{
public:
  virtual ~ZgyReaderMock();
  ZgyReaderMock(const size3i_t& size, SampleDataType dt)
    : ZgyMetaMock(size, dt)
  {
  }
  virtual void read(const size3i_t& start, const size3i_t& size, float* data, int lod = 0) const {
    const std::int64_t n = size[0]*size[1]*size[2];
    float value{0};
    for (std::int64_t ii=0; ii<n; ++ii, value += 1)
      data[ii] = value;
  }
  virtual void read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod = 0) const {
    const std::int64_t n = size[0]*size[1]*size[2];
    std::int32_t value{0};
    for (std::int64_t ii=0; ii<n; ++ii, value = (value+1) % 32768)
      data[ii] = (std::int16_t)value;
  }
  virtual void read(const size3i_t& start, const size3i_t& size, std::int8_t* data, int lod = 0) const {
    const std::int64_t n = size[0]*size[1]*size[2];
    std::int32_t value{0};
    for (std::int64_t ii=0; ii<n; ++ii, value = (value+1) % 128)
      data[ii] = (std::int8_t)value;
  }
  virtual std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, int lod = 0, bool as_float = true) const {
    return std::make_pair(false,0);
  }
  virtual void close()
  {
  }
  static std::shared_ptr<IZgyReader> mock(const size3i_t& size, SampleDataType dt) {
    return std::shared_ptr<IZgyReader>(new ZgyReaderMock(size, dt));
  }
  static std::shared_ptr<IZgyReader> mock(const std::array<std::int64_t,4>& size) {
    const size3i_t sz{size[0],size[1],size[2]};
    switch (size[3]) {
    case 1: return mock(sz, SampleDataType::int8);
    case 2: return mock(sz, SampleDataType::int16);
    case 4: return mock(sz, SampleDataType::float32);
    default:
      throw std::runtime_error("Bytes per sample must be 1, 2, or 4,");
    }
  }
};

class ZgyWriterMock : public ZgyMetaMock, virtual public IZgyWriter
{
public:
  ZgyWriterMock(const size3i_t& size, SampleDataType dt)
    : ZgyMetaMock(size, dt)
  {
  }
  virtual ~ZgyWriterMock();
  virtual void read(const size3i_t& start, const size3i_t& size, float* data) const {}
  virtual void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const {}
  virtual void read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const {}
  virtual std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float = true) const {
    return std::make_pair(false,0);
  }
  virtual void write(const size3i_t& start, const size3i_t& size, const float* data) {}
  virtual void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data) {}
  virtual void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data) {}
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const float* data) {}
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) {}
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data) {}
  virtual void finalize(
               const std::vector<DecimationType>& decimation = std::vector<DecimationType>(),
               const std::function<bool(std::int64_t,std::int64_t)>& progress = nullptr,
               FinalizeAction action = FinalizeAction::BuildDefault,
               bool force = false){}
  virtual void close_incomplete() {}
  virtual void close() {}
  virtual bool errorflag() const {return false;}
  virtual void set_errorflag(bool) {}
  static std::shared_ptr<IZgyWriter> mock(const size3i_t& size, SampleDataType dt) {
    return std::shared_ptr<IZgyWriter>(new ZgyWriterMock(size, dt));
  }
  static std::shared_ptr<IZgyWriter> mock(const ZgyWriterArgs& args) {
    return mock(args._size, args._datatype);
  }
};

}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
