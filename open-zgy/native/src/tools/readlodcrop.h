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

#include "readwriterelay.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY { namespace Tools {
#if 0
}}
#endif

/**
 * Virtual cropped and decimated seismic cube for reading.
 *
 * Intercept all calls to the public API that needs to return a
 * different result compared to the input.
 *
 * The idea is that the features of cropping the area and selecting the lod level
 * to read can be added to the existing zgycopyc tool with
 * just a few lines of code. So it won't bloat that tool even more.
 */
class ZgyReadLodCrop : public ZgyReaderRelay
{
public:
  ZgyReadLodCrop(
       const std::shared_ptr<IZgyReader>& relay_in,
       int startlod,
       const size3i_t& cropstart,
       const size3i_t& cropsize);
  virtual ~ZgyReadLodCrop();

public:

  // FUNCTIONS FROM IZgyReader:

  void read(
       const size3i_t& start, const size3i_t& size,
       float* data, int lod) const override
  {
    size3i_t realstart
      {start[0] + cropstart_[0],
       start[1] + cropstart_[1],
       start[2] + cropstart_[2]};
    // The API already expects the caller to provide start and size relative
    // to the LOD being requested. And for the caller, the actual size of
    // the survey as lod0 is what the real cube has at startlod_.
    relay().read(realstart, size, data, lod + this->startlod_);
  }

  void read(
       const size3i_t& start, const size3i_t& size,
       std::int16_t* data, int lod = 0) const override
  {
    size3i_t realstart
      {start[0] + cropstart_[0],
       start[1] + cropstart_[1],
       start[2] + cropstart_[2]};
    relay().read(realstart, size, data, lod + this->startlod_);
  }

  void read(
       const size3i_t& start, const size3i_t& size,
       std::int8_t* data, int lod = 0) const override
  {
    size3i_t realstart
      {start[0] + cropstart_[0],
       start[1] + cropstart_[1],
       start[2] + cropstart_[2]};
    relay().read(realstart, size, data, lod + this->startlod_);
  }

  std::pair<bool,double> readconst(
       const size3i_t& start, const size3i_t& size,
       int lod, bool as_float) const override
  {
    size3i_t realstart
      {start[0] + cropstart_[0],
       start[1] + cropstart_[1],
       start[2] + cropstart_[2]};
    return relay().readconst(realstart, size, lod + this->startlod_, as_float);
  }

  // FUNCTIONS FROM IZgyMeta();

  size3i_t              size()         const override { return cropsize_; }
  float32_t             zstart()       const override { return zstart_;}
  float32_t             zinc()         const override { return zinc_;}
  std::array<float32_t,2> annotstart() const override { return annotstart_; }
  std::array<float32_t,2> annotinc()   const override { return annotinc_; }
  const corners_t       corners()      const override { return worldcorners_; }
  const corners_t       indexcorners() const override { return indexcorners_; }
  const corners_t       annotcorners() const override { return annotcorners_; }
  std::vector<size3i_t> brickcount()   const override { return brickcount_; }
  int32_t               nlods()        const override { return nlods_; }
  std::string           verid()        const override { return verid_; }

  // Statistics(), histogram(), and filestats() will not be upated.
  // That would be prohibitivly expensive. Returning the real data
  // is probably better than throwing an exception. as the information
  // will be mostly correct.

private:
  void checkAnnotation(double ipos, double jpos);

private:
  int                     startlod_;
  std::int64_t            lodfactor_;
  size3i_t                cropstart_;
  size3i_t                cropsize_;
  float32_t               zstart_;
  float32_t               zinc_;
  std::array<float32_t,2> annotstart_;
  std::array<float32_t,2> annotinc_;
  corners_t               worldcorners_;
  corners_t               indexcorners_;
  corners_t               annotcorners_;
  std::vector<size3i_t>   brickcount_;
  int32_t                 nlods_;
  std::string             verid_;
};

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
