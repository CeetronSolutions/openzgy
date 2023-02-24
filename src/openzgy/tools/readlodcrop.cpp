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

#include "readlodcrop.h"
#include "readwriterelay.h"
#include "../impl/guid.h"

#include <tuple>
#include <iostream>
#include <iomanip>
#include <cmath>

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

namespace {
  static std::string format(const std::array<double,2>& in)
  {
    std::stringstream ss;
    ss << "(" << (float)in[0] << ", " << (float)in[1] << ")";
    return ss.str();
  };

  static std::string format(const std::array<std::int64_t,3>& in)
  {
    std::stringstream ss;
    ss << in[0] << "x" << in[1] << "x" << in[2];
    return ss.str();
  };
}

// From zgydumpc, TODO remove when done debugging.
void
dump_basic(OpenZGY::IZgyReader* r, const std::string& filename, std::ostream& os)
{
  using namespace OpenZGY;
  using namespace OpenZGY::Formatters;

  auto oldprec = os.precision();
  auto oldflags = os.flags();

  auto ocp = [r](int ix) {
               // Original formats %5d %9g %11.2f
               std::stringstream ss;
               ss << "["
                  << std::setw(5) << r->indexcorners()[ix][0] << ", "
                  << std::setw(5) << r->indexcorners()[ix][1] << "] {"
                  << std::setw(9) << r->annotcorners()[ix][0] << ", "
                  << std::setw(9) << r->annotcorners()[ix][1] << "} ("
                  << std::fixed << std::setprecision(2)
                  << std::setw(11) << r->corners()[ix][0] << ", "
                  << std::setw(11) << r->corners()[ix][1] << ")";
               return ss.str();
             };

#if 0
  const SampleStatistics stat = r->statistics();
  const SampleHistogram  hist = r->histogram();
  std::shared_ptr<const FileStatistics> filestats = r->filestats();
#endif

  // Note that file size in bytes, ZGY version, and projection system
  // are not available in the API. They might not be useful anyway.

  os << "" // << "File name                      = '" << filename << "'\n"
#if 0
     << "File size (bytes)              = " << filestats->fileSize() << "\n"
     << "File format and version        = " << r->datatype() << " ZGY version " << filestats->fileVersion() << "\n"
     //<< "Data identifier                = " << r->dataid() << "\n"
     << "Current data Version           = " << r->verid() << "\n"
     //<< "Previous data version          = " << r->previd() << "\n"
     << "Brick size I,J,K               = " << "(" << r->bricksize()[0] << ", " << r->bricksize()[1] << ", " << r->bricksize()[2] << ")\n"
     << "Number of bricks I,J,K         = " << "(" << r->brickcount()[0][0] << ", " << r->brickcount()[0][1] << ", " << r->brickcount()[0][2] << ")\n"
     << "Number of LODs                 = " << r->nlods() << "\n"
     << "Coding range min/max           = " << std::setprecision(6) << r->datarange()[0] << " " << r->datarange()[1] << " (raw: " << r->raw_datarange()[0] << " " << r->raw_datarange()[1] << ") " << r->size()[0] * r->size()[1] * r->size()[2] << "\n"
     << "Statistical min/max/count      = " << std::setprecision(6) << stat.min << " " << stat.max << " " << stat.cnt << "\n"
     << "Histogram range min/max/count  = " << std::setprecision(6) << hist.minvalue << " " << hist.maxvalue << " " << hist.samplecount << "\n"
     << "Inline start/increment/count   = " << r->annotstart()[0] << " " << r->annotinc()[0] << " " << r->size()[0] << "\n"
     << "Xline  start/increment/count   = " << r->annotstart()[1] << " " << r->annotinc()[1] << " " << r->size()[1] << "\n"
     << "Sample start/increment/count   = " << r->zstart() << " " << r->zinc() << " " << r->size()[2] << "\n"
     << "Horizontal projection system   = " << "?\n" // {r._accessor._metadata._ih._hprjsys}
     << "Horizontal dim/factor/name     = " << r->hunitdim() << " " << r->hunitfactor() << " '" << r->hunitname() << "'\n"
     << "Vertical dim/factor/name       = " << r->zunitdim() << " " << r->zunitfactor() << " '" << r->zunitname() << "'\n"
#endif
     << "Ordered Corner Points Legend   = " << "[  <i>,   <j>] { <inline>,   <xline>} (  <easting>,  <northing>)\n"
     << "Ordered Corner Point 1         = " << ocp(0) << "\n"
     << "Ordered Corner Point 2         = " << ocp(1) << "\n"
     << "Ordered Corner Point 3         = " << ocp(2) << "\n"
     << "Ordered Corner Point 4         = " << ocp(3) << "\n"
     << std::flush;

  os.flags(oldflags);
  os.precision(oldprec);
}

/**
 * Shamelessly duplicated from meta.cpp
 *
 * Compute the size of the file in bricks, for all LOD levels
 * and all 3 dimensions. Indirectly also compute the number
 * of LOD levels, since by definition the last LOD level
 * is the one that has size (1, 1, 1) i.e. the entire cube
 * is decimated enought to fit inside a single brick.
 */
static std::vector<std::array<std::int64_t,3>>
calcLodSizes(const std::array<std::int64_t,3>& size_in,
	     const std::array<std::int64_t,3>& bricksize)
{
  std::array<std::int64_t,3> size = size_in;
  std::vector<std::array<std::int64_t,3>> result;
  if (size[0] < 1 || size[1] < 1 || size[2] < 1) {
    // Empty file is the only case where the highest lod isn't one brick.
    // There will be one lod and its size will be zero bricks.
    result.push_back(std::array<std::int64_t,3>{0,0,0});
  }
  else {
    // number of full and partial bricks:
    for (int ii=0; ii<3; ++ii)
      size[ii] = (size[ii] + bricksize[ii] - 1) / bricksize[ii];
    result.push_back(size);
    // Continue to add new levels until we get a level with just one brick.
    while (size[0] > 1 || size[1] > 1 || size[2] > 1) {
      for (int ii=0; ii<3; ++ii)
        size[ii] = (size[ii] + 1) / 2;
      result.push_back(size);
    }
  }
  return result;
}

/**
 * Adjust the cropping specified by the user.
 *
 * Clip invalid values.
 * If just one of --cropstart of --cropsize is specified,
 * default the other to "center".
 *
 * \param cropstart As requested by the user
 * \param cropsize  As requested by the user
 * \param size      Size of underlying cube, adjusted for lod level.
 */
static std::tuple<std::int64_t, std::int64_t>
adjustCropping(std::int64_t cropstart, std::int64_t cropsize, std::int64_t size)
{
  if (cropstart >= 0 && cropsize > 0) {
    // Both specified.
    cropstart = std::min(cropstart, size - 1);
    cropsize = std::min(cropsize, size - cropstart);
  }
  else if (cropstart < 0 && cropsize > 0) {
    // Only --cropsize. Crop leaving the center.
    cropsize = std::min(cropsize, size);
    cropstart = (size - cropsize) / 2;
  }
  else if (cropstart >= 0 && cropsize <= 0) {
    // Only --cropstart. Crop leaving the center.
    // Start must be less than half the survey size,
    // otherwise the result will be from start to end of cube.
    cropstart = std::min(cropstart, size - 1);
    cropsize = size - 2 * cropstart;
    if (cropsize < 1)
      cropsize = size - cropstart;
  }
  else {
    // Neither --cropstart nor --cropsize, or something seriously wrong.
    // No cropping in this dimension.
    cropsize = size;
    cropstart = 0;
  }
  return std::make_tuple(cropstart, cropsize);
}

ZgyReadLodCrop::~ZgyReadLodCrop()
{
}

ZgyReadLodCrop::ZgyReadLodCrop(
     const std::shared_ptr<IZgyReader>& relay_in,
     int startlod,
     const size3i_t& cropstart,
     const size3i_t& cropsize)
  : ZgyReaderRelay(relay_in)
  , startlod_(startlod)
  , lodfactor_(std::int64_t(1) << startlod)
  , cropstart_(cropstart)
  , cropsize_(cropsize)
{
  // Note: Apply lod first, so cropstart and cropsize is relative to startlod.

  for (int dim = 0; dim < 3; ++dim) {
    const std::int64_t real_size_at_lod =
      (relay().size()[dim] + lodfactor_ - 1) / lodfactor_;
    std::tie(cropstart_[dim], cropsize_[dim]) =
      adjustCropping(cropstart_[dim], cropsize_[dim], real_size_at_lod);
  }
  std::cerr << "fullsize "     << format(relay().size())
            << " --lod "       << startlod_
            << " --cropstart " << format(cropstart_)
            << " --cropsize "  << format(cropsize_)
            << std::endl;

  // Annotation numbers both in XY and in Z are not relative to
  // a particular lod, but the increments increase at higher lod.
  zinc_   = relay().zinc() * lodfactor_;
  annotinc_[0] = relay().annotinc()[0] * lodfactor_;
  annotinc_[1] = relay().annotinc()[1] * lodfactor_;

  zstart_ = relay().zstart() + zinc_ * cropstart_[2];
  annotstart_[0] = relay().annotstart()[0] + cropstart_[0] * annotinc_[0];
  annotstart_[1] = relay().annotstart()[1] + cropstart_[1] * annotinc_[1];

  // Note that annotToIndex() and its 5 friends are implemented in
  // term of the virtual corners(), annotcorners(), and indexcorners().
  // Once they are set correctly, all 6 transforms ought to work.
  // At this point in the code they will *not* work.

  const double beg[2] {(double)cropstart_[0], (double)cropstart_[1]};
  const double siz[2] {(double)cropsize_[0],  (double)cropsize_[1]};
  const double end[2] {beg[0] + siz[0] - 1, beg[1] + siz[1] - 1};
  // Ordinal corners relative to relay() and the real lod0.
  this->indexcorners_ = corners_t
    {{{beg[0]*lodfactor_, beg[1]*lodfactor_},
      {end[0]*lodfactor_, beg[1]*lodfactor_},
      {beg[0]*lodfactor_, end[1]*lodfactor_},
      {end[0]*lodfactor_, end[1]*lodfactor_}}};
  for (int ii=0; ii<4; ++ii)
    this->annotcorners_[ii] = relay().indexToAnnot(this->indexcorners_[ii]);
  for (int ii=0; ii<4; ++ii)
    this->worldcorners_[ii] = relay().indexToWorld(this->indexcorners_[ii]);
  // Annotation and world are the same for all lods, while ordinals are not.
  this->indexcorners_ = corners_t
    {{{0, 0},
      {(double)cropsize_[0]-1,                      0},
      {                     0, (double)cropsize_[1]-1},
      {(double)cropsize_[0]-1, (double)cropsize_[1]-1}}};

  // brickcount a.k.a. lodsizes, and nlods. Re-calculate from size.
  this->brickcount_ = calcLodSizes(cropsize_, relay().bricksize());
  this->nlods_ = static_cast<std::int32_t>(this->brickcount_.size());

  // Make the file appear different than its source. As if anybody cares.
  this->verid_ = InternalZGY::GUID::makeGUID();

  checkAnnotation(0, 0);
  checkAnnotation(1, 1);
  checkAnnotation(5, 3);
}

/**
 * annotstart, annotinc are redundant with annotcorners.
 * The followig consistency check should technically be
 * in some unit test. But it is cheap enough to do here.
 * Caveat: If called from a contructor, virtual function
 * calls might not do what you expect.
 */
void
ZgyReadLodCrop::checkAnnotation(double ipos, double jpos)
{
  const std::array<float64_t,2> test
    {annotstart_[0] + ipos*annotinc_[0], annotstart_[1] + jpos*annotinc_[1]};
  const std::array<float64_t,2> expect{ipos, jpos};
  const std::array<float64_t,2> actual = annotToIndex(test);
  const std::array<float64_t,2> round  = indexToAnnot(actual);
  bool error = (std::abs(actual[0] - expect[0]) > 0.001 ||
                std::abs(actual[1] - expect[1]) > 0.001);
  if (error) {
    static const auto fmt = [](const std::array<float64_t,2>& in)
      {
        std::stringstream ss;
        ss << "(" << (float)in[0] << ", " << (float)in[1] << ")";
        return ss.str();
      };
    std::cerr << "ZgyReadLodCrop:"
              << " Annot " << fmt(test)
              << " index " << fmt(actual)
              << " expected " << fmt(expect)
              << " roundtrip " << fmt(round)
              << std::endl;
  }
  if (error) {
    throw std::runtime_error("ZgyReadLodCrop: Error computing annotation.");
  }
}

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
