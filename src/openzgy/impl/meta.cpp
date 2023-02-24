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

#include "meta.h"
#include "types.h"
#include "file.h"
#include "cornerpoints.h"
#include "structaccess.h"
#include "transform.h"
#include "guid.h"
#include "logger.h"
#include "file_smallcache.h"
#include "environment.h"
#include "../exception.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#define USE_NEW_LATTICE 0 // TODO-Low enable this after more testing.
#define TEST_NEW_LATTICE 0

namespace InternalZGY {
#if 0
}
#endif

// This is to bring several operator<< methods into scope. Note that
// placing the declaration in the global namespace doesn't seem to
// work. Placing it here gives it a narrower scope and (I speculate)
// somewhat higher precedence. This problem is only for the <<
// overload. A hypothetical InternalZGY::Formatters::foo() would be
// found just fine either way.
using namespace InternalZGY::Formatters;

/** file: meta.cpp
 *
 * TODO-Low: Move most of this code to a new meta_storage.cpp/h.
 * This file should only implement ZgyInternalMeta.
 *
 * Most of the existing meta.cpp deals with version-specific layout on
 * the ZGY file. All that code should be private and moved to the new
 * source file. Exports from meta_storage.h would be the factories.
 * One for each header type. The interface returned by each factory is
 * defined in the existing meta.h. That would not change.
 *
 * I cannot do this today because ZgyInternalMeta::initFromScratch()
 * that creates a new file does not use factories. It doesn't really
 * need to since only the latest version of the file layout can be
 * created. I could of course just move the entire initFromScratch()
 * into a single factory method in the new meta_storage.h. But then
 * there would be almost nothing left in this file and I would be back
 * where I started. The current situation actually isn't that bad
 * since all the version-specific classes are defined here in the .cpp
 * source. So even today they aren't visible outside.
 */
namespace {
  template<typename T, std::size_t  N>
  std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a)
  {
    os << "[";
    for (std::size_t ii=0; ii<N; ++ii)
      os << a[ii] << (ii == N-1 ? "" : ", ");
    os << "]";
    return os;
  }
}

namespace {
/**
 * Technically I should have mapped these format specific
 * codes to am enum defined in the impl layer. But since
 * the enum is really unlikely to change I'll just use casts
 * and expose the RawXXX enums to the layers above.
 * They will need to be remapped in the public api anyway.
 * RawGridDefinition and RawCoordType are even less of a
 * problem as those should not leave this layer.
 */
static RawDataType            DecodeDataType(std::uint8_t e)            { return static_cast<RawDataType>(e); }
static RawHorizontalDimension DecodeHorizontalDimension(std::uint8_t e) { return static_cast<RawHorizontalDimension>(e); }
static RawVerticalDimension   DecodeVerticalDimension(std::uint8_t e)   { return static_cast<RawVerticalDimension>(e); }
static RawCoordType           DecodeCoordType(std::uint8_t e)           { return static_cast<RawCoordType>(e); }
//static RawGridDefinition      DecodeGridDefinition(std::uint8_t e)      { return static_cast<RawGridDefinition>(e); }

/**
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
 * Compute the offset into the lookup tables by LOD level.
 * Return an array of offsets indexed by LOD. Also return
 * (appended to the end of the result) the lookup table size.
 * The result differs depending on whether this is the alpha
 * or brick LUT.
 * The highest LOD is at the start of each table so it by
 * definition has offset 0. Since the highest LOD level
 * just has a single brick, the second highest will start
 * at offset 1.
 */
static std::vector<std::int64_t>
calcLutOffsets(const std::vector<std::array<std::int64_t,3>>& lods_in,
               bool isalpha)
{
  std::vector<std::array<std::int64_t,3>> lods(lods_in);
  std::vector<std::int64_t> result;
  std::reverse(lods.begin(), lods.end());
  std::int64_t pos = 0;
  for (const auto& e : lods) {
    result.push_back(pos);
    pos += (e[0] * e[1] * (isalpha ? 1 : e[2]));
  }
  std::reverse(result.begin(), result.end());
  result.push_back(pos);
  return result;
}

static std::int64_t
calcEntriesInLut(const std::array<std::int64_t,3>& size,
                 const std::array<std::int64_t,3>& bricksize,
                 bool is_alpha)
{
  return calcLutOffsets(calcLodSizes(size, bricksize), is_alpha).back();
}

template<typename T>
static IHeaderAccess::podbytes_t
pod2bytes(const T& pod)
{
  IHeaderAccess::podbytes_t tmp(sizeof(pod));
  memcpy(tmp.data(), &pod, sizeof(pod));
  return tmp;
}

/**
 * Convert three arbitrary control points to 4 ordered corner points ordered
 * first il, first xl; last il, first xl, first il, last xl, last il, last xl.
 * Also return the same 4 corners in annotation- and ordinal coordinates.
 * The last two are redundant with size and annot but are provided as a
 * convenience for the api layer.
 *
 * The geometry stored in the ZGY file can be specified either as 4 ordered
 * corner points (_gdef = GridDefinition.FourPoint) or 3 arbitrary points
 * (_gdef = GridDefinition.ThreePoint) where both annotation and world
 * coordinates are given. There is also GridDefinition..Parametric which
 * specifies azimuth and spacing, but that is no longer supported. If it
 * ever was. With FourPoint the writer is still required to store the
 * annotation coordinates just in case there is some disagreement over how
 * the corners are ordered. So, we treat FourPoint the same as ThreePoint.
 * Ignore the last point (nominally this is the corner opposite the origin)
 * and treat the remaining 3 points as if they were arbitrary points
 * instead of the grid corners. V1 files are usually ThreePoint and V2 files
 * are usually FourPoint but don't count on that.
 *
 * TODO-Medium not fully implemented:
 * If the file contains world corners but no annotation, assume the writer
 * used GridDefinition.FourPoint but without storing the apparently redundant
 * annotation corners. This is contrary to spec but we will handle that case
 * for robustness.
 *
 * If the conversion fails for another reason then return all zeros because
 * in that case coordinate conversion will not be possible.
 */
static std::tuple<std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>>
calcOrderedCornersOld(const std::array<float,3>& orig,
                   const std::array<float,3>& inc,
                   const std::array<std::int64_t,3>& size,
                   const std::array<float,4>& gpiline,
                   const std::array<float,4>& gpxline,
                   const std::array<double,4>& gpx,
                   const std::array<double,4>& gpy)
{
  OrderedCornerPoints ocp(
      orig[0], inc[0], size[0], // il
      orig[1], inc[1], size[1], // xl
      gpiline[0], gpxline[0], gpx[0], gpy[0],
      gpiline[1], gpxline[1], gpx[1], gpy[1],
      gpiline[2], gpxline[2], gpx[2], gpy[2]);
  // Testing OrderedCornerPoints...
  std::array<std::array<double,2>,4> check_index{{
      { 0, 0 },
      { (double)size[0]-1, 0 },
      { 0, (double)size[1]-1 },
      { (double)size[0]-1, (double)size[1]-1}}};
  if (ocp.index_coords() != check_index)
    throw OpenZGY::Errors::ZgyInternalError("Messed up corner point computation.");
  return std::make_tuple(ocp.index_coords(),
                         ocp.annot_coords(),
                         ocp.world_coords());
}

#if USE_NEW_LATTICE
/**
 * Alternative implementation. The generalTransform() method it calls
 * is a lot simpler than OrderedCornerPoints but I worry about testing.
 * OrderedCornerPoints was cribbed from the old ZGY accessor so it
 * should be ok... but I might still be using it incorrectly.
 * TODO-Test export both calcOrderedCorners() methods and run some
 * exhaustive tests. Then switch.
 */
static std::tuple<std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>>
calcOrderedCornersNew(const std::array<float,3>& orig,
                      const std::array<float,3>& inc,
                      const std::array<std::int64_t,3>& size,
                      const std::array<float,4>& gpiline,
                      const std::array<float,4>& gpxline,
                      const std::array<double,4>& gpx,
                      const std::array<double,4>& gpy)
{
  typedef std::array<std::array<double,2>,4> coords_t;
  const double last[2]{orig[0] + inc[0] * (size[0]-1),
                       orig[1] + inc[1] * (size[1]-1)};
  double I[4] {0,    (double)size[0]-1,                 0, (double)size[0]-1};
  double J[4] {0,                    0, (double)size[1]-1, (double)size[1]-1};
  double X[4] {orig[0], last[0], orig[0], last[0]};
  double Y[4] {orig[1], orig[1], last[1], last[1]};
  coords_t index{{{I[0],J[0]}, {I[1],J[1]}, {I[2],J[2]}, {I[3],J[3]}}};
  coords_t annot{{{X[0],Y[0]}, {X[1],Y[1]}, {X[2],Y[2]}, {X[3],Y[3]}}};
  generalTransform(gpiline[0], gpxline[0],
                   gpiline[1], gpxline[1],
                   gpiline[2], gpxline[2],
                   gpx[0], gpy[0],
                   gpx[1], gpy[1],
                   gpx[2], gpy[2],
                   &X[0], &Y[0], 4);
  coords_t world{{{X[0],Y[0]}, {X[1],Y[1]}, {X[2],Y[2]}, {X[3],Y[3]}}};
  return std::make_tuple(index, annot, world);
}
#endif

static std::tuple<std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>,
                  std::array<std::array<double,2>,4>>
calcOrderedCorners(const std::array<float,3>& orig,
                   const std::array<float,3>& inc,
                   const std::array<std::int64_t,3>& size,
                   const std::array<float,4>& gpiline,
                   const std::array<float,4>& gpxline,
                   const std::array<double,4>& gpx,
                   const std::array<double,4>& gpy)
{
  auto oldresult = calcOrderedCornersOld(orig, inc, size, gpiline, gpxline, gpx, gpy);
#if USE_NEW_LATTICE && TEST_NEW_LATTICE
  auto newresult = calcOrderedCornersNew(orig, inc, size, gpiline, gpxline, gpx, gpy);
  static auto fmt = [](const std::array<std::array<double,2>,4>& c) {
                      std::stringstream ss;
                      ss << "(" << c[0][0] << ", " << c[0][1] << "), "
                         << "(" << c[1][0] << ", " << c[1][1] << "), "
                         << "(" << c[2][0] << ", " << c[2][1] << "), "
                         << "(" << c[3][0] << ", " << c[3][1] << ")";
                      return ss.str();
                    };
  std::array<std::array<double,2>,4> old_index, old_annot, old_world;
  std::array<std::array<double,2>,4> new_index, new_annot, new_world;
  std::tie(old_index, old_annot, old_world) = oldresult;
  std::tie(new_index, new_annot, new_world) = newresult;
  std::cout << "calcOrderedCorners() debugging\n"
            << "   Lattice as stored on file\n"
            << "      il beg = " << orig[0]
            << " step = " << inc[0]
            << " size = " << size[0] << "\n"
            << "      xl beg = " << orig[1]
            << " step = " << inc[1]
            << " size = " << size[1] << "\n"
            << "      point 0 annot = (" << gpiline[0] << ", " <<  gpxline[0]
            << "), world = (" << gpx[0] << ", " << gpy[0] << ")\n"
            << "      point 1 annot = (" << gpiline[1] << ", " <<  gpxline[1]
            << "), world = (" << gpx[1] << ", " << gpy[1] << ")\n"
            << "      point 2 annot = (" << gpiline[2] << ", " <<  gpxline[2]
            << "), world = (" << gpx[2] << ", " << gpy[2] << ")\n"
            << "   Old converter:\n"
            << "      index " << fmt(old_index) << "\n"
            << "      annot " << fmt(old_annot) << "\n"
            << "      world " << fmt(old_world) << "\n"
            << "   New converer:\n"
            << "      index " << fmt(new_index) << "\n"
            << "      annot " << fmt(new_annot) << "\n"
            << "      world " << fmt(new_world) << "\n"
            << std::flush;

  // The test might not be appropriate when testing corner cases.
  // E.g. for missing annotation but existing world coords
  // the new code might be doing a better job.
  for (int ii=0; ii<4; ++ii)
    for (int jj=0; jj<2; ++jj)
      if (std::abs(old_index[ii][jj] - new_index[ii][jj]) > 0.0)
        throw new std::runtime_error("calcOrderedCorners mismatch in index");
  for (int ii=0; ii<4; ++ii)
    for (int jj=0; jj<2; ++jj)
      if (std::abs(old_annot[ii][jj] - new_annot[ii][jj]) > 0.001)
        throw new std::runtime_error("calcOrderedCorners mismatch in annot");
  for (int ii=0; ii<4; ++ii)
    for (int jj=0; jj<2; ++jj)
      if (std::abs(old_world[ii][jj] - new_world[ii][jj]) > 0.1)
        throw new std::runtime_error("calcOrderedCorners mismatch in index");
#endif
  return oldresult;
}

/**
 * Sanity check. If the codingrange for an int cube is bad, silently
 * use a range that causes no conversion between storage and float.
 *
 * This avoids several corner cases both inside OpenZGY and in applications.
 *
 * Rationale: A non-finite range is always bad. A range with min==max
 * is technically valid when reading, as all storage values would map
 * to the same float value. But this is almost certainly not what the
 * writer intended. Similarly a range with min>max technically means
 * that increasing storage values correspond to decreasing float values.
 * Again, the values are more likely to be completely bogus.
 *
 * Leave the range alone if this is a float cube. There is no conversion
 * involved, and for float the codingrange is ignored by the API anyway.
 *
 * For files written by the OpenZGY library an exception wouls be thrown
 * on create. So the codingrange should always be ok for those.
 *
 * Note: The sanity check could also be applied to the histogram range.
 * That fix would also apply to float cubes. The histogram is less
 * important though, and it should be ok to let the application worry
 * about a bad histogram.
 */
static void
_fix_codingrange(float *lo, float *hi, std::uint8_t datatype)
{
  // The datatype codes used in the file map 1:1 to the RawDataType enum.
  const RawDataTypeDetails details(static_cast<RawDataType>(datatype));
  if (details.is_integer) {
    if (!std::isfinite(*lo) || !std::isfinite(*hi) || *hi <= *lo) {
      if (false)
        std::cout << "Bad codingrange "
                  << "(" << *lo << ", " << *hi << ") "
                  << "-> use "
                  << "(" << details.lowest << ", " << details.highest << ")"
                  << std::endl;
      *lo = static_cast<float>(details.lowest);
      *hi = static_cast<float>(details.highest);
    }
  }
}

} // namespace

IHeaderAccess::~IHeaderAccess()
{
}

/////////////////////////////////////////////////////////////////////////////
//    FileHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Python layout is: <4s I
// Size in storage: 8 bytes
// Thread safety: None. The user of this type is responsible for that.
#pragma pack(1)
class FileHeaderPOD
{
public:
  std::uint8_t      _magic[4];      // Always VBS\0.
  std::uint32_t     _version;       // Current version is 3, or 4 if ZFP used.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class FileHeaderAccess : public IFileHeaderAccess
{
public:
  FileHeaderPOD _pod;
  FileHeaderAccess() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset) override;
  void byteswap() override;
  void dump(std::ostream& out, const std::string& prefix = "") override;
public:
  std::array<std::uint8_t,4> magic() const override { return ptr_to_array<std::uint8_t,4>(_pod._magic); }
  std::uint32_t version() const override { return _pod._version; }
  void set_version(std::uint32_t value) override { _pod._version = value; }
};

void
FileHeaderAccess::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset)
{
  file->xx_read(&this->_pod, offset, sizeof(this->_pod));
  byteswap();
}

void
FileHeaderAccess::byteswap()
{
  // byteswap not needed for std::uint8_t _magic[4] because of its type.
  byteswapT(&_pod._version);
}

void
FileHeaderAccess::dump(std::ostream& out, const std::string& prefix)
{
  out
    << prefix << "magic():            " << array_to_string<std::uint8_t,4>(magic()) << "\n"
    << prefix << "version():          " << version() << "\n"
  ;
}

std::shared_ptr<IFileHeaderAccess>
HeaderAccessFactory::createFileHeader()
{
  return std::shared_ptr<IFileHeaderAccess>(new FileHeaderAccess());
}

/////////////////////////////////////////////////////////////////////////////
//    OffsetHeader   ////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Thread safety: None. The user of this type is responsible for that.
class OffsetHeaderAccess : public IOffsetHeaderAccess
{
  void dump(std::ostream& out, const std::string& prefix = "") override;
};

void
OffsetHeaderAccess::dump(std::ostream& out, const std::string& prefix)
{
  out
    << std::hex
    << prefix << "infoff():           " << infoff() << "\n"
    << prefix << "stroff():           " << "N/A\n"
    << prefix << "alphalupoff():      " << alphalupoff() << "\n"
    << prefix << "bricklupoff():      " << bricklupoff() << "\n"
    << prefix << "histoff():          " << histoff() << "\n"
    << std::dec
  ;
}

// Python layout is: <8I
// Size in storage: 32 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Offset Header version 1. */
class OffsetHeaderV1POD
{
public:
  std::int64_t      _infoff;        // InfoHeader position in file.
  //std::int64_t      _stroff;        // String table position, N/A in V1 and pulled from InfoHeader in V2.
  std::int64_t      _alphalupoff;   // Alpha tile lookup table position in file.
  std::int64_t      _bricklupoff;   // Brick data lookup table position in file.
  std::int64_t      _histoff;       // Histogram position in file.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class OffsetHeaderV1Access : public OffsetHeaderAccess
{
public:
  OffsetHeaderV1POD _pod;
  OffsetHeaderV1Access() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  // Set by calculate().
  std::int64_t _derived_infsize;
  std::int64_t _derived_histsize;
  std::int64_t _derived_alphalupsize;
  std::int64_t _derived_bricklupsize;

  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset) override;
  void byteswap() override;
  void calculate(const std::shared_ptr<IInfoHeaderAccess>& ih) override;
public:
  std::int64_t      infoff()      const override { return align(_pod._infoff); }
  std::int64_t      stroff()      const override { return 0; }
  std::int64_t      alphalupoff() const override { return align(_pod._alphalupoff); }
  std::int64_t      bricklupoff() const override { return align(_pod._bricklupoff); }
  std::int64_t      histoff()     const override { return align(_pod._histoff); }
  std::int64_t      infsize()      const override { return _derived_infsize; }
  std::int64_t      histsize()     const override { return _derived_histsize; }
  std::int64_t      alphalupsize() const override { return _derived_alphalupsize; }
  std::int64_t      bricklupsize() const override { return _derived_bricklupsize; }
};

void
OffsetHeaderV1Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset)
{
  // In V1, V2, V3 all the headers are fixed size, so we might as well check.
  if (sizeof(this->_pod) != 32)
    throw OpenZGY::Errors::ZgyFormatError("Wrong OffsetHeaderV1 size, expected 32 bytes");
  file->xx_read(&this->_pod, offset, sizeof(this->_pod));
  byteswap();
}

void
OffsetHeaderV1Access::byteswap()
{
  byteswapV1Long(&_pod._infoff);
  // byteswap not needed for std::int64_t _stroff because it is not stored.
  byteswapV1Long(&_pod._alphalupoff);
  byteswapV1Long(&_pod._bricklupoff);
  byteswapV1Long(&_pod._histoff);
}

// Python layout is: <B
// Size in storage: 1 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Offser Header version 2 and 3 and 4 (empty). */
class OffsetHeaderV2POD
{
public:
  //std::int64_t      _infoff;        // InfoHeader position in file.
  //std::int64_t      _stroff;        // String table position, N/A in V1 and pulled from InfoHeader in V2.
  //std::int64_t      _alphalupoff;   // Alpha tile lookup table position in file.
  //std::int64_t      _bricklupoff;   // Brick data lookup table position in file.
  //std::int64_t      _histoff;       // Histogram position in file.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class OffsetHeaderV2Access : public OffsetHeaderAccess
{
public:
  // Set by calculate().
  std::int64_t _derived_infoff;        // InfoHeader position in file.
  std::int64_t _derived_stroff;        // String table position, N/A in V1 and pulled from InfoHeader in V2.
  std::int64_t _derived_alphalupoff;   // Alpha tile lookup table position in file.
  std::int64_t _derived_bricklupoff;   // Brick data lookup table position in file.
  std::int64_t _derived_histoff;       // Histogram position in file.
  std::int64_t _derived_infsize;
  std::int64_t _derived_strsize;
  std::int64_t _derived_histsize;
  std::int64_t _derived_alphalupsize;
  std::int64_t _derived_bricklupsize;
public:
  // There is no _pod, and read() and byteswap() do nothing
  // Because V2 and V3 have no offset header on the file.
  // The data members declared here are computed from other headers.
  // For historical reasons the header still occupies one byte on file.
  podbytes_t podbytes() const override { return podbytes_t(1); }
  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset) override;
  void byteswap() override;
  void calculate(const std::shared_ptr<IInfoHeaderAccess>& ih) override;
  std::int64_t      infoff()      const override { return sizeof(FileHeaderPOD)+1; }
  std::int64_t      stroff()      const override { return align(_derived_stroff); }
  std::int64_t      alphalupoff() const override { return align(_derived_alphalupoff); }
  std::int64_t      bricklupoff() const override { return align(_derived_bricklupoff); }
  std::int64_t      histoff()     const override { return align(_derived_histoff); }
  std::int64_t      infsize()      const override { return align(_derived_infsize); }
  std::int64_t      histsize()     const override { return align(_derived_histsize); }
  std::int64_t      alphalupsize() const override { return align(_derived_alphalupsize); }
  std::int64_t      bricklupsize() const override { return align(_derived_bricklupsize); }
};

void
OffsetHeaderV2Access::read(const std::shared_ptr<IFileADT>&, std::int64_t)
{
  // Not stored on file in V1 and V2.
}

void
OffsetHeaderV2Access::byteswap()
{
  // Not stored on file in V1 and V2.
}

std::shared_ptr<IOffsetHeaderAccess>
HeaderAccessFactory::createOffsetHeader(std::uint32_t version)
{
  switch (version) {
  case 1: return std::shared_ptr<IOffsetHeaderAccess>(new OffsetHeaderV1Access());
  case 2:
  case 3:
  case 4: return std::shared_ptr<IOffsetHeaderAccess>(new OffsetHeaderV2Access());
  default: throw OpenZGY::Errors::ZgyFormatError("Unsupported ZGY version");
  }
}

/////////////////////////////////////////////////////////////////////////////
//    InfoHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Thread safety: None. The user of this type is responsible for that.
class InfoHeaderAccess : public IInfoHeaderAccess
{
public:
  void dump(std::ostream& out, const std::string& prefix = "") override;
  // Derived information computed from other virtual members, which means
  // they can be implemented in the base class. Note that all of these
  // could be cached if desired. They cannot change after a file is created.
  std::int64_t bytesperalpha() const override;
  std::int64_t bytesperbrick() const override;
  std::int64_t bytespersample() const override;
  std::array<double,2> storagetofloat() const override;
  double storagetofloat_slope() const override;
  double storagetofloat_intercept() const override;
  double defaultstorage() const override;
  double defaultvalue() const override;
};

void
InfoHeaderAccess::dump(std::ostream& out, const std::string& prefix)
{
  out
    << prefix << "bricksize():        " << array_to_string(bricksize()) << "\n"
    << prefix << "datatype():         " << (int)datatype() << "\n"
    << prefix << "safe_codingrange(): " << array_to_string(safe_codingrange()) << "\n"
    << prefix << "dataid():           " << array_to_hex(dataid()) << "\n"
    << prefix << "verid():            " << array_to_hex(verid()) << "\n"
    << prefix << "previd():           " << array_to_hex(previd()) << "\n"
    << prefix << "srcname():          " << srcname() << "\n"
    << prefix << "srcdesc():          " << srcdesc() << "\n"
    //<< prefix << "srctype():          " << srctype() << "\n"
    << prefix << "orig():             " << array_to_string(orig()) << "\n"
    << prefix << "inc():              " << array_to_string(inc()) << "\n"
    << prefix << "size():             " << array_to_string(size()) << "\n"
    //<< prefix << "curorig():          " << array_to_string(curorig()) << "\n"
    //<< prefix << "cursize():          " << array_to_string(cursize()) << "\n"
    << prefix << "scnt():             " << scnt() << "\n"
    << prefix << "ssum():             " << ssum() << "\n"
    << prefix << "sssq():             " << sssq() << "\n"
    << prefix << "smin():             " << smin() << "\n"
    << prefix << "smax():             " << smax() << "\n"
    //<< prefix << "srvorig():          " << array_to_string(srvorig()) << "\n"
    //<< prefix << "srvsize():          " << array_to_string(srvsize()) << "\n"
    //<< prefix << "gdef():             " << gdef() << "\n"
    //<< prefix << "gazim():            " << array_to_string(gazim()) << "\n"
    //<< prefix << "gbinsz():           " << array_to_string(gbinsz()) << "\n"
    << prefix << "gpiline():          " << array_to_string(gpiline()) << "\n"
    << prefix << "gpxline():          " << array_to_string(gpxline()) << "\n"
    << prefix << "gpx():              " << array_to_string(gpx()) << "\n"
    << prefix << "gpy():              " << array_to_string(gpy()) << "\n"
    << prefix << "hprjsys():          " << hprjsys() << "\n"
    << prefix << "hdim():             " << (int)hdim() << "\n"
    << prefix << "hunitfactor():      " << hunitfactor() << "\n"
    << prefix << "hunitname():        " << hunitname() << "\n"
    << prefix << "vdim():             " << (int)vdim() << "\n"
    << prefix << "vunitfactor():      " << vunitfactor() << "\n"
    << prefix << "vunitname():        " << vunitname() << "\n"
    << prefix << "slbufsize():        " << slbufsize() << "\n"
    << prefix << "nlods():            " << nlods() << "\n"
    << prefix << "bytesperalpha():    " << bytesperalpha() << "\n"
    << prefix << "bytesperbrick():    " << bytesperbrick() << "\n"
    << prefix << "bytespersample():   " << bytespersample() << "\n"
    << prefix << "storagetofloat():   " << array_to_string(storagetofloat()) << "\n"
    << prefix << "defaultstorage():   " << defaultstorage() << "\n"
    << prefix << "defaultvalue():     " << defaultvalue() << "\n"
    ;
  //out << prefix << "ocp_world():        " << "(derived)\n";
  out << prefix << "lodsizes():       ";
  for (const auto& it : lodsizes())
    out << "  " << array_to_string(it);
  out << "\n";
  out << prefix << "alphaoffsets():   " << std::hex;
  for (const auto& it : alphaoffsets())
    out << "  " << std::setw(4) << it;
  out << std::dec << "\n";
  out << prefix << "brickoffsets():   " << std::hex;
  for (const auto& it : brickoffsets())
    out << "  " << std::setw(4) << it;
  out << std::dec << "\n";
}

std::int64_t
InfoHeaderAccess::bytesperalpha() const
{
  std::array<std::int64_t,3> bs = this->bricksize();
  return bs[0] * bs[1] * sizeof(std::uint8_t);
}

std::int64_t
InfoHeaderAccess::bytesperbrick() const
{
  std::array<std::int64_t,3> bs = this->bricksize();
  return bs[0] * bs[1] * bs[2] * this->bytespersample();
}

std::int64_t
InfoHeaderAccess::bytespersample() const
{
  return (std::int64_t)RawDataTypeDetails(this->datatype()).size;
}

/**
 * Get the linear transform y = a*x + b for converting from
 * storage values to actual floating point values.
 * The smallest value in storage (e.g. -128) should map to codingrange[0]
 * and the largest value (e.g. +127) should map to codingrange[1].
 * If file_dtype is a floating point type there is never any conversion.
 *
 * Note that the codingrange stored in file is two float32 numbers i.e.
 * 24 bits mantissa. This is ok for our purpose because it is only used
 * for converting int8 and int16 data so the user doesn't expect anything
 * more than 16 bits of precision. But to avoid numerical issues with
 * the transform the codingrange should immediately be upcast to float64
 * befor being used. Case in point: (double)(hi-lo) is less accurate than
 * ((double)hi - (double)lo) when hi,lo are float32. Admittedly only by
 * 1 ULPS but the error might be amplified later. Enough to raise errors
 * when testing even if the signal still has 16 bits of precision intact.
 *
 * See also ZgyInternalBulk._scaleToFloat in OpenZGY/Python.
 */
std::array<double,2>
InfoHeaderAccess::storagetofloat() const
{
  double slope = 1.0, intercept = 0.0;
  const RawDataTypeDetails info(this->datatype());
  if (info.is_integer) {
    const std::array<float,2> range = this->safe_codingrange();
    slope = ((double)range[1]-(double)range[0]) / (info.highest - info.lowest);
    intercept = range[0] - slope * info.lowest;
#if 0 // Debug issues with numerical inaccuracy.
    static double debug_old_slope = 0, debug_old_intercept = 0;
    if (debug_old_slope != slope || debug_old_intercept != intercept) {
      std::cerr << "@@ storagetofloat " << std::setprecision(14)
                << " slope " << slope
                << " intercept " << intercept
                << " range " << range[0] << " " << range[1]
                << " test " << slope*-32768+intercept
                << " " << slope*32767+intercept
                << std::endl;
      debug_old_slope = slope;
      debug_old_intercept = intercept;
    }
#endif
  }
  return std::array<double,2>{slope,intercept};
}

double
InfoHeaderAccess::storagetofloat_slope() const
{
  return storagetofloat()[0];
}

double
InfoHeaderAccess::storagetofloat_intercept() const
{
  return storagetofloat()[1];
}

/**
 * Return the storage value to be used for data that has never been written.
 * This is supposed to be the storage value that, when converted to float,
 * is as close to zero as possible. I.e. it is simply zero converted
 * from float to storage values.
 *
 * Implementation note: Should have called RoundAndClip<T>() here,
 * but that method is templated which makes it awkward to use here.
 * Try duplicating the logic in RoundAndClip as close as possible.
 */
double InfoHeaderAccess::defaultstorage() const
{
  const RawDataTypeDetails info(this->datatype());
  double value = (0 - storagetofloat_intercept()) / storagetofloat_slope();
  if (!info.is_integer)
    return 0;
  else if (value <= info.lowest)
    return info.lowest;
  else if (value >= info.highest)
    return info.highest;
  else if (value < 0)
    return static_cast<double>(static_cast<std::int64_t>(value - 0.5));
  else
    return static_cast<double>(static_cast<std::int64_t>(value + 0.5));
}

/**
 * Return the float value to be used for data that has never been written.
 * This is supposed to be zero afyer a float/storage/float round trip.
 */
double InfoHeaderAccess::defaultvalue() const
{
  return defaultstorage() * storagetofloat_slope() + storagetofloat_intercept();
}

// Python layout is: <3i 3i 3i 3f 4i 4i 4d 4d B B
// Size in storage: 146 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Info Header version 1. */
class InfoHeaderV1POD
{
public:
  std::int32_t      _size[3];       // Integer size in inline, crossline, vertical directions.
  std::int32_t      _orig[3];       // First inline, crossline, time/depth. Only integral values allowed.
  std::int32_t      _inc[3];        // Integer increment in inline, crossline, vertical directions.
  float             _incfactor[3];  // Unused. Write as (1,1,1), ignore on read.
  std::int32_t      _gpiline[4];    // Inline component of 4 control points.
  std::int32_t      _gpxline[4];    // Crossline component of 4 control points.
  double            _gpx[4];        // X coordinate of 4 control points.
  double            _gpy[4];        // Y coordinate of 4 control points.
  std::uint8_t      _datatype;      // Type of samples in each brick: int8 = 0, int16 = 2, float32 = 6.
  std::uint8_t      _coordtype;     // Coordinate type: unknown = 0, meters = 1, feet = 2, degrees*3600 = 3, degrees = 4, DMS = 5.
  // Note, _datatype and _coordtype are V1 specific and need special handling in accessor.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class InfoHeaderV1Access : public InfoHeaderAccess
{
public:
  InfoHeaderV1POD _pod;
  float _cached_file_sample_min;
  float _cached_file_sample_max;
  float _cached_safe_codingrange_min;
  float _cached_safe_codingrange_max;
  std::int32_t _cached_nlods;
  std::vector<std::array<std::int64_t,3>> _cached_lodsizes;
  std::vector<std::int64_t> _cached_alphaoffsets;
  std::vector<std::int64_t> _cached_brickoffsets;
  std::array<std::array<double,2>,4> _cached_index;
  std::array<std::array<double,2>,4> _cached_annot;
  std::array<std::array<double,2>,4> _cached_world;
  InfoHeaderV1Access() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  void       read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) override;
  void       byteswap() override;
  void       calculate_cache() override;
  void       calculate_read(const podbytes_t& slbuf, const std::shared_ptr<IHistHeaderAccess>& hh) override;
  podbytes_t calculate_write() override;
public:
  std::array<std::int64_t,3> size() const override { return array_cast<std::int64_t,std::int32_t,3>(ptr_to_array<std::int32_t,3>(_pod._size)); }
  std::array<float,3>        orig() const override { return array_cast<float,std::int32_t,3>(ptr_to_array<std::int32_t,3>(_pod._orig)); }
  std::array<float,3>        inc()  const override { return array_cast<float,std::int32_t,3>(ptr_to_array<std::int32_t,3>(_pod._inc)); }
  //unused: virtual std::array<float,3>        incfactor() const override { return ptr_to_array<float,3>(_pod._incfactor); }
  std::array<float,4>        gpiline() const override { return array_cast<float,std::int32_t,4>(ptr_to_array<std::int32_t,4>(_pod._gpiline)); }
  std::array<float,4>        gpxline() const override { return array_cast<float,std::int32_t,4>(ptr_to_array<std::int32_t,4>(_pod._gpxline)); }
  std::array<double,4>       gpx() const override { return ptr_to_array<double,4>(_pod._gpx); }
  std::array<double,4>       gpy() const override { return ptr_to_array<double,4>(_pod._gpy); }
  RawDataType                datatype() const override { return DecodeDataType(_pod._datatype); }
  //special: virtual std::uint8_t               coordtype() const override { return _pod._coordtype; }
  std::array<std::int64_t,3> bricksize() const override { return std::array<std::int64_t,3>{64, 64, 64}; }
  std::array<float,2>        safe_codingrange() const override { return std::array<float,2>{_cached_safe_codingrange_min, _cached_safe_codingrange_max}; }
  std::array<float,2>        raw_codingrange() const override { return std::array<float,2>{_cached_file_sample_min, _cached_file_sample_max}; }
  std::array<std::uint8_t,16>dataid() const override { return std::array<std::uint8_t,16>{0}; }
  std::array<std::uint8_t,16>verid() const override { return std::array<std::uint8_t,16>{0}; }
  std::array<std::uint8_t,16>previd() const override { return std::array<std::uint8_t,16>{0}; }
  std::string                srcname() const override { return std::string(); }
  std::string                srcdesc() const override { return std::string(); }
  RawDataType                srctype() const override { return datatype(); }
  //unused: virtual std::array<std::int32_t,3> curorig() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  //unused: virtual std::array<std::int32_t,3> cursize() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  std::int64_t               scnt() const override { return 0; }
  double                     ssum() const override { return 0; }
  double                     sssq() const override { return 0; }
  float                      smin() const override { return align(_cached_file_sample_min); }
  float                      smax() const override { return align(_cached_file_sample_max); }
  //unused: std::array<float,3>        srvorig() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  //unused: std::array<float,3>        srvsize() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  //unused: std::uint8_t               gdef() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  //unused: std::array<double,2>       gazim() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  //unused: std::array<double,2>       gbinsz() const override { throw OpenZGY::Errors::ZgyInternalError("Not implemented"); }
  std::string                hprjsys() const override { return std::string(); }
  RawHorizontalDimension     hdim() const override {
    switch (DecodeCoordType(_pod._coordtype)) {
    default:
    case RawCoordType::Unknown:      return RawHorizontalDimension::Unknown;
    case RawCoordType::Meters:       return RawHorizontalDimension::Length;
    case RawCoordType::Feet:         return RawHorizontalDimension::Length;
    case RawCoordType::ArcSec:       return RawHorizontalDimension::ArcAngle;
    case RawCoordType::ArcDeg:       return RawHorizontalDimension::ArcAngle;
    //case RawCoordType::ArcDegMinSec: return RawHorizontalDimension::ArcAngle;
    // ArcDegMinSec is unsupported. Does not work in the old code either.
    }
  }
  double                     hunitfactor() const override {
    switch (DecodeCoordType(_pod._coordtype)) {
    default:
    case RawCoordType::Unknown:      return 1.0;
    case RawCoordType::Meters:       return 1.0;
    case RawCoordType::Feet:         return 0.3048;
    case RawCoordType::ArcSec:       return 3600.0;
    case RawCoordType::ArcDeg:       return 1.0;
    case RawCoordType::ArcDegMinSec: return 1.0;
    }
  }
  std::string                hunitname() const override {
    switch (DecodeCoordType(_pod._coordtype)) {
    default:
    case RawCoordType::Unknown:      return "";
    case RawCoordType::Meters:       return "m";
    case RawCoordType::Feet:         return "ft";
    case RawCoordType::ArcSec:       return "arcsec";
    case RawCoordType::ArcDeg:       return "deg";
    case RawCoordType::ArcDegMinSec: return "DMS";
    }
  }
  // V1 did not store the vertical unit at all.
  RawVerticalDimension       vdim() const override { return RawVerticalDimension::Unknown; }
  double                     vunitfactor() const override { return 1.0; }
  std::string                vunitname() const override { return std::string(); }
  std::uint32_t              slbufsize() const override { return 0; }
  const std::array<std::array<double,2>,4>& ocp_index() const override { return _cached_index; }
  const std::array<std::array<double,2>,4>& ocp_annot() const override { return _cached_annot; }
  const std::array<std::array<double,2>,4>& ocp_world() const override { return _cached_world; }
  const std::vector<std::array<std::int64_t,3>>& lodsizes() const override { return _cached_lodsizes; }
  std::int32_t nlods() const override { return _cached_nlods; }
  const std::vector<std::int64_t>& alphaoffsets() const override { return _cached_alphaoffsets; }
  const std::vector<std::int64_t>& brickoffsets() const override { return _cached_brickoffsets; }
  // Write support.
  void setstats(std::int64_t scnt, double ssum, double sssq,
                        double smin, double smax) override
  {
    throw OpenZGY::Errors::ZgyInternalError("Writing InfoHeader is only supported for the latest version.");
  }
};

void
InfoHeaderV1Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size)
{
  // In V1, V2, V3 all the headers are fixed size, so we might as well check.
  if (size != 146 || size != sizeof(this->_pod))
    throw OpenZGY::Errors::ZgyFormatError("Wrong InfoHeaderV1 size, expected 146 bytes");
  memset(&this->_pod, 0, sizeof(this->_pod));
  file->xx_read(&this->_pod, offset, std::min(sizeof(this->_pod), (size_t)size));
  byteswap();
}

void
InfoHeaderV1Access::byteswap()
{
  byteswapT(&_pod._size[0], 3);
  byteswapT(&_pod._orig[0], 3);
  byteswapT(&_pod._inc[0], 3);
  byteswapT(&_pod._incfactor[0], 3);
  byteswapT(&_pod._gpiline[0], 4);
  byteswapT(&_pod._gpxline[0], 4);
  byteswapT(&_pod._gpx[0], 4);
  byteswapT(&_pod._gpy[0], 4);
}

// Python layout is: <3i B 2f 16s 16s 16s   B 3f 3f 3i 3i 3i q d d f f 3f 3f B 2d 2d 4f 4f 4d 4d  B d  B d  I
// Size in storage: 337 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Info Header version 2 and 3 and 4. */
class InfoHeaderV2POD
{
public:
  std::int32_t      _bricksize[3];  // Brick size. Values other than (64,64,64) will likely not work.
  std::uint8_t      _datatype;      // Type of samples in each brick: int8 = 0, int16 = 2, float32 = 6.
  float             _file_codingrange[2]; // If datatype is integral, this is the value range samples will be scaled to when read as float. In this case it must be specified on file creation. If datatype is float then this is the value range of the data and should be set automatically when writing the file.
  std::uint8_t      _dataid[16];    // GUID set on file creation.
  std::uint8_t      _verid[16];     // GUID set each time the file is changed.
  std::uint8_t      _previd[16];    // GUID before last change.
  //char*             _srcname;       // Optional name of this data set. Rarely used.
  //char*             _srcdesc;       // Optional description of this data set. Rarely used.
  std::uint8_t      _srctype;       // Optional datatype the samples had before being stored in this file.
  float             _orig[3];       // First inline, crossline, time/depth. Unlike v1 these are now floating point.
  float             _inc[3];        // Increment in inline, crossline, vertical directions.
  std::int32_t      _size[3];       // Size in inline, crossline, vertical directions.
  std::int32_t      _curorig[3];    // Unused. Set to (0,0,0) on write and ignore on read.
  std::int32_t      _cursize[3];    // Unused. Set to size on write and ignore on read.
  std::int64_t      _scnt;          // Count of values used to compute statistics.
  double            _ssum;          // Sum of all "scnt" values.
  double            _sssq;          // Sum of squared "scnt" values.
  float             _smin;          // Statistical (computed) minimum value.
  float             _smax;          // Statistical (computed) maximum value.
  float             _srvorig[3];    // Unused. Set equal to orig on write. Ignore on read.
  float             _srvsize[3];    // Unused. Set to inc*size on write. Ignore on read.
  std::uint8_t      _gdef;          // Grid definition type. Set to 3 (enum: "FourPoint") on write. Ignored on read. See notes for a longer explanation.
  double            _gazim[2];      // Unused.
  double            _gbinsz[2];     // Unused.
  float             _gpiline[4];    // Inline component of 4 control points.
  float             _gpxline[4];    // Crossline component of 4 control points.
  double            _gpx[4];        // X coordinate of 4 control points.
  double            _gpy[4];        // Y coordinate of 4 control points.
  //char*             _hprjsys;       // Free form description of the projection coordinate system. Usually not parseable into a well known CRS.
  std::uint8_t      _hdim;          // Horizontal dimension. Unknown = 0, Length = 1, ArcAngle = 2. Few applications support ArcAngle.
  double            _hunitfactor;   // Multiply by this factor to convert from storage units to SI units. Applies to gpx, gpy.
  //char*             _hunitname;     // For annotation only. Use hunitfactor, not the name, to convert to or from SI.
  std::uint8_t      _vdim;          // Vertical dimension. Unknown = 0, Depth = 1, SeismicTWT = 1, SeismicOWT = 3.
  double            _vunitfactor;   // Multiply by this factor to convert from storage units to SI units. Applies to orig[2], inc[2].
  //char*             _vunitname;     // For annotation only. Use vunitfactor, not the name, to convert to or from SI.
  std::uint32_t     _slbufsize;     // Size of the StringList section.
  //enum              _ocp_world;     // Ordered corner points: ((i0,j0),(iN,j0),(i0,jM),(iN,jM)
  //std::int32_t      _lodsizes[lod]; // Size of the survey at reduced level of detail.
  //std::int32_t      _nlods;         // How many levels of details. 1 means only full resolution. Currently nlods will always be just enough to make the highest LOD (i.e. lowest resolution) fit in a single brick.
  //std::int64_t      _brickoffsets[lod]; // How many entries in the lookup table to skip when dealing with level N.
  //std::int64_t      _alphaoffsets[lod]; // How many entries in the lookup table to skip when dealing with level N.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class InfoHeaderV2Access : public InfoHeaderAccess
{
public:
  InfoHeaderV2POD _pod;
  std::string _srcname;
  std::string _srcdesc;
  std::string _hprjsys;
  std::string _hunitname;
  std::string _vunitname;
  std::int32_t _cached_nlods;
  std::vector<std::array<std::int64_t,3>> _cached_lodsizes;
  std::vector<std::int64_t> _cached_alphaoffsets;
  std::vector<std::int64_t> _cached_brickoffsets;
  std::array<std::array<double,2>,4> _cached_index;
  std::array<std::array<double,2>,4> _cached_annot;
  std::array<std::array<double,2>,4> _cached_world;
  float _cached_safe_codingrange[2];

  InfoHeaderV2Access() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  void       read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) override;
  void       byteswap() override;
  void       calculate_read(const podbytes_t& slbuf, const std::shared_ptr<IHistHeaderAccess>& hh) override;
  void       calculate_cache() override;
  podbytes_t calculate_write() override;
public:
  std::array<std::int64_t,3> bricksize() const override { return array_cast<std::int64_t,std::int32_t,3>(ptr_to_array<std::int32_t,3>(_pod._bricksize)); }
  RawDataType                datatype() const override { return DecodeDataType(_pod._datatype); }
  std::array<float,2>        safe_codingrange() const override { return ptr_to_array<float,2>(_cached_safe_codingrange); }
  std::array<float,2>        raw_codingrange() const override { return ptr_to_array<float,2>(_pod._file_codingrange); }
  std::array<std::uint8_t,16>dataid() const override { return ptr_to_array<std::uint8_t,16>(_pod._dataid); }
  std::array<std::uint8_t,16>verid() const override { return ptr_to_array<std::uint8_t,16>(_pod._verid); }
  std::array<std::uint8_t,16>previd() const override { return ptr_to_array<std::uint8_t,16>(_pod._previd); }
  std::string                srcname() const override { return _srcname; }
  std::string                srcdesc() const override { return _srcdesc; }
  RawDataType                srctype() const override { return DecodeDataType(_pod._srctype); }
  std::array<float,3>        orig() const override { return ptr_to_array<float,3>(_pod._orig); }
  std::array<float,3>        inc()  const override { return ptr_to_array<float,3>(_pod._inc); }
  std::array<std::int64_t,3> size() const override { return array_cast<std::int64_t,std::int32_t,3>(ptr_to_array<std::int32_t,3>(_pod._size)); }
  //unused: std::array<std::int32_t,3> curorig() const override { return ptr_to_array<std::int32_t,3>(_pod._curorig); }
  //unused: std::array<std::int32_t,3> cursize() const override { return ptr_to_array<std::int32_t,3>(_pod._cursize); }
  std::int64_t               scnt() const override { return align(_pod._scnt); }
  double                     ssum() const override { return align(_pod._ssum); }
  double                     sssq() const override { return align(_pod._sssq); }
  float                      smin() const override { return align(_pod._smin); }
  float                      smax() const override { return align(_pod._smax); }
  //unused: std::array<float,3> srvorig() const override { return ptr_to_array<float,3>(_pod._srvorig); }
  //unused: std::array<float,3> srvsize() const override { return ptr_to_array<float,3>(_pod._srvsize); }
  //unused: std::uint8_t gdef() const override { return _pod._gdef; }
  //unused: std::array<double,2> gazim() const override { return ptr_to_array<double,2>(_pod._gazim); }
  //unused: std::array<double,2> gbinsz() const override { return ptr_to_array<double,2>(_pod._gbinsz); }
  std::array<float,4>        gpiline() const override { return ptr_to_array<float,4>(_pod._gpiline); }
  std::array<float,4>        gpxline() const override { return ptr_to_array<float,4>(_pod._gpxline); }
  std::array<double,4>       gpx() const override { return ptr_to_array<double,4>(_pod._gpx); }
  std::array<double,4>       gpy() const override { return ptr_to_array<double,4>(_pod._gpy); }
  std::string                hprjsys() const override { return _hprjsys; }
  RawHorizontalDimension     hdim() const override { return DecodeHorizontalDimension(_pod._hdim); }
  double                     hunitfactor() const override { return align(_pod._hunitfactor); }
  std::string                hunitname() const override { return _hunitname; }
  RawVerticalDimension       vdim() const override { return DecodeVerticalDimension(_pod._vdim); }
  double                     vunitfactor() const override { return align(_pod._vunitfactor); }
  std::string                vunitname() const override { return _vunitname; }
  std::uint32_t              slbufsize() const override { return align(_pod._slbufsize); }
  const std::array<std::array<double,2>,4>& ocp_index() const override { return _cached_index; }
  const std::array<std::array<double,2>,4>& ocp_annot() const override { return _cached_annot; }
  const std::array<std::array<double,2>,4>& ocp_world() const override { return _cached_world; }
  std::int32_t nlods() const override { return _cached_nlods; }
  const std::vector<std::array<std::int64_t,3>>& lodsizes() const override { return _cached_lodsizes; }
  const std::vector<std::int64_t>& alphaoffsets() const override { return _cached_alphaoffsets; }
  const std::vector<std::int64_t>& brickoffsets() const override { return _cached_brickoffsets; }
  void setstats(std::int64_t scnt, double ssum, double sssq,
                        double smin, double smax) override
  {
    _pod._scnt = scnt;
    _pod._ssum = ssum;
    _pod._sssq = sssq;
    // TODO-Low: After the next ZGY version the casts must be removed.
    _pod._smin = static_cast<float>(smin);
    _pod._smax = static_cast<float>(smax);
  }
};

void
InfoHeaderV2Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size)
{
  // In V1, V2, V3 all the headers are fixed size, so we might as well check.
  if (size != 337 || size != sizeof(this->_pod))
    throw OpenZGY::Errors::ZgyFormatError("Wrong InfoHeaderV2 size, expected 337 bytes");
  memset(&this->_pod, 0, sizeof(this->_pod));
  file->xx_read(&this->_pod, offset, std::min(sizeof(this->_pod), (size_t)size));
  byteswap();
}

void
InfoHeaderV2Access::byteswap()
{
  byteswapT(&_pod._bricksize[0], 3);
  // byteswap not needed for std::uint8_t _datatype because of its type.
  byteswapT(&_pod._file_codingrange[0], 2);
  // byteswap not needed for std::uint8_t _dataid[16] because of its type.
  // byteswap not needed for std::uint8_t _verid[16] because of its type.
  // byteswap not needed for std::uint8_t _previd[16] because of its type.
  // byteswap not needed for char* _srcname because it is not stored.
  // byteswap not needed for char* _srcdesc because it is not stored.
  // byteswap not needed for std::uint8_t _srctype because of its type.
  byteswapT(&_pod._orig[0], 3);
  byteswapT(&_pod._inc[0], 3);
  byteswapT(&_pod._size[0], 3);
  byteswapT(&_pod._curorig[0], 3);
  byteswapT(&_pod._cursize[0], 3);
  byteswapT(&_pod._scnt);
  byteswapT(&_pod._ssum);
  byteswapT(&_pod._sssq);
  byteswapT(&_pod._smin);
  byteswapT(&_pod._smax);
  byteswapT(&_pod._srvorig[0], 3);
  byteswapT(&_pod._srvsize[0], 3);
  // byteswap not needed for std::uint8_t _gdef because of its type.
  byteswapT(&_pod._gazim[0], 2);
  byteswapT(&_pod._gbinsz[0], 2);
  byteswapT(&_pod._gpiline[0], 4);
  byteswapT(&_pod._gpxline[0], 4);
  byteswapT(&_pod._gpx[0], 4);
  byteswapT(&_pod._gpy[0], 4);
  // byteswap not needed for char* _hprjsys because it is not stored.
  // byteswap not needed for std::uint8_t _hdim because of its type.
  byteswapT(&_pod._hunitfactor);
  // byteswap not needed for char* _hunitname because it is not stored.
  // byteswap not needed for std::uint8_t _vdim because of its type.
  byteswapT(&_pod._vunitfactor);
  // byteswap not needed for char* _vunitname because it is not stored.
  byteswapT(&_pod._slbufsize);
  // byteswap not needed for enum _ocp_world because it is not stored.
  // byteswap not needed for std::int32_t _lodsizes[lod] because it is not stored.
  // byteswap not needed for std::int32_t _nlods because it is not stored.
  // byteswap not needed for std::int64_t _brickoffsets[lod] because it is not stored.
  // byteswap not needed for std::int64_t _alphaoffsets[lod] because it is not stored.
}

std::shared_ptr<IInfoHeaderAccess>
HeaderAccessFactory::createInfoHeader(std::uint32_t version)
{
  switch (version) {
  case 1: return std::shared_ptr<IInfoHeaderAccess>(new InfoHeaderV1Access());
  case 2:
  case 3:
  case 4: return std::shared_ptr<IInfoHeaderAccess>(new InfoHeaderV2Access());
  default: throw OpenZGY::Errors::ZgyFormatError("Unsupported ZGY version");
  }
}

/////////////////////////////////////////////////////////////////////////////
//    HistHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Thread safety: None. The user of this type is responsible for that.
class HistHeaderAccess : public IHistHeaderAccess
{
  void dump(std::ostream& out, const std::string& prefix) override;
};

void
HistHeaderAccess::dump(std::ostream& out, const std::string& prefix)
{
  out
    << prefix << "samplecount():      " << samplecount() << "\n"
    << prefix << "minvalue():         " << minvalue() << "\n"
    << prefix << "maxvalue():         " << maxvalue() << "\n"
    << prefix << "bins():             (array)\n"
  ;
}

// Python layout is: <f f 256I
// Size in storage: 1032 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Histogram Header version 1. */
class HistHeaderV1POD
{
public:
  float             _max;           // Center point of first bin.
  float             _min;           // Center point of last bin.
  std::uint32_t     _bin[256];      // Histogram.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class HistHeaderV1Access : public HistHeaderAccess
{
public:
  HistHeaderV1POD _pod;
private:
  std::vector<std::int64_t> _converted_bins;
public:
  HistHeaderV1Access() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) override;
  void byteswap() override;
  void calculate() override
  {
    _converted_bins.clear();
    for (int ii=0; ii<bincount(); ++ii)
      _converted_bins.push_back(_pod._bin[ii]);
  }
  void sethisto(double minvalue, double maxvalue,
                        const std::int64_t* bins, std::int64_t bincount) override
  {
    throw OpenZGY::Errors::ZgyInternalError("Writing HistHeader is only supported for the latest version.");
  }

public:
  std::int64_t bincount()    const override { return 256; }
  std::int64_t samplecount() const override { return 0; }
  double       minvalue()    const override { return align(_pod._min); }
  double       maxvalue()    const override { return align(_pod._max); }
  // The returned pointer is invalidated by a call to calculate().
  // This should not be an issue since V1 only supports reading,
  // so there isn't really any reason to call calculate() more
  // than once.
  const std::int64_t* bins() const override { return _converted_bins.data(); }
};

void
HistHeaderV1Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size)
{
  // In V1, V2, V3 all the headers are fixed size, so we might as well check.
  if (size != 1032 || size != sizeof(this->_pod))
    throw OpenZGY::Errors::ZgyFormatError("Wrong HistHeaderV1 size, expected 1032 bytes");
  memset(&this->_pod, 0, sizeof(this->_pod));
  file->xx_read(&this->_pod, offset, sizeof(this->_pod));
  byteswap();
}

void
HistHeaderV1Access::byteswap()
{
  byteswapT(&_pod._max);
  byteswapT(&_pod._min);
  byteswapT(&_pod._bin[0], 256);
}

// Python layout is: <q f f 256q
// Size in storage: 2064 bytes
#pragma pack(1)
// Thread safety: None. The user of this type is responsible for that.
/** \brief Physical layout of Histogram Header version 2 and 3 and 4. */
class HistHeaderV2POD
{
public:
  std::int64_t      _cnt;           // Total number of samples.
  float             _min;           // Center point of first bin.
  float             _max;           // Center point of last bin.
  std::int64_t      _bin[256];      // Histogram.
};
#pragma pack()

// Thread safety: None. The user of this type is responsible for that.
class HistHeaderV2Access : public HistHeaderAccess
{
public:
  HistHeaderV2POD _pod;
  HistHeaderV2Access() { memset(&_pod, 0, sizeof(_pod)); }
  podbytes_t podbytes() const override { return pod2bytes(_pod); }
  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) override;
  void byteswap() override;
  void calculate() override { }

public:
  std::int64_t bincount()    const override { return 256; }
  std::int64_t samplecount() const override { return align(_pod._cnt); }
  double       minvalue()    const override { return align(_pod._min); }
  double       maxvalue()    const override { return align(_pod._max); }
  // Should be safely aligned, checked by hand.
  const std::int64_t* bins() const override { return _pod._bin; }
  // Write support.
  void sethisto(double minvalue, double maxvalue,
                        const std::int64_t* bins, std::int64_t bincount) override;
};

void
HistHeaderV2Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size)
{
  // In V1, V2, V3 all the headers are fixed size, so we might as well check.
  if (size != 2064 || size != sizeof(this->_pod))
    throw OpenZGY::Errors::ZgyFormatError("Wrong HistHeaderV2 size, expected 2064 bytes");
  memset(&this->_pod, 0, sizeof(this->_pod));
  file->xx_read(&this->_pod, offset, std::min(sizeof(this->_pod), (size_t)size));
  byteswap();
}

void
HistHeaderV2Access::byteswap()
{
  byteswapT(&_pod._cnt);
  byteswapT(&_pod._min);
  byteswapT(&_pod._max);
  byteswapT(&_pod._bin[0], 256);
}

void HistHeaderV2Access::sethisto(
    double minvalue, double maxvalue,
    const std::int64_t* bins, std::int64_t bincount)
{
  const std::size_t nbins = sizeof(_pod._bin) / sizeof(_pod._bin[0]);
  if (bins == nullptr || bincount == 0) {
    // Missing bins means clear the histogram.
    _pod._cnt = 0;
    for (std::size_t ii = 0; ii < nbins; ++ii) {
      _pod._bin[ii] = 0;
    }
  }
  else {
    // If the in-memory representation has a different size then the caller
    // is expected to do the expand or shrink. Not us.
    // TODO-Low: Robustness: Unfortunately the caller needs to hard code 256
    // as the size of the histogram stored on file. Had this code been able
    // to resize the histogram then that problem goes away.
    if (bincount != (std::int64_t)nbins)
      throw OpenZGY::Errors::ZgyInternalError("Histogram resize not implemented");
    _pod._cnt = 0;
    for (std::size_t ii = 0; ii < nbins; ++ii) {
      _pod._bin[ii] = bins[ii];
      _pod._cnt += bins[ii];
    }
  }
  // TODO-Low: Some kind of static_assert making the following cast safe.
  // Note reason for cast: ZGY only deals with int8, int16, float
  // which means that converting the value range to float will not
  // lose any precision. If we later decide to support int32 then
  // this and many similar places need to be updated.
  _pod._min = static_cast<float>(minvalue);
  _pod._max = static_cast<float>(maxvalue);
}

std::shared_ptr<IHistHeaderAccess>
HeaderAccessFactory::createHistHeader(std::uint32_t version)
{
  switch (version) {
  case 1: return std::shared_ptr<IHistHeaderAccess>(new HistHeaderV1Access());
  case 2:
  case 3:
  case 4: return std::shared_ptr<IHistHeaderAccess>(new HistHeaderV2Access());
  default: throw OpenZGY::Errors::ZgyFormatError("Unsupported ZGY version");
  }
}

/////////////////////////////////////////////////////////////////////////////
//    AlphaLUP, BrickLUP   //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Both the Alpha lookup table and the Brick lookup table hold a 64-bit
 * file offset for each tile or brick in the file. Alpha tiles are bitmaps
 * used to flag dead traces, and only have (i,j) coordinates. Bricks
 * contain the actual samples and are indexed with (i, j, k).
 *
 * The size of the lookup tables depend on the survey size. The first
 * entry in the lookup table is for the brick or tile (always just one)
 * holding the lowest resolution. This is immediately followed by one
 * or more entries for the bricks or tiles at level of detail N-1, and
 * so on until the entries for lod 0. Within one lod level the first
 * entry is for the lowest numbered i,j,[k]. For subsequent entries the
 * i numbers vary fastest and the k (or j in the alpha case) varies
 * slowest. Note that this is somewhat non intuitive as it is the
 * opposite of the ordering of samples within a tile.
 *
 * In version 1 of the lookup tables the file offsets are stored in a
 * somewhat quirky manner. The high 32 bits and the low 32 bits are
 * both stored as little-endian integers, but the high part is stored
 * first. So it is part big-endian, part little-endian.
 *
 * An offset of 0 means the corresponding brick or tile does not exist.
 * An offset of 1 means the brick or tile contains all zeros and does
 * not take any space on the file. An offset with the most significant
 * bit set also means the brick or tile has a constant value. In this
 * case the actual value is encoded in the least significant 8/16/32
 * bits (depending on valuetype) of the stored offset.
 * Offsets 0x8000000000000000 and 0x0000000000000001 are equivalent.
 * Actually, v2 and later uses the first form while v1 used the second.
 * For robustness both forms should be accepted regardless of version.
 *
 * Thread safety: None. The user of this type is responsible for that.
 */
class LookupTableV0Access : public ILookupTableAccess
{
private:
  std::vector<std::uint64_t> _lookup;  // This is the POD data from disk
  std::vector<std::uint64_t> _lookend; // Other data members are derived.
  bool _mustflip;
  bool _isalpha;
  //std::int64_t _lupsize;
public:
  explicit LookupTableV0Access(bool mustflip, bool isalpha, std::int64_t num_bricks = 0);
  podbytes_t podbytes() const override;
  void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) override;
  void byteswap() override;
  void dump(std::ostream& out, const std::string& prefix = "") override;
public:
  std::uint64_t lookupLinearIndex(std::int64_t index) const override;
  std::vector<std::uint64_t>& lup() override { return _lookup; }
  std::vector<std::uint64_t>& lupend() override { return _lookend; }
  const std::vector<std::uint64_t>& lup() const override { return _lookup; }
  const std::vector<std::uint64_t>& lupend() const override { return _lookend; }
};

LookupTableV0Access::LookupTableV0Access(bool mustflip, bool isalpha, std::int64_t num_bricks)
  : _lookup(num_bricks)
  , _lookend(num_bricks)
  , _mustflip(mustflip)
  , _isalpha(isalpha)
{
}

LookupTableV0Access::podbytes_t
LookupTableV0Access::podbytes() const
{
  std::int64_t bytesize = _lookup.size() * sizeof(std::uint64_t);
  IHeaderAccess::podbytes_t tmp(bytesize);
  memcpy(tmp.data(), _lookup.data(), bytesize);
  return tmp;
}

void
LookupTableV0Access::read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size)
{
  this->_lookup.clear();
  this->_lookup.resize(size/sizeof(std::int64_t), 0);
  file->xx_read(this->_lookup.data(), offset, size);
  byteswap();
  // TODO-Low can I get maxsize here as well?
  bool file_truncated = false;
  bool bricks_overlap = false;
  this->_lookend = LookupTable::calcLookupSize(this->_lookup, file->xx_eof(), 0,
                                               &file_truncated, &bricks_overlap);
  if (file_truncated)
    throw OpenZGY::Errors::ZgyFormatError("The ZGY file is truncated.");
  if (bricks_overlap)
    throw OpenZGY::Errors::ZgyFormatError("The ZGY file is corrupt. Bricks overlap.");
}

void
LookupTableV0Access::byteswap()
{
  if (_mustflip)
    byteswapV1Long(this->_lookup.data(), this->_lookup.size());
  else
    byteswapT(this->_lookup.data(), this->_lookup.size());
}

void
LookupTableV0Access::dump(std::ostream& out, const std::string& prefix)
{
  out << prefix << (_isalpha ? "alpha" : "brick")
      << " lookup table  "
      << this->_lookup.size() << " entries"
      << (_mustflip ? " (flipped)" : "")
      << "\n";
  for (size_t ii=0; ii<this->_lookup.size(); ++ii)
    out << prefix << std::hex
        << std::setw(16) << this->_lookup[ii] << " "
        << std::setw(16) << this->_lookend[ii] // - this->_lookup[ii]
        << std::dec << "\n";
}

std::uint64_t
LookupTableV0Access::lookupLinearIndex(std::int64_t index) const
{
  throw OpenZGY::Errors::ZgyInternalError("Not implemented");
}

std::shared_ptr<ILookupTableAccess>
HeaderAccessFactory::createAlphaLookup(std::uint32_t version)
{
  return std::shared_ptr<ILookupTableAccess>(new LookupTableV0Access(version==1, true));
}

std::shared_ptr<ILookupTableAccess>
HeaderAccessFactory::createBrickLookup(std::uint32_t version)
{
  return std::shared_ptr<ILookupTableAccess>(new LookupTableV0Access(version==1, false));
}

/////////////////////////////////////////////////////////////////////////////
//    Methods that logically belong earlier in this file but need to be   ///
//    here because of dependencies.                                       ///
/////////////////////////////////////////////////////////////////////////////

void
OffsetHeaderV1Access::calculate(const std::shared_ptr<IInfoHeaderAccess>& ih)
{
  _derived_infsize = sizeof(InfoHeaderV1POD);
  _derived_histsize = sizeof(HistHeaderV1POD);
  if (ih) {
    _derived_alphalupsize =
      calcEntriesInLut(ih->size(),ih->bricksize(),true) * sizeof(std::int64_t);
    _derived_bricklupsize =
      calcEntriesInLut(ih->size(),ih->bricksize(), false) * sizeof(std::int64_t);
  }
}

/**
 * Calculate offsets and sizes for the various headers and tables.
 * Some information requires the InfoHeader to be already known.
 * If it isn't we will just calculate as much as we can.
 *
 * In general the size of a header as written to file might be
 * larger than the size that the header expects to unpack.
 * This allows adding more data fields at the end of the header.
 * Older readers will just unpack the fields they know about.
 *
 * For ZGY V2 and V3 this is moot, as all the offsets are implicit
 * with all the headers written sequentially. So the size needs
 * to match exactly or the headers following this will be corrupt.
 */
void
OffsetHeaderV2Access::calculate(const std::shared_ptr<IInfoHeaderAccess>& ih)
{
  _derived_infoff = sizeof(FileHeaderPOD) + 1;
  _derived_infsize = sizeof(InfoHeaderV2POD);
  if (ih) {
    _derived_strsize = ih->slbufsize();
    _derived_histsize = sizeof(HistHeaderV2POD);
    _derived_alphalupsize =
      calcEntriesInLut(ih->size(), ih->bricksize(),true) * sizeof(std::int64_t);
    _derived_bricklupsize =
      calcEntriesInLut(ih->size(), ih->bricksize(), false) * sizeof(std::int64_t);
     _derived_stroff = _derived_infoff + _derived_infsize;
     _derived_histoff = _derived_stroff + _derived_strsize;
    _derived_alphalupoff = _derived_histoff + _derived_histsize;
    _derived_bricklupoff = _derived_alphalupoff + _derived_alphalupsize;
  }
}

/**
 * Calculate derived information that is too expensive to compute on
 * the fly. The code here is needed both after reading an existing
 * file and after creating a new one.
 */
void
InfoHeaderV1Access::calculate_cache()
{
  // See FileVersion<1>::InfoHeader::Interpret in the old library.

  // This information that cannot go stale because it only depends
  // on size and bricksize which are both immutable.
  // These are handled the same way in V1 and V2.
  this->_cached_lodsizes = calcLodSizes(this->size(),this->bricksize());
  this->_cached_nlods = static_cast<std::int32_t>(this->_cached_lodsizes.size());
  this->_cached_alphaoffsets = calcLutOffsets(this->_cached_lodsizes, true);
  this->_cached_brickoffsets = calcLutOffsets(this->_cached_lodsizes, false);
  std::tie(this->_cached_index, this->_cached_annot, this->_cached_world) =
    calcOrderedCorners(this->orig(),    this->inc(),    this->size(),
                       this->gpiline(), this->gpxline(),
                       this->gpx(),     this->gpy());
}

/**
 * Fix up information that logically belong in this header but is stored
 * elsewhere. This is not quite the same as calculate_cache() and it should
 * only be called after reading an existing fille.
 */
void
InfoHeaderV1Access::calculate_read(const podbytes_t& /*stringlist_in*/, const std::shared_ptr<IHistHeaderAccess>& hh)
{
  if (hh) {
    // Note reason for cast: ZGY only deals with int8, int16, float.
    _cached_file_sample_min = static_cast<float>(hh->minvalue());
    _cached_file_sample_max = static_cast<float>(hh->maxvalue());
  }
  else {
    _cached_file_sample_min = 0;
    _cached_file_sample_max = 0;
  }
  // Sanity check. If codingrange is bad, silently use a range
  // that causes no conversion between storage and float.
  _cached_safe_codingrange_min = _cached_file_sample_min;
  _cached_safe_codingrange_max = _cached_file_sample_max;
  _fix_codingrange(&_cached_safe_codingrange_min, &_cached_safe_codingrange_max, _pod._datatype);
  calculate_cache();
}

/**
 * Prepare to write this header to disk.
 * If this is a newly created header then the cached derived information
 * needs to be computed as well. It might be safer to do that unconditionally.
 */
InfoHeaderV1Access::podbytes_t
InfoHeaderV1Access::calculate_write()
{
  throw OpenZGY::Errors::ZgyInternalError("Writing is only supported for the latest version.");
}

/**
 * Calculate derived information that is too expensive to compute on
 * the fly. The code here is needed both after reading an existing
 * file and after creating a new one.
 */
void
InfoHeaderV2Access::calculate_cache()
{
  // More information that cannot go stale because it only depends
  // on size and bricksize which are both immutable.
  // These are handled the same way in V1 and V2.
  this->_cached_lodsizes = calcLodSizes(this->size(),this->bricksize());
  this->_cached_nlods = static_cast<std::int32_t>(this->_cached_lodsizes.size());
  this->_cached_alphaoffsets = calcLutOffsets(this->_cached_lodsizes, true);
  this->_cached_brickoffsets = calcLutOffsets(this->_cached_lodsizes, false);
  std::tie(this->_cached_index, this->_cached_annot, this->_cached_world) =
    calcOrderedCorners(this->orig(),    this->inc(),    this->size(),
                       this->gpiline(), this->gpxline(),
                       this->gpx(),     this->gpy());
}

/**
 * Fix up information that logically belong in this header but is stored
 * elsewhere. This is not quite the same as calculate_cache() and it should
 * only be called after reading an existing fille.
 */
void
InfoHeaderV2Access::calculate_read(const podbytes_t& stringlist_in, const std::shared_ptr<IHistHeaderAccess>& /*hh*/)
{
  // See FileVersion<2>::InfoHeader::Interpret in the old library.

  // Set the 5 strings if a StringHeader is provided, else clear them.
  // Normally we get called exactly once without a StringList, followed
  // by exactly one call that has the strings.
  // If strings are provided then there are supposed to be 5 consecutive
  // strings, all of them including the last should be null terminated.
  // Don't assume this is true, though.
  podbytes_t stringlist(stringlist_in);
  if (!stringlist.empty() && stringlist.back() != '\0')
    stringlist.push_back('\0');
  const unsigned char *cp = stringlist.data();
  const unsigned char *end = cp + stringlist.size();
  std::vector<std::string> strings;
  while (cp < end) {
    strings.push_back(std::string((char*)cp));
    cp += (strlen((char*)cp) + 1);
  }
  while (strings.size() < 5)
    strings.push_back(std::string());
  this->_srcname = strings[0];
  this->_srcdesc = strings[1];
  this->_hprjsys = strings[2];
  this->_hunitname = strings[3];
  this->_vunitname = strings[4];
  // Sanity check. If codingrange is bad, silently use a range
  // that causes no conversion between storage and float.
  _cached_safe_codingrange[0] = this->_pod._file_codingrange[0],
  _cached_safe_codingrange[1] = this->_pod._file_codingrange[1],
  _fix_codingrange(&this->_cached_safe_codingrange[0],
                   &this->_cached_safe_codingrange[1],
                   this->_pod._datatype);
  calculate_cache();
}

/**
 * Prepare for writing out this info header.
 *   1) Pack the 5 strings into a string header, which isn't a separate type.
 *   2) Record the size of the string header inside this info header.
 *   3) Return the bytes of the string header.
 *
 * If this is a newly created header then the cached derived information
 * needs to be computed as well. It might be safer to do that unconditionally.
 *
 * TODO-WARNING: If I at some point add a method to change the corner
 * coordinates of an existing file then either that function is responsible
 * for updating both the pod and the derved corners or this function
 * needs to compute one from the other. Depending on which was actually
 * changed.
 */
InfoHeaderV2Access::podbytes_t
InfoHeaderV2Access::calculate_write()
{
  std::string nullbyte(1,0);
  std::string stringlist =
    this->srcname() + nullbyte +
    this->srcdesc() + nullbyte +
    this->hprjsys() + nullbyte +
    this->hunitname() + nullbyte +
    this->vunitname() + nullbyte;
  podbytes_t slbuf(stringlist.begin(), stringlist.end());
  this->_pod._slbufsize = static_cast<std::uint32_t>(slbuf.size());
  calculate_cache();
  return slbuf;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

ZgyInternalMeta::ZgyInternalMeta()
  : _is_bad(false)
{
}

/**
 * Duplicated between impl/bulk.cpp and impl/meta.cpp but sets
 * different flags. See ZgyInternalMeta::ErrorsWillCorruptFile
 * for details, and ZgyWriter::errorflag().
 *
 * Thread safety:
 * The class itself is not thread safe but this should not be an issue,
 * because instances are meant to be short lived. Typically used inside
 * a method and not acceesible outside.
 */
class ZgyInternalMeta::ErrorsWillCorruptFile
{
  ZgyInternalMeta *_owner;
public:
  explicit ErrorsWillCorruptFile(ZgyInternalMeta* owner) : _owner(owner)
  {
  }
  ~ErrorsWillCorruptFile()
  {
    if (_owner)
      _owner->set_errorflag(true);
    _owner = nullptr;
  }
  void disarm()
  {
    _owner = nullptr;
  }
};

/**
 * Open for read
 */
ZgyInternalMeta::ZgyInternalMeta(const std::shared_ptr<IFileADT>& file, const LoggerFn& logger)
  : _loggerfn(logger ? logger : LoggerBase::standardCallback(LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"), "openzgy-meta: ", ""))
  , _is_bad(false)
{
  initFromOpenFile(file);
}

/**
 * Open for write
 */
ZgyInternalMeta::ZgyInternalMeta(const ZgyInternalWriterArgs& args_in, bool compress, const LoggerFn& logger)
  : _loggerfn(logger ? logger : LoggerBase::standardCallback(LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"), "openzgy-meta: ", ""))
  , _is_bad(false)
{
  initFromScratch(args_in, compress);
}

/**
 * Open for update
 */
ZgyInternalMeta::ZgyInternalMeta(const std::shared_ptr<IFileADT>& file, const ZgyInternalWriterArgs& args_in, bool compress, const LoggerFn& logger)
  : _loggerfn(logger ? logger : LoggerBase::standardCallback(LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"), "openzgy-meta: ", ""))
  , _is_bad(false)
{
  ZgyInternalWriterArgs args = validateOpenForUpdate(args_in);
  initFromOpenFile(file);
  if (this->_fh->version() < 3)
    throw OpenZGY::Errors::ZgyUpdateRules("File version is too old for this library to update");
  // No need to test for file too new, because that is already checked.
  initFromReopen(args, compress);

  constexpr float nan = std::numeric_limits<float>::quiet_NaN();
  if (args.datatype == RawDataType::Float32) {
    args.have_datarange = true;
    args.datarange = std::array<float,2>{nan, nan};
  }
}

/**
 * Convenience for invoking _loggerfn with a simple message.
 * This isn't much different from invoking the callback directly.
 * But it makes debugging slightly simpler to have an easy place
 * to set a breakpoint. It also adds more symmetry with respect
 * to the stringstream version, which does add value.
 */
bool
ZgyInternalMeta::_logger(int priority, const std::string& message) const
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
ZgyInternalMeta::_logger(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return _logger(priority, sstream ? sstream->str() : std::string());
}

/**
 * Read all the headers and save pointers to each of them in this instance.
 * Some care is needed to read them in the correct order, as there are
 * several dependencies between them. The next version of the file format
 * will hopefully be a bit simpler in this respect.
 *
 * Corresponds to ZgyInternalMeta._init_from_open_file() in Python.
 */
void
ZgyInternalMeta::initFromOpenFile(const std::shared_ptr<IFileADT>& file_in)
{
  // Wrapper to cache the first 256 KB
  std::shared_ptr<IFileADT> file(new FileWithSmallCache(file_in, 1024*256));

  // The file header contains only a magic string and the file version.
  // In V2 and V3 and V4 the headers are stored consecutively:
  //     FileHeader OffsetHeader InfoHeader StringList HistHeader
  //     AlphaLUT BrickLUT
  // In V1 and most likely in the upcoming V4 only FileHeader and
  // OffsetHeader are known to be consecutive.
  this->_fh = HeaderAccessFactory::createFileHeader();
  // TODO-Low, catch ZgyEndOfFile and map it to FormatError()
  // But EndOfFile usually means just that anyway. Maybe rephrase the message,
  // to "End of File, possibly a corrupted file, ..."
  this->_fh->read(file, 0);
  if (this->_fh->magic() == std::array<std::uint8_t,4>{'V','C','S',0})
    throw OpenZGY::Errors::ZgyNeedOldLibrary("Old ZGY compressed files are not supported in OpenZGY");
  else if (this->_fh->magic() != std::array<std::uint8_t,4>{'V','B','S',0})
    throw OpenZGY::Errors::ZgyFormatError("This is not a ZGY file");
  else if (this->_fh->version() < 1 || this->_fh->version() > 4)
    throw OpenZGY::Errors::ZgyFormatError("Unsupported ZGY version " +
                                          std::to_string(this->_fh->version()));

  int version = this->_fh->version();

  // The offset header immediately follows the file header.
  // Changed in v2 and v3: Offset header is no longer used
  // (all offsets are now implicit), but due to a quirk in the
  // implementation it still occupies one byte on the file.
  // This is actually a bug waiting to happen, because that
  // byte (which is the size of a class with no members) is
  // probably compiler dependant.
  // Removing the OffsetHeader actually made the files tricker
  // to read, as the size of some sections depend on contents
  // of other sections.
  this->_oh = HeaderAccessFactory::createOffsetHeader(version);
  this->_oh->read(file, 8);
  // Set infoff() and infsize(), the rest must wait.
  this->_oh->calculate(nullptr);

  // 'oh' at this point may be incomplete (V2, V3) but the offset
  // to the InfoHeader should be known by now. In V1 is is mostly
  // complete but is missing a few section sizes.
  this->_ih = HeaderAccessFactory::createInfoHeader(version);
  this->_ih->read(file, _oh->infoff(), _oh->infsize());

  // For V2 and V3, fill in the rest of the offsets now that
  // the InfoHeader is known.
  this->_oh->calculate(this->_ih);

  // Variable length strings are stored in a separate header.
  // Most (currently all) of these logically belong to InfoHeader.
  //sl = StringListFactory(fh._version).read(f, oh, ih);
  IHeaderAccess::podbytes_t sl(this->_ih->slbufsize());
  if (!sl.empty())
    file->xx_read(sl.data(), this->_oh->stroff(), sl.size());

  this->_hh = HeaderAccessFactory::createHistHeader(version);
  this->_hh->read(file, _oh->histoff(), _oh->histsize());
  this->_hh->calculate();

  // For V2 and later, fill in the missing strings in InfoHeader
  // now that the StringList is known.
  // For V1 we also need to copy the value range (used for scaling)
  // from the histogram header.
  this->_ih->calculate_read(sl, this->_hh);

  this->_alup = HeaderAccessFactory::createAlphaLookup(version);
  this->_blup = HeaderAccessFactory::createBrickLookup(version);

  this->_alup->read(file, _oh->alphalupoff(), _oh->alphalupsize());
  this->_blup->read(file, _oh->bricklupoff(), _oh->bricklupsize());
}

namespace {
  // Not efficient, but really portable and really not used much.
  static bool is_power_of_two(std::int64_t n)
  {
    for (int shift=0; shift<63; ++shift)
      if (n == static_cast<std::int64_t>(1) << shift)
        return true;
    return false;
  }
} // namespace

/**
 * Check that the arguments are consistent. Apply default values where
 * appropriate. The method is logically a part of the ZgyInternalMeta
 * constructor. The usual caveats regarding virtual methods apply.
 */
ZgyInternalWriterArgs
ZgyInternalMeta::validateCreateNewFile(const ZgyInternalWriterArgs& args_in)
{
  ZgyInternalWriterArgs args(args_in);
  if (args.size[0] <= 0 || args.size[1] <= 0 || args.size[2] <= 0)
    throw OpenZGY::Errors::ZgyUserError("size must be at least 1 in each dimension.");
  if (args.bricksize[0] <= 0 || args.bricksize[1] <= 0 || args.bricksize[2] <= 0)
    args.bricksize = std::array<std::int64_t,3>{64,64,64};
  // TEMPORARY. Officially there is no support for 2d data and it probably
  // doesn't work so this should only be used for experiments.
  std::int64_t minsize =
    Environment::getNumericEnv("OPENZGY_ALLOW_2D", 0) > 0 ? 1 : 4;
  for (std::size_t ii=0; ii<args.bricksize.size(); ++ii)
    if (args.bricksize[ii] < minsize || !is_power_of_two(args.bricksize[ii]))
      throw OpenZGY::Errors::ZgyUserError("bricksize must be >= " + std::to_string(minsize) + " and a power of 2.");
  constexpr std::int64_t maxsize = std::numeric_limits<std::int32_t>::max();
  for (std::size_t ii=0; ii<args.bricksize.size(); ++ii)
    if (args.bricksize[ii] > maxsize || args.size[ii] > maxsize)
      throw OpenZGY::Errors::ZgyUserError("size must be < " + std::to_string(maxsize));

  // The codingrange for floating point data is special. The user is
  // not allowed to set it, and its value is not used internally.
  // To ensure that any user supplied range is really ignored we set
  // the range to NaN. In _calculate_write it will be set to the
  // statistical range before being written out. As an additional
  // bullet-proofing, to avoid surprises with older files, this
  // rule can also be enforced in api.ZgyMeta.datarange.
  // Note: For integral types I might have defaulted the datarange
  // to no conversion (-128..+127 or -32768..+32767) and also
  // silently re-ordered min and max if needed. But in both cases
  // the application is buggy. So, make a fuss.
  // A data range for integral data covering just a single value
  // (presumably the user will just write that single value to
  // the file) is also forbidden because it just isn't useful
  // and it triggers several corner cases.
  constexpr float nan = std::numeric_limits<float>::quiet_NaN();
  if (args.datatype == RawDataType::Float32)
    args.datarange = std::array<float,2>{nan, nan};
  else if (args.datarange[0] >= args.datarange[1])
    throw OpenZGY::Errors::ZgyUserError("datarange must be specified and have min < max.");
  else if (!std::isfinite(args.datarange[0]) || !std::isfinite(args.datarange[1]))
    throw OpenZGY::Errors::ZgyUserError("datarange must be finite");
  if (!args.zunitfactor) args.zunitfactor = 1.0;
  if (!args.hunitfactor) args.hunitfactor = 1.0;
  return args;
}

/**
 * Check arguments for a file opened for update.
 *
 *  - size, bricksize, datatype cannot change.
 *
 *  - datarange cannot change for an integral file and is
 *    ignored for float files so update is blocked in both cases.
 *
 *  - Horizontal and vertical units are problematic because
 *    changing the unit name may resize the string table.
 *    Can probably be made to work but it probably isn't
 *    very useful.
 *
 *  - Explicitly specifying any of these values is disallowed,
 *    even if the old and new values match.
 *    To make things easier for the application it would be
 *    possible to define a ZgyUpdateArgs as a supertype of
 *    ZgyWriterArgs where the banned attributes are absent.
 *    I don't think I'll bother.
 */
ZgyInternalWriterArgs
ZgyInternalMeta::validateOpenForUpdate(const ZgyInternalWriterArgs& args_in)
{
  ZgyInternalWriterArgs args(args_in);
  if (args.have_size      ||
      args.have_bricksize ||
      args.have_datatype  ||
      args.have_datarange ||
      args.have_hunit     ||
      args.have_zunit) {
    throw OpenZGY::Errors::ZgyUserError("Cannot change size, bricksize, datatype, datarange or units of an existing file");
  }
  if (args.zunitfactor == 0)
    args.zunitfactor = 1;
  if (args.hunitfactor == 0)
    args.hunitfactor = 1;
  return args;
}

/**
 * Create a new ZGY file in memory with all metadata filled in and no
 * bulk written yet.
 *
 * This method is logically a part of the ZgyInternalMeta constructor.
 * The usual caveats regarding virtual methods apply.
 *
 * Corresponds to ZgyInternalMeta._init_from_scratch() in Python.
 * Arguments here are passed as a struct, since very long argument
 * lists are a problem for languages that don't allow keyword arguments.
 *
 * The complete list of ZgyInternalWriterArgs we need to handle:
 *
 *     (filename), size, bricksize, datatype, datarange,
 *     zunitdim, hunitdim, zunitname, hunitname, zunitfactor, hunitfactor,
 *     zstart, zinc, annotstart, annotinc, corners.
 *
 * We also need to know whether ZFP or other compression is in use,
 * which means that the user provided a compressor or we are updating
 * a file already flagged as potentially compressed (i.e. version >= 4)
 *
 * The api currently does not allow setting hprjsys, srcname, srcdesc.
 * If this is desired then those need to be added as well.
 */
void
ZgyInternalMeta::initFromScratch(const ZgyInternalWriterArgs& args_in, bool compress)
{
  if (_logger(2))
    _logger(2, std::stringstream()
            << "\nCreating new file \"" << args_in.filename << "\""
            << (compress ? " with compression" : "")
            << "\n");

  ZgyInternalWriterArgs args = validateCreateNewFile(args_in);

  std::shared_ptr<FileHeaderAccess> fh(new FileHeaderAccess());
  std::uint8_t magic[4]{'V', 'B', 'S', 0};
  memcpy(fh->_pod._magic, magic, 4);
  fh->_pod._version = compress ? 4 : 3;

  std::shared_ptr<OffsetHeaderV2Access> oh(new OffsetHeaderV2Access());

  std::shared_ptr<InfoHeaderV2Access> ih(new InfoHeaderV2Access());
  InfoHeaderV2POD& pod(ih->_pod);
  // I will try to fill in everything below; this is just for robustness.
  memset(&pod, 0, sizeof(pod));

  // Fill in the parts of InfoHeader that affects headers elsewhere.
  // Validate checks that args did not exceed INT_MAX, which is a
  // limitation in current ZGY. With that check the casts are safe.
  pod._bricksize[0] = static_cast<std::int32_t>(args.bricksize[0]);
  pod._bricksize[1] = static_cast<std::int32_t>(args.bricksize[1]);
  pod._bricksize[2] = static_cast<std::int32_t>(args.bricksize[2]);
  pod._size[0]      = static_cast<std::int32_t>(args.size[0]);
  pod._size[1]      = static_cast<std::int32_t>(args.size[1]);
  pod._size[2]      = static_cast<std::int32_t>(args.size[2]);
  pod._datatype     = static_cast<std::uint8_t>(args.datatype);
  pod._slbufsize    = 0;

  // Meta information caller is allowed to specify.
  // Naming is not 100% consistent; this is because
  // the parameters to this function wree set to match
  // the existing Python wrapper for the old ZGY.

  pod._file_codingrange[0] = args.datarange[0];
  pod._file_codingrange[1] = args.datarange[1];
  ih->_cached_safe_codingrange[0] = pod._file_codingrange[0];
  ih->_cached_safe_codingrange[1] = pod._file_codingrange[1];

  pod._vdim = static_cast<std::uint8_t>(args.zunitdim);
  pod._hdim = static_cast<std::uint8_t>(args.hunitdim);

  pod._vunitfactor = args.zunitfactor;
  pod._hunitfactor = args.hunitfactor;

  pod._orig[0] = args.annotstart[0];
  pod._orig[1] = args.annotstart[1];
  pod._orig[2] = args.zstart;

  pod._inc[0] = args.annotinc[0];
  pod._inc[1] = args.annotinc[1];
  pod._inc[2] = args.zinc;

  // The annotation corners stored in _gpiline and _gpxline are redundant
  // when using the four-point method but the ZGY file format requires
  // them to be written and also uses them on read.
  float beg[2] {pod._orig[0], pod._orig[1]};
  float end[2] {pod._orig[0] + pod._inc[0] * (pod._size[0]-1),
                pod._orig[1] + pod._inc[1] * (pod._size[1]-1)};
  float acorners[4][2] {
      {beg[0], beg[1]},
      {end[0], beg[1]},
      {beg[0], end[1]},
      {end[0], end[1]}};
  for (size_t ii=0; ii<4; ++ii) {
    pod._gpx[ii] = args.corners[ii][0];
    pod._gpy[ii] = args.corners[ii][1];
    pod._gpiline[ii] = acorners[ii][0];
    pod._gpxline[ii] = acorners[ii][1];
  }

  // Meta information that might be updated after creation.
  // Except for dataid.

  GUID().copyTo(&pod._dataid[0], sizeof(pod._dataid));
  GUID().copyTo(&pod._verid[0],  sizeof(pod._verid));
  GUID(nullptr).copyTo(&pod._previd[0], sizeof(pod._previd));

  if (_logger(5))
    _logger(5, std::stringstream()
            << "dataid: "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._dataid[0])).toString()
            << "\nverid:  "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._verid[0])).toString()
            << "\nprevid: "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._previd[0])).toString()
            << "\n");

  pod._srctype = pod._datatype;

  // Statistics.
  pod._scnt = 0;
  pod._ssum = 0.0;
  pod._sssq = 0.0;
  pod._smin = 0.0;
  pod._smax = 0.0;

  // Unused fields, required to be set this way for compatibility.
  for (int ii=0; ii<3; ++ii) {
    pod._curorig[ii] = 0;
    pod._cursize[ii] = pod._size[ii];
    pod._srvorig[ii] = pod._orig[ii];
    pod._srvsize[ii] = pod._size[ii] * pod._inc[ii];
  }
  pod._gdef = (std::uint8_t)RawGridDefinition::FourPoint;
  pod._gazim[0] = 0.0;
  pod._gazim[1] = 0.0;
  pod._gbinsz[0] = 0.0;
  pod._gbinsz[1] = 0.0;

  // Set derived members that cannot change later.
  //    _cached_lodsizes, _cached_nlods, _cached_{alpha,brick}offsets
  //    _cached_{index,annot,world}
  // TODO-Worry: codingrange for float data cannot be set until finalize,
  // make sure I remember to do that. Also, should I reset
  // _cached_sample_min and _cached_sample_max here, just in case?

  // Note that this will clear the 5 attributes from the string table.
  // So make sure this is done *before* setting those.

  // STRING TABLE. These logically belong in ih->_pod but since they are
  // variable length strings they are put in a separate section.
  // Variable length also means they cannot be changed (or at least made
  // longer) once set. Even if we might support updates in general, later.
  ih->_srcname = "";
  ih->_srcdesc = "";
  ih->_hprjsys = "";
  ih->_hunitname = args.hunitname;
  ih->_vunitname = args.zunitname;

  // We don't need the stringlist header until flushing, but some
  // derived fields may need to be set.
  IHeaderAccess::podbytes_t discard = ih->calculate_write();

  // Fill in the rest of the offsets now that the InfoHeader is known.
  // Requires _slbufsize, ih._alphaoffsets, _brickoffsets
  oh->calculate(ih);

  // Histogram gets initialized to empty
  std::shared_ptr<HistHeaderV2Access> hh(new HistHeaderV2Access());
  // Belt & suspenders (the constructor also does this)
  memset(&hh->_pod, 0, sizeof(hh->_pod));

  // Lookup tables get initialized to empty
  std::shared_ptr<LookupTableV0Access> alup
    (new LookupTableV0Access(false, true,  ih->_cached_alphaoffsets.back()));
  std::shared_ptr<LookupTableV0Access> blup
    (new LookupTableV0Access(false, false, ih->_cached_brickoffsets.back()));

  this->_fh = fh;
  this->_oh = oh;
  this->_ih = ih;
  this->_hh = hh;
  this->_alup = alup;
  this->_blup = blup;

  // If "assert contig" in _flush_meta fails, this is also wrong:
  // ZGY aligns data to the basic block size, which depends on
  // data type and bricksize. This simplifies caching data.
  // If all headers are written sequentially from the start
  // of the file then it is simpler to add the padding as part
  // of the header data instead of before the first data block.
  // Forget about using leftover space in the header to store the
  // first few alpha tiles. We probably won't be writing those anyway.

  //code = impl.enum._map_DataTypeToStructFormatCode(self._ih._datatype);
  //bs = self._ih._bricksize;
  //bs = bs[0] * bs[1] * bs[2] * struct.calcsize(code);
  //hdrsize = self._oh._bricklupoff + self._oh._bricklupsize;
  //padsize = (((hdrsize+bs-1)/bs)*bs)-hdrsize;
  //self._data_beg = hdrsize + padsize;
  //self._data_end = hdrsize + padsize;
}

void
ZgyInternalMeta::initFromReopen(const ZgyInternalWriterArgs& args_in, bool compress)
{
  if (_logger(2))
    _logger(2, std::stringstream()
            << "\nOpen for update \"" << args_in.filename << "\""
            << (compress ? " with compression" : "")
            << "\n");

  ZgyInternalWriterArgs args = args_in;

  // All the information the user is allowed to update is in the InfoHeader.
  // And there has already been a check on the file version number so the
  // cast below should not throw. But use a dynamic cast just to be safe.

  InfoHeaderV2Access& ih(dynamic_cast<InfoHeaderV2Access&>(*this->_ih));
  InfoHeaderV2POD& pod(ih._pod);

  // TODO-Low: Refactor: It is possible to partly consolidate initFromScratch()
  // and initFromReopen() if I implement ZgyInternalWriterArgs::merge()
  // so I can fill in args with values from file for each of the arguments
  // the user doesn't want to override. Probably doesn't make much difference.

  if (args.have_zunit) {
    pod._vdim = static_cast<std::uint8_t>(args.zunitdim);
    pod._vunitfactor = args.zunitfactor;
  }

  if (args.have_hunit) {
    pod._hdim = static_cast<std::uint8_t>(args.hunitdim);
    pod._hunitfactor = args.hunitfactor;
  }

  if (args.have_ilstart)
    pod._orig[0] = args.annotstart[0];
  if (args.have_xlstart)
    pod._orig[1] = args.annotstart[1];
  if (args.have_zstart)
    pod._orig[2] = args.zstart;

  if (args.have_ilinc)
    pod._inc[0] = args.annotinc[0];
  if (args.have_xlinc)
    pod._inc[1] = args.annotinc[1];
  if (args.have_zinc)
    pod._inc[2] = args.zinc;

  if (args.have_corners) {
    float beg[2] {pod._orig[0], pod._orig[1]};
    float end[2] {pod._orig[0] + pod._inc[0] * (pod._size[0]-1),
        pod._orig[1] + pod._inc[1] * (pod._size[1]-1)};
    float acorners[4][2] {
      {beg[0], beg[1]},
      {end[0], beg[1]},
      {beg[0], end[1]},
      {end[0], end[1]}};
    for (size_t ii=0; ii<4; ++ii) {
      pod._gpx[ii] = args.corners[ii][0];
      pod._gpy[ii] = args.corners[ii][1];
      pod._gpiline[ii] = acorners[ii][0];
      pod._gpxline[ii] = acorners[ii][1];
    }
  }

  // Unused fields, required to be set this way for compatibility.
  // COPY/PASTE from ZgyInternalMeta::initFromScratch().
  for (int ii=0; ii<3; ++ii) {
    pod._curorig[ii] = 0;
    pod._cursize[ii] = pod._size[ii];
    pod._srvorig[ii] = pod._orig[ii];
    pod._srvsize[ii] = pod._size[ii] * pod._inc[ii];
  }
  pod._gdef = (std::uint8_t)RawGridDefinition::FourPoint;
  pod._gazim[0] = 0.0;
  pod._gazim[1] = 0.0;
  pod._gbinsz[0] = 0.0;
  pod._gbinsz[1] = 0.0;
  // End COPY/PASTE.

  memcpy(&pod._previd[0], &pod._verid[0], sizeof(pod._previd));
  GUID().copyTo(&pod._verid[0], sizeof(pod._verid));

  if (_logger(5))
    _logger(5, std::stringstream()
            << "dataid: "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._dataid[0])).toString()
            << "\nverid:  "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._verid[0])).toString()
            << "\nprevid: "
            << GUID(ptr_to_array<std::uint8_t,16>(&pod._previd[0])).toString()
            << "\n");

  // calculate_read and calculate_cache have already been called when
  // opening the file. As far as I can see there is very little that
  // needs to be updated here, and nothing from calculate_write which
  // mostly deals with the string table. Which currently doesn't allow
  // changes. Copy/paste what is needed from calculate_cache. Yes this
  // is a bad code smell. My excuse: The rest of calculate_cache ought
  // to be a no-op but it is safer to not try to touch it.
  std::tie(ih._cached_index, ih._cached_annot, ih._cached_world) =
    calcOrderedCorners(ih.orig(),    ih.inc(),    ih.size(),
                       ih.gpiline(), ih.gpxline(),
                       ih.gpx(),     ih.gpy());

  if (!has_finalized_feature()) {
    // Clear histogram, statistics, and low resolution data because
    // everything will be rebuilt on finalize. The only thing kept
    // is the previous histogram range. See ZgyInternalBulk ctor.
    // Which hasn't been run yet, so take care not to clobber all.
    // Code in ZgyInternalBulk() and elsewhere won't need the
    // convoluted test for havelods==nlods; it just looks at _scnt.
    std::shared_ptr<HistHeaderV2Access> hh(new HistHeaderV2Access());
    memset(&hh->_pod, 0, sizeof(hh->_pod));
    if (this->_hh->minvalue() <= this->_hh->maxvalue()) {
      // TODO-Worry: These casts to hide compiler warnings MUST be
      // removed once the file format is changed to use wider (double)
      // types. Otherwise the added accuracy won't help.
      hh->_pod._min = static_cast<float>(this->_hh->minvalue());
      hh->_pod._max = static_cast<float>(this->_hh->maxvalue());
      this->_hh = hh;
    }
    pod._scnt = 0;
    pod._ssum = 0;
    pod._sssq = 0;
    pod._smin = 0;
    pod._smax = 0;
    LookupTable::setNoBrickLOD
      (ih.lodsizes(), ih.brickoffsets(),
       &this->_blup->lup(), &this->_blup->lupend());
  }

  if (!can_append_bulk(compress)) {
    // Disallow opening a finalized file if it was or will be compressed.
    // Because we must assume that the application is going to finalize
    // in this session and that would not be allowed, and the attempt
    // would leave the file corrupted. To be really pedantic I could
    // allow the compressor to be set as long as the lodcompressor is
    // not, but this is getting ridiculous.
    // TODO-@@@-Medium: It should still be possible to open a compressed
    // file for the sole purpose of changing annotation and world coords.
    // That would not require a finalize. So ideally the test should be
    // deferred to the first write. But there are other problems such as
    // also deferring the re-open segment.
    if (compress)
      throw OpenZGY::Errors::ZgyUserError
        ("A finalized file cannot have compressed data appended.");
    else
      throw OpenZGY::Errors::ZgyUserError
        ("A finalized compressed file cannot be opened for update.");
  }
}

/**
 * Store metadata on file and return the number of bytes written,
 * including padding. Pass a null file for a dry run that only reports
 * the size. The dry run may be used as a consistency check when
 * opening a file on the cloud for update.
 */
std::int64_t
ZgyInternalMeta::flushMeta(const std::shared_ptr<IFileADT>& file)
{
  std::int64_t bytes_written{0};
  const std::int64_t alignto = this->_ih->bytesperbrick();

  IHeaderAccess::podbytes_t slbuf = this->_ih->calculate_write(); // sets _pod._slbufsize
  this->_oh->calculate(this->_ih); // offsets change with slbuf.size()

  std::vector<IHeaderAccess::podbytes_t> bytes
    { this->_fh->podbytes(),
      this->_oh->podbytes(),
      this->_ih->podbytes(),
      slbuf,
      this->_hh->podbytes(),
      this->_alup->podbytes(),
      this->_blup->podbytes() };

  if (_logger(5)) {
    std::stringstream ss;
    this->_oh->dump(ss);
    _logger(5, ss);
  }

  if ((std::int64_t)bytes[2].size() != this->_oh->infsize() ||
      //(std::int64_t)bytes[3].size() != this->_oh->strsize() ||
      (std::int64_t)bytes[3].size() != this->_ih->slbufsize() ||
      (std::int64_t)bytes[4].size() != this->_oh->histsize() ||
      (std::int64_t)bytes[5].size() != this->_oh->alphalupsize() ||
      (std::int64_t)bytes[6].size() != this->_oh->bricklupsize())
    throw OpenZGY::Errors::ZgyInternalError("Header size mismatch on write");

  const IOffsetHeaderAccess& oh = *this->_oh;

  if ( (oh.infoff()      == (std::int64_t)(bytes[0].size() + bytes[1].size()) &&
        oh.stroff()      == oh.infoff() + oh.infsize() &&
        //oh.histoff()     == oh.stroff() + oh.strsize() &&
        oh.alphalupoff() == oh.histoff() + oh.histsize() &&
        oh.bricklupoff() == oh.alphalupoff() + oh.alphalupsize())) {
    // All headers are written sequentially from the start of the file.
    // In this situation it is simpler to add any required padding between
    // the last header and the (alignment needed) first data brick to the
    // end of the header. Especially if we are writing to the cloud.
    // Forget about using leftover space in the header to store the
    // first few alpha tiles. We probably won't be writing those anyway.
    // TODO-Low: Should the padding be removed for the compressed case?
    // Not strictly needed there, but if writing directly to the cloud
    // the padding ensures that all segments have nice sizes.
    // This can make life easier for a future OpenZGY with a cache module
    // when it wants to read the file that we are creating now.
    IHeaderAccess::podbytes_t allbytes;
    for (const auto& it : bytes)
      std::copy(it.begin(), it.end(), std::back_inserter(allbytes));
    const std::int64_t paddedsize = ((allbytes.size() + alignto - 1) / alignto) * alignto;
    allbytes.resize(paddedsize, 0);
    if (file) {
      ErrorsWillCorruptFile watchdog(this);
      file->xx_write(allbytes.data(), 0, allbytes.size(), UsageHint::Header);
      watchdog.disarm();
    }
    bytes_written = allbytes.size();
  }
  else {
    // Headers are not contiguous. This should never happen for
    // uncompressed V2 and V3. So I won't add code for that case
    // (which I will never be able to test).
    throw OpenZGY::Errors::ZgyInternalError("Headers not contiguous on write");
  }
  return bytes_written;
}

void
ZgyInternalMeta::dump(std::ostream& out, const std::string& prefix)
{
  out << prefix << "FileHeader" << (_fh ? "" : " (null)") << "\n";
  if (_fh) _fh->dump(out, prefix + "  ");

  out << prefix << "OffsetHeader" << (_oh ? "" : " (null)") << "\n";
  if (_oh) _oh->dump(out, prefix + "  ");

  out << prefix << "InfoHeader" << (_ih ? "" : " (null)") << "\n";
  if (_ih) _ih->dump(out, prefix + "  ");

  out << prefix << "HistHeader" << (_hh ? "" : " (null)") << "\n";
  if (_hh) _hh->dump(out, prefix + "  ");

  out << prefix << "AlphaLup" << (_alup ? "" : " (null)") << "\n";
  if (_alup) _alup->dump(out, prefix + "  ");

  out << prefix << "BrickLup" << (_blup ? "" : " (null)") << "\n";
  if (_blup) _blup->dump(out, prefix + "  ");
}

/*
 * Caveat: Adding bloat to this class. Find a better design?
 *
 * The has_xxx_feature() and can_xxx() tests will respectively
 * collect information about the file and inform the caller about
 * what is allowed.
 *
 * Initially there were just "was/is compressed" and "was finalized"
 * cheks but things get more complicated when a file can be finalized
 * with respect to low resolution data but not statistics or vice
 * versa. Or the file may be small enough to not even have low
 * resolution bricks.
 *
 * The tests that determine the version number are similar to the tests
 * that decide whether incremental builds are allowed and whether a
 * compressed file may be re-opened.
 *
 * Why does this matter?
 *
 * Petrel creates an empty file that will be written to later. This
 * file must not be finalized because that would prevent it from being
 * later written with compressed data. And if the finalize was done
 * with compression it would prevent even ubcompressed writes. If no
 * compression at all is involved then a finalize at this point "only"
 * a serious performance issue.
 */

/**
 * \brief Statistics and histogram are good.
 *
 * \details
 * Good statistics and histogram has both sane limits (min<max) and
 * nonzero counts. Note a few implementation details: An empty
 * histogram might have limits (0,0) if never updated, and it might
 * have (min<max) coupled with zero count because even if the file
 * didn't get finalized it might need to know what the range would
 * have been.
 */
bool
ZgyInternalMeta::has_stathist_feature() const
{
  return (_hh->samplecount() != 0 &&
          _ih->scnt() != 0 &&
          _ih->smin() <= _ih->smax());
}

/**
 * \brief Entire file is just a single brick.
 *
 * \details
 * A survey small enough to fit inside a single brick triggers several
 * corner cases. Such as never having any low resolution bricls.
 */
bool
ZgyInternalMeta::has_tiny_feature() const
{
  const std::int64_t nlods = static_cast<std::int64_t>(_ih->lodsizes().size());
  return nlods == 1;
}

/**
 * \brief Has low resolution. Always false for single-brick files.
 *
 * \details
 * Returns true if finalize has been run and low resolution bricks
 * have been stored. For tiny files where all samples fit into a
 * single brick, i.e. has_tiny_feature() == true, this will return
 * false because the file still won't have (and won't need) lowres.
 */
bool
ZgyInternalMeta::has_lowres_feature() const
{
  const std::int64_t nlods = static_cast<std::int64_t>(_ih->lodsizes().size());
  const std::int64_t havelods = LookupTable::usableBrickLOD
    (_ih->lodsizes(), _ih->brickoffsets(), _blup->lup(), _blup->lupend());
  return nlods == havelods && nlods > 1;
}

/**
 * \brief Contains compressed bricks.
 *
 * \details
 * TODO-@@@-Low: Performance: Might need to cache the result.
 * In that case, marking the cache dirty later on is a hassle.
 */
bool
ZgyInternalMeta::has_compression_feature() const
{
  return InternalZGY::LookupTable::hasBrickCompression
    (_blup->lup(), _blup->lupend());
}

/**
 * \brief File is fully finalized.
 *
 * \details
 * Both low resolution bricks (if needed), statistics, and histogram
 * are present. If this test fails we might as well remove all the
 * results of a previous finalize.
 */
bool
ZgyInternalMeta::has_finalized_feature() const
{
  return (has_stathist_feature() &&
          (has_tiny_feature() || has_lowres_feature()));
}

/**
 * \brief Readable by the old library.
 *
 * \details
 * When OpenZGY saves a file that isn't finalized then it is flagged
 * as v4 because the old reader cannot handle that. The same happens
 * if the file contains compressed data.
 *
 * If histogram or statistics are missing due to the file not being
 * finalized then this does not in itself prevent the old library
 * from using it. It is assumed that application code will react
 * properly when it sees a range with min>max or count==0.
 *
 * A tiny uncompressed file that was not finalized will still be
 * usable by ZGY-Public because it never has lowres.
 *
 * The version number (3 or 4) should reflect the current state of the
 * file, not its history. So e.g. a file that was v4 only due to
 * missing low resolution data should be changed to v3 when finalized.
 * The test thus needs to be done right before closing the file.
 */
bool
ZgyInternalMeta::can_old_library_read() const
{
  return (!has_compression_feature() &&
          (has_tiny_feature() || has_lowres_feature()));
}

/**
 * \brief Eligible for incremental finalize.
 *
 * \details
 * The test is similar to can_old_library_read(). If low resolution data
 * is missing then not only can't the file be opened by old zgy, it can't
 * be incrementally finalized either because there is no starting point.
 * The test here is stricter because even if only statistics and/or
 * histogram are bad then incremental finalize os still disallowed.
 *
 * The user can still elect to do a full finalize. In fact, that is
 * currently the default.
 *
 * Compressed files can't be finalized more than once. Nor can
 * uncompressed files if the latest reopen asked for compression.
 */
bool
ZgyInternalMeta::can_finalize_incremental() const
{
  return !has_compression_feature() && has_finalized_feature();
}

/**
 * \brief Allowed to append (not necessarily update) bulk data.
 *
 * \details
 * Compressed and finalized files are not allowed to have data appended
 * to them. Not because the write itself would be a problem, but because
 * re-finalizing the file counts an an update not an append so that would
 * not be allowed.
 *
 * In the special single-brick case, appending is technically allowed
 * even for compressed finalized data because no low resolution bricks
 * need to be updated. This is of academic interest only (or it might
 * change the wording of error messages) because a compressed file
 * must have at least one brick of real data. And since the entire
 * file is just one brick, there are no more empty bricks. So there
 * is no space to append anything anyway.
 *
 * If the file is currently uncompressed but this current session will
 * be adding compressed data then it still needs to be considered to
 * be compressed. The caller needs to provide the will_compress flag
 * because ibky the accessor knows what it is planning to do.
 *
 * Caveat: Ideally the test should be made on the first write and not
 * the file reopen, because maybe the application just wanted to update
 * the metadata. Which ought to still be allowed. But that special case
 * also triggers other issues. How to prevent (for cloud access) the
 * last segment from being re-opened needlessly.
 */
bool
ZgyInternalMeta::can_append_bulk(bool will_compress) const
{
  return ((!has_compression_feature() && !will_compress) ||
          (has_tiny_feature() || !has_lowres_feature()));
}

} // namespace
