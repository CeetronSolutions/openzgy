// Copyright 2017-2024, Schlumberger
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

#include "../declspec.h"
#include "enum.h" // Need full definition of struct InternalLodMode.

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <array>
#include <limits>

namespace OpenZGY {
  class ZgyWriterArgs;
  class ZgyWriterArgsV2;
  class ZgyWriterArgsV3;
  enum class DecimationType;
  class IOContext;
  class IZgyWriter;
}

namespace OpenZGY { namespace Impl {
  class ZgyWriter;
}}

namespace InternalZGY {
  enum class LodAlgorithm;
  class ZgyInternalMeta;
}

namespace Test {
  class ZgyWriterMock;
}

namespace InternalZGY {
#if 0
}
#endif

/////////////////////////////////////////////////////////////////////////////
/// LEGACY TYPES DO NOT USE -- See V3 version below.                      ///
/////////////////////////////////////////////////////////////////////////////

/**
 * There are no setters in the Impl class. The owner (an instance of
 * OpenZGY::ZgyWriterArgsV2 instance) is a friend and is allowed to
 * set data members directly. Inlined and/or virtual methods
 * are not allowed in the API facing class.
 *
 * There are getters and they are not in-line, i.e. they need to be
 * both declared here and implemented in writerargsimpl.h. This is in
 * addition to the user visible setter that is declared in writerargs.h
 * and implemented in writerargs.cpp.
 *
 * Applications should treat the argument package as write only.
 * Applications should not access OpenZGY::ZgyWriterArgsV2::impl()
 * in spite of it being public. Maybe it shouldn't be. But even if
 * it does, it will return an opaque type unless "impl/writeargsimpl.h"
 * is included. And that file is not in the SDK.
 *
 * TODO-WIP-BrickedAPI: Improve data hiding. Eventually the Impl class
 * could provide getters for all the parameters. The user visible contructor
 * in api.cpp could then pass only the Impl class down to the lower
 * levels. As an interim solution, without breaking the ABI, the "old"
 * parameters would be retrieved from the parent that this class stores
 * a reference to. Later, all members should be moved inside the pimpl.
 * See also the discussion in iocontextimpl.h.
 *
 * Code at the impl level should not access the OpenZGY::ZgyWriterArgsV2
 * instance, and should treat the Impl class as read-only.
 */
class OPENZGY_TEST_API ZgyWriterArgsV2Impl
{
  friend OpenZGY::ZgyWriterArgsV2;
  friend OpenZGY::ZgyWriterArgsV3;

protected:
  int version_; // Concrete type is ZgyWriterArgsV{version_}Impl

private:
  InternalLodMode internal_lod_mode_;
  std::array<float,2> historange_;
  // Not clear what is cleanest, storing the API facing type or the internal.
  std::vector<OpenZGY::DecimationType> decimation_;
  std::vector<InternalZGY::LodAlgorithm> lod_algorithm_;
  std::string dataid_; // Normally set implicitly on file creation.
  std::string verid_;  // Normally set implicitly on file creation and update.
  std::string previd_; // Normally copied from verid when a new verid is set.

  bool has_lod_lowlevel_;
  bool has_lod_lodmode_;
  bool has_lod_environ_;
  bool has_dataid_;
  bool has_verid_;
  bool has_previd_;

public:
  InternalLodMode getInternalLodMode() const;
  std::array<float,2> getHistoRange() const;
  std::vector<OpenZGY::DecimationType> getDecimation() const;
  std::vector<InternalZGY::LodAlgorithm> getLodAlgorithm() const;
  std::string getDataId() const;
  std::string getVerId() const;
  std::string getPrevId() const;

public:
  explicit ZgyWriterArgsV2Impl();
  ZgyWriterArgsV2Impl(const ZgyWriterArgsV2Impl&);
  ZgyWriterArgsV2Impl(const ZgyWriterArgsV2Impl&&) = delete;
  const ZgyWriterArgsV2Impl& operator=(const ZgyWriterArgsV2Impl&) = delete;
  const ZgyWriterArgsV2Impl& operator=(const ZgyWriterArgsV2Impl&&) = delete;
};

/////////////////////////////////////////////////////////////////////////////
/// CURRENT TYPES                                                         ///
/////////////////////////////////////////////////////////////////////////////

/**
 * \brief PImpl pattern for OpenZGY::ZgyWriterArgsV3.
 * \details
 *
 * Note that there are only setters in the API visible class
 * (because it should be write-only seen from the applicaton)
 * and no setters in the implementation (owner is a friend)
 * and no getters in the implementation class (because I am lazy).

 * Due to my laziness, members in the implementation will need to be
 * public. This is simple to fix but might not be worse the additional
 * cost of maintaining the class.
 *
 * The list of members is from the old classes ZgyWriterArgs in api.h
 * and class ZgyInternalWriterArgs in meta.h. Plus what was added in
 * V2, class ZgyWriterArgsV2Impl in this file.
 *
 * Inlined and/or virtual methods are not allowed in the API facing
 * class. I won't make that mistake again.
 *
 * Applications should treat the argument package as write only.
 * Applications should not access OpenZGY::ZgyWriterArgsV2::impl() in
 * spite of it being public. Maybe it shouldn't be. But even if
 * accessed from the application, it will return an opaque type unless
 * "impl/writeargsimpl.h" is included. And that file is not in the SDK.
 *
 * Code at the impl level should not access the OpenZGY::ZgyWriterArgsV3
 * instance, and should ideally treat the Impl class as read-only. As
 * with enum types, code in the api level (api.cpp or writerargs.h)
 * knows how to create one of these instances from a ZgyWriterArgsV3.
 * Conversion in the other direction is not useful.
 *
 * This is a short lived helper class for passing arguments to the
 * functions that create a ZGY file. It uses the named parameter
 * idiom. This is needed because there are a lot of arguments.
 * Do NOT use this class for holding on to the information after
 * open() or reopen() returns to the application that called it.
 *
 * Compared to ZgyWriterArgs there were some members missing here
 * because they are processed at a higher level.
 *
 * @(#) That might needs to change now that the user visible class doesn't store a copy.
 *
 *  - template
 *  - iocontext
 *  - compressor
 *  - lodcompressor
 *
 * @(#) The following enums need to be mapped from api types to internal.
 *
 *  - datatype
 *  - zunitdim
 *  - hunitdim
 *  - internal_lod_mode
 *  - decimation
 *  - lod_algorithm_
 *
 * The have_xxx members keep track of what information has been
 * changed from the default. Currently metafrom() will set all to
 * true. That might change. The information is only used by merge().
 * Transient information such as iocontext is ignored.

 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to open().
 */
class ZgyWriterArgsV3Impl
{
public:
  typedef std::array<std::int64_t,3> index3_t;
  typedef double float64_t;
  typedef std::array<std::array<float64_t,2>,4> corners_t;
  typedef std::pair<std::shared_ptr<const void>, std::int64_t> rawdata_t;
  typedef std::function<rawdata_t(const rawdata_t&,const std::array<int64_t,3>&)> compressor_t;
  friend class OpenZGY::ZgyWriterArgsV3;
  // Grant read/write access to data members for all *internal* code that
  // needs it. Could simply have made the data members public. But doing
  // it this way makes it harder for the application code to gain access.
  friend class OpenZGY::Impl::ZgyWriter;
  friend class OpenZGY::IZgyWriter;
  friend class InternalZGY::ZgyInternalMeta;
  friend class Test::ZgyWriterMock;

protected:
  int version_; // Concrete type is ZgyWriterArgsV{version_}Impl

private:
  std::string            filename_;
  std::shared_ptr<const OpenZGY::IOContext> iocontext_;
  compressor_t           compressor_;
  compressor_t           lodcompressor_;
  index3_t               size_;
  index3_t               bricksize_;
  RawDataType            datatype_;
  std::array<float,2>    datarange_;
  RawVerticalDimension   zunitdim_;
  RawHorizontalDimension hunitdim_;
  std::string            zunitname_, hunitname_;
  double                 zunitfactor_, hunitfactor_;
  float                  zstart_, zinc_;
  std::array<float,2>    annotstart_, annotinc_;
  corners_t              corners_;
  // New in V2
  InternalLodMode        internal_lod_mode_;
  std::array<float,2>    historange_;
  std::vector<OpenZGY::DecimationType> decimation_;
  std::vector<InternalZGY::LodAlgorithm> lod_algorithm_;
  std::string dataid_; // Normally set implicitly on file creation.
  std::string verid_;  // Normally set implicitly on file creation and update.
  std::string previd_; // Normally copied from verid when a new verid is set.
  // Old
  bool have_size_;
  bool have_bricksize_;
  bool have_datatype_;
  bool have_datarange_;
  bool have_zunit_;
  bool have_hunit_;
  bool have_ilstart_;
  bool have_ilinc_;
  bool have_xlstart_;
  bool have_xlinc_;
  bool have_zstart_;
  bool have_zinc_;
  bool have_corners_;
  // New in V2
  bool has_lod_lowlevel_;
  bool has_lod_lodmode_;
  bool has_lod_environ_;
  bool has_dataid_;
  bool has_verid_;
  bool has_previd_;

public:
  ZgyWriterArgsV3Impl()
    : version_(3)
    , filename_("")
    , iocontext_(nullptr)
    , compressor_()
    , lodcompressor_()
    , size_{0,0,0}
    , bricksize_{64,64,64}
    , datatype_(RawDataType::Float32)
    , datarange_{0, -1}
    , zunitdim_(RawVerticalDimension::Unknown)
    , hunitdim_(RawHorizontalDimension::Unknown)
    , zunitname_("")
    , hunitname_("")
    , zunitfactor_(1.0)
    , hunitfactor_(1.0)
    , zstart_(0)
    , zinc_(0)
    , annotstart_{0,0}
    , annotinc_{0,0}
    , corners_{{{0,0},{0,0},{0,0},{0,0}}}
      // New in V2
    , internal_lod_mode_{}
    , historange_(std::array<float,2>
                  {std::numeric_limits<float>::infinity(),
                   -std::numeric_limits<float>::infinity()})
    , decimation_()
    , lod_algorithm_()
    , dataid_()
    , verid_()
    , previd_()
      // Old
    , have_size_(false)
    , have_bricksize_(false)
    , have_datatype_(false)
    , have_datarange_(false)
    , have_zunit_(false)
    , have_hunit_(false)
    , have_ilstart_(false)
    , have_ilinc_(false)
    , have_xlstart_(false)
    , have_xlinc_(false)
    , have_zstart_(false)
    , have_zinc_(false)
    , have_corners_(false)
      // New in V2
    , has_lod_lowlevel_(false)
    , has_lod_lodmode_(false)
    , has_lod_environ_(false)
    , has_dataid_(false)
    , has_verid_(false)
    , has_previd_(false)
  {
    // Default to LodMode::Early
    internal_lod_mode_.incr.level = -1;
    internal_lod_mode_.incr.force = 1;
    internal_lod_mode_.last.level = -1;
    internal_lod_mode_.last.force = 2;
  }

  void merge(const ZgyWriterArgsV3Impl& other);
};

} // namespace
