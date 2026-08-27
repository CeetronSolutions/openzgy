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

#include "writerargsimpl.h"
#include "environment.h"
#include <string>
#include <cstdint>
#include <iostream>
#include <limits>

namespace InternalZGY {
#if 0
}
#endif

/////////////////////////////////////////////////////////////////////////////
/// LEGACY TYPES DO NOT USE -- See V3 version below.                      ///
/////////////////////////////////////////////////////////////////////////////

ZgyWriterArgsV2Impl::ZgyWriterArgsV2Impl()
  : version_(2)
  , internal_lod_mode_{}
  , historange_(std::array<float,2>
                {std::numeric_limits<float>::infinity(),
                 -std::numeric_limits<float>::infinity()})
  , decimation_()
  , lod_algorithm_()
  , dataid_()
  , verid_()
  , previd_()
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

ZgyWriterArgsV2Impl::ZgyWriterArgsV2Impl(const ZgyWriterArgsV2Impl& other)
  : version_(2)
  , internal_lod_mode_(other.internal_lod_mode_)
  , historange_(other.historange_)
  , decimation_(other.decimation_)
  , lod_algorithm_(other.lod_algorithm_)
  , dataid_(other.dataid_)
  , verid_(other.verid_)
  , previd_(other.previd_)
  , has_lod_lowlevel_(other.has_lod_lowlevel_)
  , has_lod_lodmode_(other.has_lod_lodmode_)
  , has_lod_environ_(other.has_lod_environ_)
  , has_dataid_(other.has_dataid_)
  , has_verid_(other.has_verid_)
  , has_previd_(other.has_previd_)
{
}

InternalLodMode
ZgyWriterArgsV2Impl::getInternalLodMode() const
{
  return internal_lod_mode_;
}

std::array<float,2>
ZgyWriterArgsV2Impl::getHistoRange() const
{
  return historange_;
}

std::vector<OpenZGY::DecimationType>
ZgyWriterArgsV2Impl::getDecimation() const
{
  return decimation_;
}

std::vector<InternalZGY::LodAlgorithm>
ZgyWriterArgsV2Impl::getLodAlgorithm() const
{
  return lod_algorithm_;
}

std::string
ZgyWriterArgsV2Impl::getDataId() const
{
  return dataid_;
}

std::string
ZgyWriterArgsV2Impl::getVerId() const
{
  return verid_;
}

std::string
ZgyWriterArgsV2Impl::getPrevId() const
{
  return previd_;
}

/////////////////////////////////////////////////////////////////////////////
/// CURRENT TYPES                                                         ///
/////////////////////////////////////////////////////////////////////////////

/**
 * Update *this with settings from rhs, but only those values that
 * were explicitly set. Compared to the version in ZgyWriterArgsV2
 * this one is more consistent. Including filename, iocontext,
 * and the two compressors. Neither of those are likely to be useful,
 * but it is safer to do everything.
 *
 * Do not use in conjunction with metafrom(), because it is then
 * unspecified what happens with the have_xxx members.
 *
 * TODO: Yagni?
 * This was meant to be used when opening for update, but that code
 * now looks directly at the have_xxx members.
 */
void
ZgyWriterArgsV3Impl::merge(const ZgyWriterArgsV3Impl& rhs)
{
  if (!rhs.filename_.empty()) {
    this->filename_ = rhs.filename_;
  }
  if (rhs.iocontext_) {
    this->iocontext_ = rhs.iocontext_;
  }
  if (rhs.compressor_) {
    this->compressor_ = rhs.compressor_;
  }
  if (rhs.lodcompressor_) {
    this->lodcompressor_ = rhs.lodcompressor_;
  }
  if (rhs.have_size_) {
    this->size_ = rhs.size_;
    this->have_size_ = rhs.have_size_;
  }
  if (rhs.have_bricksize_) {
    this->bricksize_ = rhs.bricksize_;
    this->have_bricksize_ = rhs.have_bricksize_;
  }
  if (rhs.have_datatype_) {
    this->datatype_ = rhs.datatype_;
    this->have_datatype_ = rhs.have_datatype_;
  }
  if (rhs.have_datarange_) {
    this->datarange_ = rhs.datarange_;
    this->have_datarange_ = rhs.have_datarange_;
  }
  if (rhs.have_zunit_) {
    this->zunitdim_    = rhs.zunitdim_;
    this->zunitname_   = rhs.zunitname_;
    this->zunitfactor_ = rhs.zunitfactor_;
    this->have_zunit_  = rhs.have_zunit_;
  }
  if (rhs.have_hunit_) {
    this->hunitdim_    = rhs.hunitdim_;
    this->hunitname_   = rhs.hunitname_;
    this->hunitfactor_ = rhs.hunitfactor_;
    this->have_hunit_  = rhs.have_hunit_;
  }
  if (rhs.have_ilstart_) {
    this->annotstart_[0] = rhs.annotstart_[0];
    this->have_ilstart_ = rhs.have_ilstart_;
  }
  if (rhs.have_ilinc_) {
    this->annotinc_[0] = rhs.annotinc_[0];
    this->have_ilinc_ = rhs.have_ilinc_;
  }
  if (rhs.have_xlstart_) {
    this->annotstart_[1] = rhs.annotstart_[1];
    this->have_xlstart_ = rhs.have_xlstart_;
  }
  if (rhs.have_xlinc_) {
    this->annotinc_[1] = rhs.annotinc_[1];
    this->have_xlinc_ = rhs.have_xlinc_;
  }
  if (rhs.have_zstart_) {
    this->zstart_ = rhs.zstart_;
    this->have_zstart_ = rhs.have_zstart_;
  }
  if (rhs.have_zinc_) {
    this->zinc_ = rhs.zinc_;
    this->have_zinc_ = rhs.have_zinc_;
  }
  if (rhs.have_corners_) {
    this->corners_ = rhs.corners_;
    this->have_corners_ = rhs.have_corners_;
  }
  if (rhs.has_lod_lowlevel_ || rhs.has_lod_lodmode_) {
    this->internal_lod_mode_ = rhs.internal_lod_mode_;
    this->has_lod_lowlevel_  = rhs.has_lod_lowlevel_;
    this->has_lod_lodmode_   = rhs.has_lod_lodmode_;
  }
  if (rhs.historange_[0] < rhs.historange_[1]) {
    this->historange_ = rhs.historange_;
  }
  if (!rhs.decimation_.empty()) {
    this->decimation_ = rhs.decimation_;
  }
  if (!rhs.lod_algorithm_.empty()) {
    this->lod_algorithm_ = rhs.lod_algorithm_;
  }
  if (rhs.has_dataid_) {
    this->dataid_ = rhs.dataid_;
    this->has_dataid_ = rhs.has_dataid_;
  }
  if (rhs.has_verid_) {
    this->verid_ = rhs.verid_;
    this->has_verid_ = rhs.has_verid_;
  }
  if (rhs.has_previd_) {
    this->previd_ = rhs.previd_;
    this->has_previd_ = rhs.has_previd_;
  }
}

} // namespace
