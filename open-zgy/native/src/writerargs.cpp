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

#include "writerargs.h"
#include "iocontext.h"
#include "impl/writerargsimpl.h"
#include "impl/environment.h"
#include "impl/lodalgo.h"
#include "impl/compression.h"
#include "impl/enum_mapper.h"

#include <sstream>

namespace OpenZGY {
#if 0
}
#endif

/////////////////////////////////////////////////////////////////////////////
/// LEGACY TYPES DO NOT USE -- See V3 version below.                      ///
/////////////////////////////////////////////////////////////////////////////

ZgyWriterArgsV2::ZgyWriterArgsV2()
  : ZgyWriterArgs()
  , pimpl_(std::make_shared<Impl>())
{
  setFromEnv();
}

ZgyWriterArgsV2::~ZgyWriterArgsV2()
{
}

ZgyWriterArgsV2::ZgyWriterArgsV2(const ZgyWriterArgs& other, FromOldWriterArgs)
  : ZgyWriterArgs(other)
  , pimpl_(std::make_shared<Impl>())
{
  setFromEnv();
}

ZgyWriterArgsV2::ZgyWriterArgsV2(const ZgyWriterArgsV2& other)
  : ZgyWriterArgs(other)
  , pimpl_(std::make_shared<Impl>(*other.pimpl_))
{
}

const ZgyWriterArgsV2::Impl&
ZgyWriterArgsV2::impl() const
{
  return *pimpl_;
}

/**
 * For testing. Override any LodMode requested by the application.
 * Unlike most other environment variables in OpenZGY, this is not
 * a fallback. It has precedence over settings in code.
 */
void
ZgyWriterArgsV2::setFromEnv()
{
  using InternalZGY::Environment;

  std::string mode = Environment::getStringEnv("OPENZGY_LODMODE");
  if (!mode.empty()) {
    // TODO-WIP-BrickedAPI: Use logger.
    std::cerr << "Overriding LodMode defaults from environment."
              << " LodMode=" << mode << "\"\n" << std::flush;
    if (mode == "Default")
      lodmode(LodMode::Default);
    else if (mode == "Early")
      lodmode(LodMode::Early);
    else if (mode == "Early1")
      lodmode(LodMode::Early1);
    else if (mode == "Never")
      lodmode(LodMode::Never);
    else if (mode == "Rebuild")
      lodmode(LodMode::Rebuild);
    else
      throw Errors::ZgyUserError
        ("Invalid value $OPENZGY_LODMODE=\"" + mode + "\"");
    pimpl_->has_lod_environ_ = true;
  }

  // For testing only. Setting these to meaningless values can have
  // seriously bad effects. Parameter combinations that make sense
  // should be added to the LodMode enum.
  // The ability to directly set the 4 low level parameters might
  // be removed without warning, both by environment and by
  // setting in the argument package.
  if (!Environment::getStringEnv("OPENZGY_LOD_INCR_FORCE").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_INCR_LEVEL").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_LAST_FORCE").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_LAST_LEVEL").empty())
  {
    pimpl_->internal_lod_mode_.incr.force = Environment::getNumericEnv
      ("OPENZGY_LOD_INCR_FORCE", pimpl_->internal_lod_mode_.incr.force);
    pimpl_->internal_lod_mode_.incr.level = Environment::getNumericEnv
      ("OPENZGY_LOD_INCR_LEVEL", pimpl_->internal_lod_mode_.incr.level);
    pimpl_->internal_lod_mode_.last.force = Environment::getNumericEnv
      ("OPENZGY_LOD_LAST_FORCE", pimpl_->internal_lod_mode_.last.force);
    pimpl_->internal_lod_mode_.last.level = Environment::getNumericEnv
      ("OPENZGY_LOD_LAST_LEVEL", pimpl_->internal_lod_mode_.last.level);
    std::cerr << "Overriding LodMode defaults from environment."
              << "\nINCR_FORCE=" << pimpl_->internal_lod_mode_.incr.force
              << " INCR_LEVEL=" << pimpl_->internal_lod_mode_.incr.level
              << " LAST_FORCE=" << pimpl_->internal_lod_mode_.last.force
              << " LAST_LEVEL=" << pimpl_->internal_lod_mode_.last.level
              << std::endl;
    pimpl_->has_lod_environ_ = true;
  }
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodmode(LodMode mode)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lowlevel_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lodmode_ = true;

  switch (mode)
  {
  case LodMode::Default:
  case LodMode::Early:
    pimpl_->internal_lod_mode_.incr.force = 1;
    pimpl_->internal_lod_mode_.incr.level = -1;
    pimpl_->internal_lod_mode_.last.force = 2;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  case LodMode::Early1:
    pimpl_->internal_lod_mode_.incr.force = 1;
    pimpl_->internal_lod_mode_.incr.level = 1;
    pimpl_->internal_lod_mode_.last.force = 2;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  case LodMode::Never:
    pimpl_->internal_lod_mode_.incr.force = 0;
    pimpl_->internal_lod_mode_.incr.level = 0;
    pimpl_->internal_lod_mode_.last.force = 0;
    pimpl_->internal_lod_mode_.last.level = 0;
    break;

  case LodMode::Rebuild:
    pimpl_->internal_lod_mode_.incr.force = 0;
    pimpl_->internal_lod_mode_.incr.level = 0;
    pimpl_->internal_lod_mode_.last.force = 3;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  default:
    throw Errors::ZgyUserError("Invalid value for lodmode");
  }
  return *this;
}


ZgyWriterArgsV2&
ZgyWriterArgsV2::lodIncrForce(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  switch (n)
  {
  case 0:
  case 1:
    pimpl_->internal_lod_mode_.incr.force = n;
    break;
  default:
    throw Errors::ZgyUserError("Invalid value for lodIncrForce");
  }
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodIncrLevel(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  pimpl_->internal_lod_mode_.incr.level = n;
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodLastForce(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  switch (n)
  {
  case 0:
  case 2:
  case 3:
    pimpl_->internal_lod_mode_.last.force = n;
    break;
  case 1: // Might be supported in the future.
  default:
    throw Errors::ZgyUserError("Invalid value for lodLastForce");
  }
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodLastLevel(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  pimpl_->internal_lod_mode_.last.level = n;
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::historange(float lo, float hi)
{
  pimpl_->historange_ = std::array<float,2>{lo,hi};
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::decimation(const std::vector<DecimationType> decimation)
{
  pimpl_->decimation_ = decimation;
  // The call to EnumMapper::mapDecimationTypeToLodAlgorithm()
  // is currently done from ZgyWriter::enableIncrementalLOD().
  // TODO: Now that EnumMapper has been refactored into its own file,
  // it should be possible to store the internal type in both
  // the V2 and the V3 pimpl. Until then, pimpl_->lod_algorithm_
  // should really have been removed so it won't get used by
  // mistake.
  //std::vector<InternalZGY::LodAlgorithm> algorithms =
  //  OpenZGY::Impl::EnumMapper::mapDecimationTypeToLodAlgorithm(decimation);
  //pimpl_->lod_algorithm_ = algorithms;
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::guidsfrom(const std::shared_ptr<OpenZGY::IZgyReader>& reader)
{
  pimpl_->verid_     = reader->verid();
  pimpl_->has_verid_ = true;

  // dataid() and previd() are missing from the IZgyReader API.
  // Adding them would technically be a breaking change.
  // But I want to maintain the information in the file.
  //
  // It is tempting to access the impl layer, like ZgyMeta::verid() does.
  // But IZgyReader cannot access the ZgyReader implementation. And there
  // might even be a chain of IZgyReader instances. None of which will
  // propagate dataid / previd or the raw file header. Not to mention
  // that including impl/meta.h in this file is a rather gross breach
  // of encapsulation. Would also need impl/guid.h, but that one is safe.
#if 0
  pimpl_->dataid_ = InternalZGY::GUID(reader->ih().dataid()).toString();
  pimpl_->previd_ = InternalZGY::GUID(reader->ih().previd()).toString();
  pimpl_->has_dataid_ = true;
  pimpl_->has_previd_ = true;
#endif
  return *this;
}

void
ZgyWriterArgsV2::dump(std::ostream& out) const
{
  ZgyWriterArgs::dump(out);
  out << "ZgyWriterArgsV2\n"
      << "  lowres:     "
      << " incr (" << impl().getInternalLodMode().incr.level
      << ", " << impl().getInternalLodMode().incr.force << ")"
      << " last (" << impl().getInternalLodMode().last.level
      << ", " << impl().getInternalLodMode().last.force << ")\n"
      << "  historange:  " << impl().getHistoRange()[0] << " to"
      << " " << impl().getHistoRange()[1] << "\n";
  out << "  decimation: ";
  for (const auto& it : impl().getDecimation())
    out << " " << (int)it;
  out << "\n";
  out << "  algorithm:  ";
  for (const auto& it : impl().getLodAlgorithm())
    out << " " << (int)it;
  out << "\n";
  if ((!impl().getDataId().empty()) ||
      (!impl().getVerId().empty()) ||
      (!impl().getPrevId().empty())) {
    out << "  dataid:      " << impl().getDataId() << "\n"
        << "  verid:       " << impl().getVerId()  << "\n"
        << "  previd:      " << impl().getPrevId() << "\n";
  }
}

// Forward to the old arg package, to get the correct (V2) struct
// returned. So the named parameter idiom still works.
// CAVEAT: Missing a function here will introduce serious bugs.
// The wrong open() or reopen() overload will get called, causing
// all the settings in ZgyWriterArgsV2 to reset to default.
// Tip: Have the last setting in the change be a V2 parameter
// if possible. If any of the earlier functions were missing
// then this will cause a compiler error.

ZgyWriterArgsV2&
ZgyWriterArgsV2::filename(const std::string& value)
{
  ZgyWriterArgs::filename(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::iocontext(const IOContext *value)
{
  ZgyWriterArgs::iocontext(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::compressor(const compressor_t& value)
{
  ZgyWriterArgs::compressor(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::compressor(const std::string& name, const std::vector<std::string>& args)
{
  ZgyWriterArgs::compressor(name, args);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::zfp_compressor(float value)
{
  ZgyWriterArgs::zfp_compressor(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodcompressor(const compressor_t& value)
{
  ZgyWriterArgs::lodcompressor(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::lodcompressor(const std::string& name, const std::vector<std::string>& args)
{
  ZgyWriterArgs::lodcompressor(name, args);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::zfp_lodcompressor(float value)
{
  ZgyWriterArgs::zfp_lodcompressor(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::size(std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  ZgyWriterArgs::size(ni, nj, nk);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::bricksize(std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  ZgyWriterArgs::bricksize(ni, nj, nk);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::datatype(SampleDataType value)
{
  ZgyWriterArgs::datatype(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::datarange(float lo, float hi)
{
  ZgyWriterArgs::datarange(lo, hi);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::zunit(UnitDimension dimension, const std::string& name, double factor)
{
  ZgyWriterArgs::zunit(dimension, name, factor);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::hunit(UnitDimension dimension, const std::string& name, double factor)
{
  ZgyWriterArgs::hunit(dimension, name, factor);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::ilstart(float value)
{
  ZgyWriterArgs::ilstart(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::ilinc(float value)
{
  ZgyWriterArgs::ilinc(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::xlstart(float value)
{
  ZgyWriterArgs::xlstart(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::xlinc(float value)
{
  ZgyWriterArgs::xlinc(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::zstart(float value)
{
  ZgyWriterArgs::zstart(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::zinc(float value)
{
  ZgyWriterArgs::zinc(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::corners(const corners_t& value)
{
  ZgyWriterArgs::corners(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::metafrom(const std::shared_ptr<OpenZGY::IZgyReader>& value)
{
  ZgyWriterArgs::metafrom(value);
  return *this;
}

ZgyWriterArgsV2&
ZgyWriterArgsV2::merge(const ZgyWriterArgsV2& value)
{
  ZgyWriterArgs::merge(value);
  // Handle newer atributes owned by the V2 structs.
  // TODO: Does this belong in the impl layer?
  ZgyWriterArgsV2::Impl& self(*this->pimpl_);
  const ZgyWriterArgsV2::Impl& other(value.impl());
  if (other.has_dataid_)
    self.dataid_ = other.dataid_;
  if (other.has_verid_)
    self.verid_ = other.verid_;
  if (other.has_previd_)
    self.previd_ = other.previd_;
  return *this;
}

/////////////////////////////////////////////////////////////////////////////
/// CURRENT TYPES                                                         ///
/////////////////////////////////////////////////////////////////////////////

ZgyWriterArgsV3::ZgyWriterArgsV3()
  : pimpl_(std::make_shared<Impl>())
{
  setFromEnv();
}

ZgyWriterArgsV3::~ZgyWriterArgsV3()
{
}

ZgyWriterArgsV3::ZgyWriterArgsV3(const ZgyWriterArgs& other, FromOldWriterArgs)
  : pimpl_(std::make_shared<Impl>())
{
  copyFromV1(other);
  setFromEnv();
}

ZgyWriterArgsV3::ZgyWriterArgsV3(const ZgyWriterArgsV2& other, FromOldWriterArgsV2)
  : pimpl_(std::make_shared<Impl>())
{
  copyFromV2(other);
}

ZgyWriterArgsV3::ZgyWriterArgsV3(const ZgyWriterArgsV3& rhs)
  : pimpl_(std::make_shared<Impl>(*rhs.pimpl_))
{
}

ZgyWriterArgsV3::ZgyWriterArgsV3(ZgyWriterArgsV3&& rhs)
  : pimpl_(nullptr)
{
  std::swap(pimpl_, rhs.pimpl_);
}

const ZgyWriterArgsV3&
ZgyWriterArgsV3::operator=(const ZgyWriterArgsV3& rhs)
{
  pimpl_ = std::make_shared<Impl>(*rhs.pimpl_);
  return *this;
}

const ZgyWriterArgsV3&
ZgyWriterArgsV3::operator=(ZgyWriterArgsV3&& rhs)
{
  if (this != &rhs) {
    std::swap(pimpl_, rhs.pimpl_);
    rhs.pimpl_.reset();
  }
  return *this;
}

/**
 * Fill as many parameters as possible from an instance of one of the older
 * argument package. Also copy the have_xxx information used for merging.
 * Paremeters that did not exist in the old package will be left unset.
 *
 * Implementation notes: The "old" parameter is a ZgyWriterArgs, not a
 * ZgyInternalWriterArgs from meta.h. The latter is deprecated and might
 * even be removed without affecting the ABI.
 *
 * Reading from the old instance must access data members using friend access.
 * Writing to *this can either use the public setters or access data members
 * in *pimpl_ directly.
 *
 * Note that using the setters will update have_xxx while using the pimpl
 * will not. As long as this code copies all the have_xxx at the end,
 * then this should not make any difference.
 */
ZgyWriterArgsV3&
ZgyWriterArgsV3::copyFromV1(const ZgyWriterArgs& old)
{
  Impl& pi(*this->pimpl_);

  pi.filename_       = old._filename;
  pi.iocontext_      = old._iocontext; // The pointer will be shared.
  pi.compressor_     = old._compressor;
  pi.lodcompressor_  = old._lodcompressor;
  pi.size_           = old._size;
  pi.bricksize_      = old._bricksize;
#if 1
  this->datatype(old._datatype); // Use setter to get enum conversion.
  pi.datarange_      = old._datarange;
  this->zunit(old._zunitdim, old._zunitname, old._zunitfactor); // enum convert.
  this->hunit(old._hunitdim, old._hunitname, old._hunitfactor); // enum convert.
#else
  pi.datatype_       = old._datatype;
  pi.datarange_      = old._datarange;
  pi.zunitdim_       = old._zunitdim;
  pi.zunitname_      = old._zunitname;
  pi.zunitfactor_    = old._zunitfactor;
  pi.hunitdim_       = old._hunitdim;
  pi.hunitname_      = old._hunitname;
  pi.hunitfactor_    = old._hunitfactor;
#endif
  pi.zstart_         = old._zstart;
  pi.zinc_           = old._zinc;
  pi.annotstart_     = old._annotstart;
  pi.annotinc_       = old._annotinc;
  pi.corners_        = old._corners;
  pi.have_size_      = old._have_size;
  pi.have_bricksize_ = old._have_bricksize;
  pi.have_datatype_  = old._have_datatype;
  pi.have_datarange_ = old._have_datarange;
  pi.have_zunit_     = old._have_zunit;
  pi.have_hunit_     = old._have_hunit;
  pi.have_ilstart_   = old._have_ilstart;
  pi.have_ilinc_     = old._have_ilinc;
  pi.have_xlstart_   = old._have_xlstart;
  pi.have_xlinc_     = old._have_xlinc;
  pi.have_zstart_    = old._have_zstart;
  pi.have_zinc_      = old._have_zinc;
  pi.have_corners_   = old._have_corners;

  // These were not introduced until V2.
  //lodmode(LodMode);
  //lodIncrForce(int);
  //lodIncrLevel(int);
  //lodLastForce(int);
  //lodLastLevel(int);
  //historange(float lo, float hi);
  //decimation(const std::vector<DecimationType>);

  // Special, N/A for copying individual parameters.
  //metafrom(const std::shared_ptr<OpenZGY::IZgyReader>&);
  //guidsfrom(const std::shared_ptr<OpenZGY::IZgyReader>&);
  //merge(const ZgyWriterArgsV3&);
  //void dump(std::ostream&) const;

  return *this;
}

/**
 * Fill as many parameters as possible from an instance of one of the older
 * argument package. Also copy the have_xxx information used for merging.
 * Paremeters that did not exist in the old package will be left unset.
 *
 * See copyFromV1 for more details.
 *
 * Copying the V2 data members here is simpler than using the setters,
 * especially for LodMode which has three ways of setting the info.
 */
ZgyWriterArgsV3&
ZgyWriterArgsV3::copyFromV2(const ZgyWriterArgsV2& old)
{
  copyFromV1(old);

  Impl& pi(*this->pimpl_);
  const ZgyWriterArgsV2::Impl& oi(old.impl());

  pi.internal_lod_mode_ = oi.internal_lod_mode_;
  pi.historange_        = oi.historange_;
  pi.decimation_        = oi.decimation_;
  pi.lod_algorithm_     = oi.lod_algorithm_; // Currently unused -> decimation.
  pi.dataid_            = oi.dataid_;
  pi.verid_             = oi.verid_;
  pi.previd_            = oi.previd_;

  pi.has_lod_lowlevel_ = oi.has_lod_lowlevel_;
  pi.has_lod_lodmode_  = oi.has_lod_lodmode_;
  pi.has_lod_environ_  = oi.has_lod_environ_;
  pi.has_dataid_       = oi.has_dataid_;
  pi.has_verid_        = oi.has_verid_;
  pi.has_previd_       = oi.has_previd_;

  return *this;
}

const ZgyWriterArgsV3::Impl&
ZgyWriterArgsV3::impl() const
{
  return *pimpl_;
}

/**
 * For testing. Override any LodMode requested by the application.
 * Unlike most other environment variables in OpenZGY, this is not
 * a fallback. It has precedence over settings in code.
 */
void
ZgyWriterArgsV3::setFromEnv()
{
  using InternalZGY::Environment;

  std::string mode = Environment::getStringEnv("OPENZGY_LODMODE");
  if (!mode.empty()) {
    // TODO-WIP-BrickedAPI: Use logger.
    std::cerr << "Overriding LodMode defaults from environment."
              << " LodMode=" << mode << "\"\n" << std::flush;
    if (mode == "Default")
      lodmode(LodMode::Default);
    else if (mode == "Early")
      lodmode(LodMode::Early);
    else if (mode == "Early1")
      lodmode(LodMode::Early1);
    else if (mode == "Never")
      lodmode(LodMode::Never);
    else if (mode == "Rebuild")
      lodmode(LodMode::Rebuild);
    else
      throw Errors::ZgyUserError
        ("Invalid value $OPENZGY_LODMODE=\"" + mode + "\"");
    pimpl_->has_lod_environ_ = true;
  }

  // For testing only. Setting these to meaningless values can have
  // seriously bad effects. Parameter combinations that make sense
  // should be added to the LodMode enum.
  // The ability to directly set the 4 low level parameters might
  // be removed without warning, both by environment and by
  // setting in the argument package.
  if (!Environment::getStringEnv("OPENZGY_LOD_INCR_FORCE").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_INCR_LEVEL").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_LAST_FORCE").empty() ||
      !Environment::getStringEnv("OPENZGY_LOD_LAST_LEVEL").empty())
  {
    pimpl_->internal_lod_mode_.incr.force = Environment::getNumericEnv
      ("OPENZGY_LOD_INCR_FORCE", pimpl_->internal_lod_mode_.incr.force);
    pimpl_->internal_lod_mode_.incr.level = Environment::getNumericEnv
      ("OPENZGY_LOD_INCR_LEVEL", pimpl_->internal_lod_mode_.incr.level);
    pimpl_->internal_lod_mode_.last.force = Environment::getNumericEnv
      ("OPENZGY_LOD_LAST_FORCE", pimpl_->internal_lod_mode_.last.force);
    pimpl_->internal_lod_mode_.last.level = Environment::getNumericEnv
      ("OPENZGY_LOD_LAST_LEVEL", pimpl_->internal_lod_mode_.last.level);
    std::cerr << "Overriding LodMode defaults from environment."
              << "\nINCR_FORCE=" << pimpl_->internal_lod_mode_.incr.force
              << " INCR_LEVEL=" << pimpl_->internal_lod_mode_.incr.level
              << " LAST_FORCE=" << pimpl_->internal_lod_mode_.last.force
              << " LAST_LEVEL=" << pimpl_->internal_lod_mode_.last.level
              << std::endl;
    pimpl_->has_lod_environ_ = true;
  }
}

/////////////////////////////////////////////////////////////////////////////
///   NAMED PARAMETERS INTRODUCED IN ZgyWriterArgsV2                      ///
/////////////////////////////////////////////////////////////////////////////

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodmode(LodMode mode)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lowlevel_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lodmode_ = true;

  switch (mode)
  {
  case LodMode::Default:
  case LodMode::Early:
    pimpl_->internal_lod_mode_.incr.force = 1;
    pimpl_->internal_lod_mode_.incr.level = -1;
    pimpl_->internal_lod_mode_.last.force = 2;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  case LodMode::Early1:
    pimpl_->internal_lod_mode_.incr.force = 1;
    pimpl_->internal_lod_mode_.incr.level = 1;
    pimpl_->internal_lod_mode_.last.force = 2;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  case LodMode::Never:
    pimpl_->internal_lod_mode_.incr.force = 0;
    pimpl_->internal_lod_mode_.incr.level = 0;
    pimpl_->internal_lod_mode_.last.force = 0;
    pimpl_->internal_lod_mode_.last.level = 0;
    break;

  case LodMode::Rebuild:
    pimpl_->internal_lod_mode_.incr.force = 0;
    pimpl_->internal_lod_mode_.incr.level = 0;
    pimpl_->internal_lod_mode_.last.force = 3;
    pimpl_->internal_lod_mode_.last.level = -1;
    break;

  default:
    throw Errors::ZgyUserError("Invalid value for lodmode");
  }
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodIncrForce(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  switch (n)
  {
  case 0:
  case 1:
    pimpl_->internal_lod_mode_.incr.force = n;
    break;
  default:
    throw Errors::ZgyUserError("Invalid value for lodIncrForce");
  }
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodIncrLevel(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  pimpl_->internal_lod_mode_.incr.level = n;
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodLastForce(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  switch (n)
  {
  case 0:
  case 2:
  case 3:
    pimpl_->internal_lod_mode_.last.force = n;
    break;
  case 1: // Might be supported in the future.
  default:
    throw Errors::ZgyUserError("Invalid value for lodLastForce");
  }
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodLastLevel(int n)
{
  if (pimpl_->has_lod_environ_)
    return *this;
  if (pimpl_->has_lod_lodmode_)
    throw Errors::ZgyUserError
      ("Please don't mix lodmode and the lower level lod settings.");
  pimpl_->has_lod_lowlevel_ = true;
  pimpl_->internal_lod_mode_.last.level = n;
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::historange(float lo, float hi)
{
  pimpl_->historange_ = std::array<float,2>{lo,hi};
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::decimation(const std::vector<DecimationType> decimation)
{
  pimpl_->decimation_ = decimation;

  // The call to EnumMapper::mapDecimationTypeToLodAlgorithm()
  // is currently done from ZgyWriter::enableIncrementalLOD().
  // TODO: Now that EnumMapper has been refactored into its own file,
  // it should be possible to store the internal type in both
  // the V2 and the V3 pimpl. Until then, pimpl_->lod_algorithm_
  // should really have been removed so it won't get used by
  // mistake.
  //std::vector<InternalZGY::LodAlgorithm> algorithms =
  //  OpenZGY::Impl::EnumMapper::mapDecimationTypeToLodAlgorithm(decimation);
  //pimpl_->lod_algorithm_ = algorithms;
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::guidsfrom(const std::shared_ptr<OpenZGY::IZgyReader>& reader)
{
  pimpl_->verid_     = reader->verid();
  pimpl_->has_verid_ = true;

  // dataid() and previd() are missing from the IZgyReader API.
  // Adding them would technically be a breaking change.
  // But I want to maintain the information in the file.
  //
  // It is tempting to access the impl layer, like ZgyMeta::verid() does.
  // But IZgyReader cannot access the ZgyReader implementation. And there
  // might even be a chain of IZgyReader instances. None of which will
  // propagate dataid / previd or the raw file header. Not to mention
  // that including impl/meta.h in this file is a rather gross breach
  // of encapsulation. Would also need impl/guid.h, but that one is safe.
#if 0
  pimpl_->dataid_ = InternalZGY::GUID(reader->ih().dataid()).toString();
  pimpl_->previd_ = InternalZGY::GUID(reader->ih().previd()).toString();
  pimpl_->has_dataid_ = true;
  pimpl_->has_previd_ = true;
#endif
  return *this;
}

/**
 * Output the contents of an argunent package
 * for debugging and/or unit tests.
 *
 * This uses impl() to get hold of the actual data members. Normally, the user visible argumeny package is write only.
 */
void
ZgyWriterArgsV3::dump(std::ostream& out) const
{
  const Impl& pi(*pimpl_);

  // The V3 class sometimes uses the internal enums, not the public
  // ones. To simplify regression testing I will have dump() still
  // show the public types.
  using OpenZGY::Impl::EnumMapper;
  const OpenZGY::SampleDataType dt =
    EnumMapper::mapRawDataTypeToSampleDataType(pi.datatype_);
  const OpenZGY::UnitDimension zdim =
    EnumMapper::mapRawVerticalDimensionToUnitDimension(pi.zunitdim_);
  const OpenZGY::UnitDimension hdim =
    EnumMapper::mapRawHorizontalDimensionToUnitDimension(pi.hunitdim_);

  // DecimationType is stored using the public enum. Bad code smell.
  // There is also an unused pimpl_->algorithms_ that could have held
  // the converted type. See ZgyWriterArgsV3::decimation().
  // The call to EnumMapper::mapDecimationTypeToLodAlgorithm()
  // is currently done from ZgyWriter::enableIncrementalLOD().

  // Code from ZgyWriterArgs::dump() in api.h
  // First line should have been ZgyWriterArgsV3, but keep
  // the old name so the unit tests won't need updating.
  // Similarly, there is a line "ZgyWriterArgsV2\n"
  // below that logically should just have been removed.

  out << "ZgyWriterArgs\n"
      << "  filename:    \"" << pi.filename_ << "\"\n"
      << "  iocontext:   " << (pi.iocontext_ ? "*" : "(null)") << "\n"
      << "  compressor:  " << (pi.compressor_ ? "*" : "(null)") << "\n"
      << "  lodcompress: " << (pi.lodcompressor_ ? "*" : "(null)") << "\n"
      << "  " << (pi.have_size_?"*":"") << "size:        "
      << "(" << pi.size_[0]
      << "," << pi.size_[1]
      << "," << pi.size_[2] << ")\n"
      << "  " << (pi.have_bricksize_?"*":"") << "bricksize:   "
      << "(" << pi.bricksize_[0]
      << "," << pi.bricksize_[1]
      << "," << pi.bricksize_[2] << ")\n"
      << "  " << (pi.have_datatype_?"*":"")
      << "datatype:    " << int(dt) << "\n"
      << "  " << (pi.have_datarange_?"*":"")
      << "datarange:   " << pi.datarange_[0]
      << " to " << pi.datarange_[1] << "\n"
      << "  " << (pi.have_zunit_?"*":"")
      << "zunit:       " << int(zdim)
      << " \"" << pi.zunitname_
      << "\" " << pi.zunitfactor_ << "\n"
      << "  " << (pi.have_hunit_?"*":"")
      << "hunit:       " << int(hdim)
      << " \"" << pi.hunitname_
      << "\" " << pi.hunitfactor_ << "\n"
      << "  " << (pi.have_ilstart_||pi.have_ilinc_?"*":"")
      << "ilstart/inc: " << pi.annotstart_[0] << " / " << pi.annotinc_[0] << "\n"
      << "  " << (pi.have_xlstart_||pi.have_xlinc_?"*":"")
      << "xlstart/inc: " << pi.annotstart_[1] << " / " << pi.annotinc_[1] << "\n"
      << "  " << (pi.have_zstart_||pi.have_zinc_?"*":"")
      << "zstart/inc:  "  << pi.zstart_ << " / " << pi.zinc_ << "\n"
      << "  " << (pi.have_corners_?"*":"") << "corner0:     "
      << pi.corners_[0][0] << ", " << pi.corners_[0][1] << "\n"
      << "  " << (pi.have_corners_?"*":"") << "corner1:     "
      << pi.corners_[1][0] << ", " << pi.corners_[1][1] << "\n"
      << "  " << (pi.have_corners_?"*":"") << "corner2:     "
      << pi.corners_[2][0] << ", " << pi.corners_[2][1] << "\n"
      << "  " << (pi.have_corners_?"*":"") << "corner3:     "
      << pi.corners_[3][0] << ", " << pi.corners_[3][1] << "\n";

  // Code from ZgyWriterArgsV2::dump()
  // Eventually remove the next line.
  out << "ZgyWriterArgsV2\n";
  out << "  lowres:     "
      << " incr (" << pi.internal_lod_mode_.incr.level
      << ", " << pi.internal_lod_mode_.incr.force << ")"
      << " last (" << pi.internal_lod_mode_.last.level
      << ", " << pi.internal_lod_mode_.last.force << ")\n"
      << "  historange:  " << pi.historange_[0] << " to"
      << " " << pi.historange_[1] << "\n";
  out << "  decimation: ";
  for (const auto& it : pi.decimation_)
    out << " " << (int)it;
  out << "\n";
  out << "  algorithm:  ";
  for (const auto& it : pi.lod_algorithm_)
    out << " " << (int)it;
  out << "\n";
  if ((!pi.dataid_.empty()) ||
      (!pi.verid_.empty()) ||
      (!pi.previd_.empty())) {
    out << "  dataid:      " << pi.dataid_ << "\n"
        << "  verid:       " << pi.verid_  << "\n"
        << "  previd:      " << pi.previd_ << "\n";
  }
}

/////////////////////////////////////////////////////////////////////////////
///   NAMED PARAMETERS INTRODUCED IN the original ZgyWriterArgs           ///
///   First, the ones that did not inline the code (see api.cpp)          ///
/////////////////////////////////////////////////////////////////////////////

ZgyWriterArgsV3&
ZgyWriterArgsV3::metafrom(const std::shared_ptr<OpenZGY::IZgyReader>& reader)
{
  OpenZGY::IZgyReader *r = reader.get();
  size(r->size()[0], r->size()[1], r->size()[2]);
  bricksize(r->bricksize()[0], r->bricksize()[1], r->bricksize()[2]);
  datatype(r->datatype());
  datarange(r->datarange()[0], r->datarange()[1]);
  zunit(r->zunitdim(), r->zunitname(), r->zunitfactor());
  hunit(r->hunitdim(), r->hunitname(), r->hunitfactor());
  ilstart(r->annotstart()[0]);
  ilinc(r->annotinc()[0]);
  xlstart(r->annotstart()[1]);
  xlinc(r->annotinc()[1]);
  zstart(r->zstart());
  zinc(r->zinc());
  corners(r->corners());
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::merge(const ZgyWriterArgsV3& other)
{
  this->pimpl_->merge(*other.pimpl_);
  return *this;
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::compressor(const std::string& name, const std::vector<std::string>& args)
{
  return compressor(InternalZGY::CompressFactoryImpl::getCompressor(name,args));
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::zfp_compressor(float snr)
{
  return compressor("ZFP", std::vector<std::string>{"snr", std::to_string(snr)});
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::lodcompressor(const std::string& name, const std::vector<std::string>& args)
{
  return lodcompressor(InternalZGY::CompressFactoryImpl::getCompressor(name,args));
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::zfp_lodcompressor(float snr)
{
  return lodcompressor("ZFP", std::vector<std::string>{"snr", std::to_string(snr)});
}

ZgyWriterArgsV3&
ZgyWriterArgsV3::iocontext(const IOContext *value)
{
  pimpl_->iocontext_ = value ? value->clone() : nullptr;
  return *this;
}

/////////////////////////////////////////////////////////////////////////////
///   NAMED PARAMETERS INTRODUCED IN the original ZgyWriterArgs           ///
///   Here are the ones that used to inline the body in api.h             ///
/////////////////////////////////////////////////////////////////////////////

  ZgyWriterArgsV3& ZgyWriterArgsV3::filename(const std::string& value)       { pimpl_->filename_ = value;      return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::compressor(const compressor_t& value)    { pimpl_->compressor_ = value;    return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::lodcompressor(const compressor_t& value) { pimpl_->lodcompressor_ = value; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::size(std::int64_t ni, std::int64_t nj, std::int64_t nk) {
    pimpl_->size_[0] = ni;
    pimpl_->size_[1] = nj;
    pimpl_->size_[2] = nk;
    pimpl_->have_size_ = true;
    return *this;
  }
  ZgyWriterArgsV3& ZgyWriterArgsV3::bricksize(std::int64_t ni, std::int64_t nj, std::int64_t nk) {
    pimpl_->bricksize_[0] = ni;
    pimpl_->bricksize_[1] = nj;
    pimpl_->bricksize_[2] = nk;
    pimpl_->have_bricksize_ = true;
    return *this;
  }
  ZgyWriterArgsV3& ZgyWriterArgsV3::datatype(SampleDataType value) {
    pimpl_->datatype_ = OpenZGY::Impl::EnumMapper::mapSampleDataTypeToRawDataType(value);
    pimpl_->have_datatype_ = true;
    return *this;
  }
  ZgyWriterArgsV3& ZgyWriterArgsV3::datarange(float lo, float hi) {
    pimpl_->datarange_[0] = lo;
    pimpl_->datarange_[1] = hi;
    pimpl_->have_datarange_ = true;
    return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::zunit(UnitDimension dimension, const std::string& name, double factor) {
    pimpl_->zunitdim_    = OpenZGY::Impl::EnumMapper::mapUnitDimensionToRawVerticalDimension(dimension);
    pimpl_->zunitname_   = name;
    pimpl_->zunitfactor_ = factor;
    pimpl_->have_zunit_  = true;
    return *this;
  }
  ZgyWriterArgsV3& ZgyWriterArgsV3::hunit(UnitDimension dimension, const std::string& name, double factor) {
    pimpl_->hunitdim_    = OpenZGY::Impl::EnumMapper::mapUnitDimensionToRawHorizontalDimension(dimension);
    pimpl_->hunitname_   = name;
    pimpl_->hunitfactor_ = factor;
    pimpl_->have_hunit_  = true;
    return *this;
  }
  ZgyWriterArgsV3& ZgyWriterArgsV3::ilstart(float value) { pimpl_->annotstart_[0] = value; pimpl_->have_ilstart_ = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::ilinc(float value)   { pimpl_->annotinc_[0]   = value; pimpl_->have_ilinc_   = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::xlstart(float value) { pimpl_->annotstart_[1] = value; pimpl_->have_xlstart_ = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::xlinc(float value)   { pimpl_->annotinc_[1]   = value; pimpl_->have_xlinc_   = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::zstart(float value)  { pimpl_->zstart_        = value; pimpl_->have_zstart_  = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::zinc(float value)    { pimpl_->zinc_          = value; pimpl_->have_zinc_    = true; return *this; }
  ZgyWriterArgsV3& ZgyWriterArgsV3::corners(const corners_t& value) { pimpl_->corners_ = value; pimpl_->have_corners_ = true; return *this; }

} // namespace
