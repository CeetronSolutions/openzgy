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

/**
 * \file meta.h
 * \brief Meta data read/write.
 *
 * The meta data in an OpenZGY file is broken down into multiple
 * headers. In this file, a set of classes for each header type
 * is responsble for reading and decoding that header. A lot of
 * the code is inlined to make it easier to maintain. Keeping
 * everything that needs to be updated close together.
 *
 * class FooHeaderV1POD\n
 * class FooHeaderV1POD
 * etc.
 *
 * \li POD class corresponding exactly to what is stored on the file.
 *     A new class needs to be created each time the version changes.
 *     The compiler gets to decide the actual layout, but is instructed
 *     to not use any padding. Hopefully this makes the layout platform
 *     independant.
 *
 * class IFooHeaderAccess : public IHeaderAccess
 *
 * \li Pure interface, more or less matching the layout of the most recent
 *     version of the POD. The purpose is to hide differences between versions
 *     when reading. Writing and updating need not be handled the same way.
 *     Those are simpler. Creating a new file only needs to worry about
 *     the latest version, and updating meta data is only allowed on very
 *     few (currently none at all) attributes.
 *
 * class FooHeaderAccess : public IHeaderAccess
 *
 * \li Adds a concrete dump() method for debugging.
 *
 * class FooHeaderV1Access : public FooHeaderAccess\n
 * class FooHeaderV2Access : public FooHeaderAccess\n
 * etc.
 *
 * \li Aggregates (i.e. does not inherit from) the corresponding POD.
 *     Implements the interface that allows version independant access
 *     of the information on file. Where there is a mismatch the code
 *     should compute the result on the fly if possible. The last resort
 *     is to cache data in the virtual compute() method and storing
 *     the result in this instance instead if the pod. Note that this
 *     gets messy if the cached data can become stale.
 *
 * \li TODO-WARNING, the current implementation assumes that the computer
 *     architecture and the C++ compiler allows unaligned access.
 *     According to the C++ standard this gives ***UNDEFINED BEHAVIOR***
 *     even on architectures such as x86_684 that do allow unaligned loads.
 *     On newer x86_64 processors this is actually fairly efficient too.
 *     BUT the compiler is free to use SSE instructions for load/store
 *     and those do NOT allow unaligned pointers.
 *     https://pzemtsov.github.io/2016/11/06/bug-story-alignment-on-x86.html
 *
 * HeaderAccessFactory::createFoo(std::uint32_t version)
 *
 * \li Construct and return the appropriate FooHeaderV?Access instance.
 *
 * This file was initially generated with some help from openzgy.tools.cppmeta.
 * But that was only a starting point. If major changes are needed to the format
 * (will hopefully never happen) then the code generator might still be useful.
 * otherwise this file is expected to be maintained by hand.
 */

#pragma once

#include "../declspec.h"
#include "../exception.h"
#include "file.h"
#include "enum.h"
#include "lookuptable.h"

#include <cstdint>
#include <memory>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <atomic>

namespace InternalZGY {
#if 0
}
#endif

class ZgyWriterArgsV2Impl;
class ZgyWriterArgsV3Impl;

class IInfoHeaderAccess;
class IHistHeaderAccess;

/*
 * Currently I am having the access function for foo[N] return std::array
 * instead of a raw pointer. This is inefficient but simplifies the case
 * where the result is a derived value. I might change my mind about this.
 *
 * The ptr_to_array method allows writing
 *    return ptr_to_array<float,3>(_foo);
 * instead of
 *    return std::array<float,3>{_foo[0], _foo[1], _foo[2]};
 * and also tries to handle the case when _foo is misaligned.
 */

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class IHeaderAccess
{
public:
  virtual ~IHeaderAccess();
  typedef std::vector<std::uint8_t> podbytes_t;
  virtual podbytes_t podbytes() const = 0;
  virtual void dump(std::ostream& out, const std::string& prefix = "") = 0;
  // read(), byteswap(), calculate() not useful to abstract at this level.
};

/////////////////////////////////////////////////////////////////////////////
//    FileHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class IFileHeaderAccess : public IHeaderAccess
{
public:
  virtual void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset) = 0;
  virtual void byteswap() = 0;
public:
  virtual std::array<std::uint8_t,4> magic() const = 0;
  virtual std::uint32_t version() const = 0;
  virtual void set_version(std::uint32_t value) = 0;
};

/////////////////////////////////////////////////////////////////////////////
//    OffsetHeader   ////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class IOffsetHeaderAccess : public IHeaderAccess
{
public:
  virtual void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset) = 0;
  virtual void byteswap() = 0;
  virtual void calculate(const std::shared_ptr<IInfoHeaderAccess>& ih) = 0;
public:
  virtual std::int64_t infoff()       const = 0;
  virtual std::int64_t stroff()       const = 0;
  virtual std::int64_t alphalupoff()  const = 0;
  virtual std::int64_t bricklupoff()  const = 0;
  virtual std::int64_t histoff()      const = 0;
  virtual std::int64_t infsize()      const = 0;
  virtual std::int64_t histsize()     const = 0;
  virtual std::int64_t alphalupsize() const = 0;
  virtual std::int64_t bricklupsize() const = 0;
};

/////////////////////////////////////////////////////////////////////////////
//    InfoHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class IInfoHeaderAccess : public IHeaderAccess
{
public:
  virtual void       read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) = 0;
  virtual void       byteswap() = 0;
  virtual void       calculate_cache() = 0;
  virtual void       calculate_read(const podbytes_t& slbuf, const std::shared_ptr<IHistHeaderAccess>& hh) = 0;
  virtual podbytes_t calculate_write() = 0;
public:
  virtual std::array<std::int64_t,3> bricksize() const = 0;
  virtual std::array<float,2>        safe_codingrange() const = 0;
  virtual std::array<float,2>        raw_codingrange() const = 0;
  virtual std::array<std::uint8_t,16>dataid() const = 0;
  virtual std::array<std::uint8_t,16>verid() const = 0;
  virtual std::array<std::uint8_t,16>previd() const = 0;
  virtual std::string                srcname() const = 0;
  virtual std::string                srcdesc() const = 0;
  virtual RawDataType                srctype() const = 0;
  virtual std::array<std::int64_t,3> size() const = 0;
  virtual std::array<float,3>        orig() const = 0;
  virtual std::array<float,3>        inc()  const = 0;
  virtual std::int64_t               scnt() const = 0;
  virtual double                     ssum() const = 0;
  virtual double                     sssq() const = 0;
  virtual float                      smin() const = 0;
  virtual float                      smax() const = 0;
  virtual std::array<float,4>        gpiline() const = 0;
  virtual std::array<float,4>        gpxline() const = 0;
  virtual std::array<double,4>       gpx() const = 0;
  virtual std::array<double,4>       gpy() const = 0;
  virtual RawDataType                datatype() const = 0;
  virtual std::string                hprjsys() const = 0;
  virtual RawHorizontalDimension     hdim() const = 0;
  virtual double                     hunitfactor() const = 0;
  virtual std::string                hunitname() const = 0;
  virtual RawVerticalDimension       vdim() const = 0;
  virtual double                     vunitfactor() const = 0;
  virtual std::string                vunitname() const = 0;
  virtual std::uint32_t              slbufsize() const = 0;
  virtual const std::array<std::array<double,2>,4>& ocp_index() const = 0;
  virtual const std::array<std::array<double,2>,4>& ocp_annot() const = 0;
  virtual const std::array<std::array<double,2>,4>& ocp_world() const = 0;
  virtual std::int32_t               nlods() const = 0;
  // The next three are derived but cannot possibly change after being set,
  // so returning a reference to cached data shold be safe.
  virtual const std::vector<std::array<std::int64_t,3>>& lodsizes() const = 0;
  virtual const std::vector<std::int64_t>& alphaoffsets() const = 0;
  virtual const std::vector<std::int64_t>& brickoffsets() const = 0;
  // More derived stuff, TODO-Performance may be cached.
  virtual std::int64_t bytesperalpha() const = 0;
  virtual std::int64_t bytesperbrick() const = 0;
  virtual std::int64_t bytespersample() const = 0;
  virtual std::array<double,2>       storagetofloat() const = 0;
  virtual double                     storagetofloat_slope() const = 0;
  virtual double                     storagetofloat_intercept() const = 0;
  virtual double                     defaultstorage() const = 0;
  virtual double                     defaultvalue() const = 0;
  // Write support.
  virtual void setstats(std::int64_t scnt, double ssum, double sssq,
                        double smin, double smax) = 0;
};

/////////////////////////////////////////////////////////////////////////////
//    HistHeader   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class IHistHeaderAccess : public IHeaderAccess
{
public:
  virtual void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) = 0;
  virtual void byteswap() = 0;
  virtual void calculate() = 0;
public:
  virtual std::int64_t bincount()    const = 0;
  virtual std::int64_t samplecount() const = 0;
  virtual double       minvalue()    const = 0;
  virtual double       maxvalue()    const = 0;
  virtual const std::int64_t* bins() const = 0;
  // Write support.
  virtual void sethisto(double minvalue, double maxvalue,
                        const std::int64_t* bins, std::int64_t bincount) = 0;
};

/////////////////////////////////////////////////////////////////////////////
//    AlphaLUP, BrickLUP   //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: Interfaces do not have race conditions.
 */
class ILookupTableAccess : public IHeaderAccess
{
public:
  virtual void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size) = 0;
  virtual void byteswap() = 0;
public:
  virtual std::uint64_t lookupLinearIndex(std::int64_t index) const = 0;
  virtual std::vector<std::uint64_t>& lup() = 0;
  virtual std::vector<std::uint64_t>& lupend() = 0;
  virtual const std::vector<std::uint64_t>& lup() const = 0;
  virtual const std::vector<std::uint64_t>& lupend() const = 0;
};

/////////////////////////////////////////////////////////////////////////////
//    HeaderAccessFactory   (for all known header types)   //////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \details Thread safety: No issues. The only way there can be
 * concurrent access is when two threads are each opening a file. And
 * that is ok.
 */
class HeaderAccessFactory
{
  HeaderAccessFactory() = delete;
  HeaderAccessFactory(const HeaderAccessFactory&) = delete;
public:
  static std::shared_ptr<IFileHeaderAccess>   createFileHeader();
  static std::shared_ptr<IOffsetHeaderAccess> createOffsetHeader(std::uint32_t);
  static std::shared_ptr<IInfoHeaderAccess>   createInfoHeader(std::uint32_t);
  static std::shared_ptr<IHistHeaderAccess>   createHistHeader(std::uint32_t);
  static std::shared_ptr<ILookupTableAccess>  createAlphaLookup(std::uint32_t);
  static std::shared_ptr<ILookupTableAccess>  createBrickLookup(std::uint32_t);
};

/////////////////////////////////////////////////////////////////////////////
//    ZgyInternalMeta   (this is the top level)   ///////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Holds references to all the individual headers needed to access ZGY.
 *
 * Thread safety: By design, OpenZGY is expected to allow concurrent reads.
 * Which means that no read operation may change any data member. This
 * includes lazy evaluation. Alternatively the data members that might be
 * updated need to be protected with a lock.
 *
 * The aggregated data members are all private and exposed via const-correct
 * member functions. This means that if the instance is const then callers
 * can only access const versions of aggregated data.
 */
class OPENZGY_TEST_API ZgyInternalMeta
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;

private:
  std::shared_ptr<IFileHeaderAccess>   _fh;
  std::shared_ptr<IOffsetHeaderAccess>   _oh;
  std::shared_ptr<IInfoHeaderAccess>     _ih;
  std::shared_ptr<IHistHeaderAccess>     _hh;
  std::shared_ptr<ILookupTableAccess>   _alup;
  std::shared_ptr<ILookupTableAccess>   _blup;
  LoggerFn                              _loggerfn;
  std::atomic<bool>                     _is_bad;
  bool                                  _cached_hasBrickLOD;
  std::int32_t                          _cached_nlods;

public:
  // Actually internal. Used by OpenZGY::ZgyMeta, ZgyInternalBulk, GenLodC.
  const IFileHeaderAccess&  fh()   const { return *_fh; }
  const IOffsetHeaderAccess& oh() const { return *_oh; }
  const IInfoHeaderAccess&  ih()   const { return *_ih; }
  const IHistHeaderAccess&  hh()   const { return *_hh; }
  const ILookupTableAccess& alup() const { return *_alup; }
  const ILookupTableAccess& blup() const { return *_blup; }
  bool  cached_hasBrickLOD()       const { return _cached_hasBrickLOD; }
  std::int32_t cached_nlods()      const { return _cached_nlods; }
  // Writing a file needs writable lookup tables and
  // writable statistics and histogram. Stats are inside
  // the info header so this grants a bit more than needed.
  // File header must be writable because close_incomplete()
  // needs to change the version number to 4.
  IFileHeaderAccess&  fh()   { return *_fh; }
  IInfoHeaderAccess&  ih()   { return *_ih; }
  IHistHeaderAccess&  hh()   { return *_hh; }
  ILookupTableAccess& alup() { return *_alup; }
  ILookupTableAccess& blup() { return *_blup; }

public:
  ZgyInternalMeta();
  ZgyInternalMeta(const ZgyInternalMeta&) = delete;
  ZgyInternalMeta& operator=(const ZgyInternalMeta&) = delete;
  explicit ZgyInternalMeta(const std::shared_ptr<IFileADT>& file,
                           const LoggerFn& logger = LoggerFn());
  explicit ZgyInternalMeta(const ZgyWriterArgsV3Impl& v3args,
                           bool compress,
                           const LoggerFn& logger = LoggerFn());
  explicit ZgyInternalMeta(const std::shared_ptr<IFileADT>& file,
                           const ZgyWriterArgsV3Impl& v3args,
                           bool compress,
                           const LoggerFn& logger = LoggerFn());
public:
  void dump(std::ostream& out, const std::string& prefix = "");
  bool errorflag() const        { return _is_bad.load(); }
  void set_errorflag(bool flag) { _is_bad.store(flag); }
  LoggerFn set_logger(const LoggerFn& logger) {
    LoggerFn old = _loggerfn;
    _loggerfn = logger;
    return old;
  }
  std::int64_t flushMeta(const std::shared_ptr<IFileADT>& file);
private:
  // Caveat: Adding bloat to this class. Find a better design?
  // Add methods to get derived information from one or more of the headers.
  // Add methods that essentially forward the request to the static class
  // LookupTable. Being static, all context needed by LookupTable must be
  // passed in every call. And most of what it needs comes from our headers.
  bool has_stathist_feature()     const; // Statistics and histogram are good.
  bool has_tiny_feature()         const; // Survey size one brick.
  bool has_lowres_feature()       const; // Has low resolution. True if tiny.
  bool has_compression_feature()  const; // Contains compressed bricks.
  bool has_empty_feature()        const; // There is no allocated data.
public:
  bool can_old_library_read()     const; // Readable by the old library.
  bool can_append_bulk(bool)      const;

private:
  bool _logger(int priority, const std::string& ss = std::string()) const;
  bool _logger(int priority, const std::ios& ss) const;
  void initFromOpenFile(const std::shared_ptr<IFileADT>&);
  void initFromScratch(const ZgyWriterArgsV3Impl&, bool compress);
  void initFromReopen(const ZgyWriterArgsV3Impl&, bool compress);
  static ZgyWriterArgsV3Impl validateCreateNewFile(const ZgyWriterArgsV3Impl&);
  static ZgyWriterArgsV3Impl validateOpenForUpdate(const ZgyWriterArgsV3Impl&);
  class ErrorsWillCorruptFile;
};

} // namespace
