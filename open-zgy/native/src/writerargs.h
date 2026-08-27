// Copyright 2017-2023, SLB
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

#include "api.h" // @@@@@@ TEMPORARY, until class ZgyWriterArgs is moved here.
#include "declspec.h"
#include "exception.h"

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <map>
#include <functional>
#include <memory>

#define OPENZGY_WRITERARGS_MAJOR_VERSION 2 // Increment when a new struct is added
#define OPENZGY_WRITERARGS_MINOR_VERSION 2 // Increment when new members are added
// 2.2: Added guidsfrom()

namespace InternalZGY
{
  class ZgyWriterArgsV2Impl;
  class ZgyWriterArgsV3Impl;
}

namespace OpenZGY {
#if 0
}
#endif

/**
 * \brief When to compute low resolution computation. Set on file open.
 *
 * \details
 *
 * This is a user visible enum used for the new single-pass
 * implementation. It replaces FinalizeAction in the two-pass code.
 *
 * The enum is only used by ZgyWriterArgs and derived types.
 * The corresponding information in the impl layer is stored
 * as 4 integers.
 *
 * TODO-WIP-BrickedAPI: Consider adding "Late" which will do all the
 * generating in finalize(), but won't process bricks never written.
 * Not really useful for single session writes, because Rebuild is
 * very efficiant for never written bricks anyway. Might be useful for
 * multi session write because then it is more efficient than Rebuild
 * while still not interleaving lodN bricks. Note that multi session
 * compressed files only Never and (a single) Rebuild is allowed.
 *
 * TODO-WIP-BrickedAPI: Check when opening a file for update:
 *   if compressed and lowres does or will exist, this is an error.
 *   if lowres does not exist, only Never and Rebuild are allowed.
 *   if lowres does or will exist, only forbid Never.
 *
 * Thread safety: Enums do not have race conditions.
 */
enum class OPENZGY_API LodMode {
  Default = 4000, ///< Default and recommended setting, currently "Early".
  Early   = 4001, ///< Generates as soon as possible.
  Early1  = 4002, ///< Generates as soon as possible, only lod1.
  Never   = 4003, ///< Don't generate. Might leave stale lowres behind.
  Rebuild = 4004, ///< Do a full rebuild when the file is finalized.
};

/////////////////////////////////////////////////////////////////////////////
/// LEGACY TYPES DO NOT USE -- See V3 version below.                      ///
/////////////////////////////////////////////////////////////////////////////

/**
 * \brief Additions to the argument package for creating a ZGY file.
 * \details
 * As with the base class, use the named parameter idiom.
 * The new members are stored in a pimpl and the accessors are
 * not inlined. This simplifies adding new settings in the future.
 * Too bad I didn't do it that way originally.
 * See SeismicStoreIOContextV2 which had the same problem.
 *
 * Note this minor inconvenience: If chaining together several
 * setters, which is one of the main reasons for the named parameter
 * idiom, the settings from V2 must be listed before those from V1.
 */
class OPENZGY_API ZgyWriterArgsV2 : public ZgyWriterArgs
{
public:
  typedef InternalZGY::ZgyWriterArgsV2Impl Impl;
  struct FromOldWriterArgs{};
  friend class ZgyWriterArgsV3;

  ZgyWriterArgsV2();
  ~ZgyWriterArgsV2(); // TODO-WIP-BrickedAPI: Virtual or not?

  /**
   * A copy constructor from the base type, filling in the V2 part
   * with default values. The FromOldWriterArgs is used to be more
   * explicit about what is going on. And to avoid invoking the
   * constructor by mistake when a V2 instance is to be copied.
   * Regular copy and assign is also needed because of the pimpl
   * which will require a deep copy. Or if not needed, delete them.
   */
  ZgyWriterArgsV2(const ZgyWriterArgs&, FromOldWriterArgs);
  ZgyWriterArgsV2(const ZgyWriterArgsV2&);
  ZgyWriterArgsV2(const ZgyWriterArgsV2&&) = delete;
  const ZgyWriterArgsV2& operator=(const ZgyWriterArgsV2&);
  const ZgyWriterArgsV2& operator=(const ZgyWriterArgsV2&&) = delete;

public: // or rely on friends?
  const Impl& impl() const; // TODO-WIP-BrickedAPI: Virtual or not?

private:
  std::shared_ptr<Impl> pimpl_;

private:
  /**
   * Set parameters from environment variables, if any. Currently
   * these have precedence after hard coded settings in the Impl
   * constructor and before any explicit settings in ZgyWriterArgsV2.
   *
   * Not polymorphic. Actually, the entire class isn't.
   * Called from the constructor.
   */
  void setFromEnv();

public:

  /**
   * \brief Declare when to compute low resolution bricks.
   * \param mode
   *        - Default:    Let the library choose. Currently "Early".
   *        - Early:      Generates as soon as possible.
   *        - Early1:     Generates as soon as possible, only lod1.
   *        - Never:      Don't generate. Might leave stale lowres behind.
   *        - Rebuild:    Do a full rebuild when the file is finalized.
   * \details
   *
   * If a file is going to be reopened and appended to at a later
   * time, it needs to be opened with "Never" on each open except the
   * last one that should specify "Force". Otherwise the "Default" mode,
   * which is the default, might work best both performance wise and
   * to help the application to make a correct progress bar.
   *
   * This is a hard requirement if the file is compressed. Currently
   * it applies to uncompressed files as well. Lifting the restriction
   * for uncompressed files might be too much work to be feasible.
   *
   * Calling lodmode() is a convenience for setting the four lower
   * level flags lod{Incr|Last}{Force|Level}. The rules for what
   * those should be set to are rather obscure.
   *
   * Most lodmode() will set both "Levels" to -1 meaning that there is no
   * additional limitation on how many lowres levels can be generated now.
   *
   * lodmode(Early) uses Force = 1 (only generate bricks where all
   * input bricks are ready) while writing, and 2 (generate for all
   * bricks that are still dirty) when finalizing. Early1 is similar
   * but defers lod2 and above to the finalize step. At that point the
   * histogram will be complete. Which could give better results for
   * the weighted average algorithm used in lod2 and above.
   *
   * lodmode(Never/Force) uses Force = 0 (i.e. never) during writes
   * and 0 or 3 respectively during finalize.
   *
   * See enableIncrementalLOD() where the settings are used.
   */
  ZgyWriterArgsV2& lodmode(LodMode);

  /**
   * Interleaved with writing, decide whether to compute low resolution data.
   *  * 0: Not.
   *  * 1: Only for lod bricks where all input is ready.
   *  * 2: Forbidden: All bricks that are still dirty.
   *  * 3: Forbidden: All bricks. Lowres gets build ridiculously many times.
   */
  ZgyWriterArgsV2& lodIncrForce(int);

  /**
   * Interleaved with writing, decide how many low resolution levels are
   * allowed to be generated. This limitation is in addition to the Force
   * setting. Levels up to and including the supplied level can be
   * generated. "-1" means all (the default) and "0" means none.
   * Setting "1" might concievably give higher quality lowres data
   * at the cost of more work to be done during finalize.
   */
  ZgyWriterArgsV2& lodIncrLevel(int);

  /**
   * When finalizing, decide whether to compute low resolution data.
   *  * 0: Not.
   *  * 1: Only for lod bricks where all input is ready.
   *       Forbidden until (if at all) the code allows to re-open
   *       a file and continue writing incremental lods.
   *  * 2: All bricks that are still dirty.
   *  * 3: All bricks. Do not trust dirty state.
   */
  ZgyWriterArgsV2& lodLastForce(int);

  /**
   * When finalizing, decide how many low resolution levels are
   * allowed to be generated. This limitation is in addition to the
   * Force setting. Levels up to and including the supplied level can
   * be generated. "-1" means all (the default) and "0" means none.
   * Other settings don't make much sense during finalize except for
   * some obscure performance measurements.
   */
  ZgyWriterArgsV2& lodLastLevel(int);

  /**
   * Force a specific histogram range instead of letting the library
   * figure it out by itself. In the latter case the range typically
   * ends up a little too wide. I.e. with the first few and last few
   * buckets being empty.
   *
   * Currently this setting is ignored for integral types.
   * It is used only when a float file is initially created, or when
   * the entire survey is set to a constant value. The latter might
   * happen both on create and during an update.
   */
  ZgyWriterArgsV2& historange(float lo, float hi);

  /**
   * Choose the decimation algorithm to be used. The first entry
   * in the vector is used for LOD 1, the second for LOD 2, etc.
   * Levels higher than the vector length use the last entry.
   */
  ZgyWriterArgsV2& decimation(const std::vector<DecimationType>);

  /**
   * Copy guids from an existing file.
   *
   * This should ONLY be done if the target file has, or will have,
   * the exact same contents.
   *
   * See doc/keepguids.md for a very long documentation.
   */
  ZgyWriterArgsV2& guidsfrom(const std::shared_ptr<OpenZGY::IZgyReader>&);

  /**
   * Output the argument package as text, for debugging only.
   */
  void dump(std::ostream& out) const;

  // Forward to the old arg package, to get the correct (V2) struct
  // returned. So the named parameter idiom still works.
  // CAVEAT: Missing a function here will introduce serious bugs.
  // The wrong open() or reopen() overload will get called, causing
  // all the settings in ZgyWriterArgsV2 to reset to default.
  // Tip: Have the last setting in the change be a V2 parameter
  // if possible. If any of the earlier functions were missing
  // then this will cause a compiler error.

  ZgyWriterArgsV2& filename(const std::string& value);
  ZgyWriterArgsV2& iocontext(const IOContext *value);
  ZgyWriterArgsV2& compressor(const compressor_t& value);
  ZgyWriterArgsV2& compressor(const std::string& name, const std::vector<std::string>& args);
  ZgyWriterArgsV2& zfp_compressor(float snr);
  ZgyWriterArgsV2& lodcompressor(const compressor_t& value);
  ZgyWriterArgsV2& lodcompressor(const std::string& name, const std::vector<std::string>& args);
  ZgyWriterArgsV2& zfp_lodcompressor(float snr);
  ZgyWriterArgsV2& size(std::int64_t ni, std::int64_t nj, std::int64_t nk);
  ZgyWriterArgsV2& bricksize(std::int64_t ni, std::int64_t nj, std::int64_t nk);
  ZgyWriterArgsV2& datatype(SampleDataType value);
  ZgyWriterArgsV2& datarange(float lo, float hi);
  ZgyWriterArgsV2& zunit(UnitDimension dimension, const std::string& name, double factor);
  ZgyWriterArgsV2& hunit(UnitDimension dimension, const std::string& name, double factor);
  ZgyWriterArgsV2& ilstart(float value);
  ZgyWriterArgsV2& ilinc(float value);
  ZgyWriterArgsV2& xlstart(float value);
  ZgyWriterArgsV2& xlinc(float value);
  ZgyWriterArgsV2& zstart(float value);
  ZgyWriterArgsV2& zinc(float value);
  ZgyWriterArgsV2& corners(const corners_t& value);
  ZgyWriterArgsV2& metafrom(const std::shared_ptr<OpenZGY::IZgyReader>&);
  ZgyWriterArgsV2& merge(const ZgyWriterArgsV2&);
};

/////////////////////////////////////////////////////////////////////////////
/// CURRENT TYPES                                                         ///
/////////////////////////////////////////////////////////////////////////////

/**
 * \brief Additions to the argument package for creating a ZGY file.
 * \details
 * As with the base class, use the named parameter idiom.
 * The new members are stored in a pimpl and the accessors are
 * not inlined. This simplifies adding new settings in the future.
 */
class OPENZGY_API ZgyWriterArgsV3
{
public:
  typedef InternalZGY::ZgyWriterArgsV3Impl Impl;
  struct FromOldWriterArgs{};
  struct FromOldWriterArgsV2{};
  typedef double float64_t;
  typedef std::array<std::array<float64_t,2>,4> corners_t;
  typedef std::pair<std::shared_ptr<const void>, std::int64_t> rawdata_t;
  typedef std::function<rawdata_t(const rawdata_t&,const std::array<int64_t,3>&)> compressor_t;

  ZgyWriterArgsV3();
  ~ZgyWriterArgsV3();

  /**
   * Regular copy constructor and assign operator are needed because
   * of the pimpl which will require a deep copy.
   *
   * Move constructor and move assign are added because they can steal
   * the pimpl instead of cloning it. Currently the performance
   * difference is hardly measurable. But if it turns out that the the
   * move becomes useful later, it might be tricky to add without
   * breaking the ABI.
   */
  ZgyWriterArgsV3(const ZgyWriterArgsV3&);
  ZgyWriterArgsV3(ZgyWriterArgsV3&&);
  const ZgyWriterArgsV3& operator=(const ZgyWriterArgsV3&);
  const ZgyWriterArgsV3& operator=(ZgyWriterArgsV3&&);

  /**
   * Copy constructors from the older, deprecated writer args.
   * Members not present in the older argument packages will get
   * default values.
   *
   * A tag is used to be more explicit about which of the constructors
   * we want. E.g. the ZgyWriterArgs& version and the ZgyWriterArgsV2&
   * version can both handle a ZgyWriterArgsV2 input. C++ will
   * probably do the right thing and choose the second one because it
   * is a better match. But it won't hurt to be explicit.
   */
  ZgyWriterArgsV3(const ZgyWriterArgs&, FromOldWriterArgs);
  ZgyWriterArgsV3(const ZgyWriterArgsV2&, FromOldWriterArgsV2);

public: // or rely on friends?
  const Impl& impl() const;

private:
  std::shared_ptr<Impl> pimpl_;

private:
  /**
   * Set parameters from environment variables, if any. Currently
   * these have precedence after hard coded settings in the Impl
   * constructor and before any explicit settings in ZgyWriterArgsV3.
   *
   * Not polymorphic. Actually, the entire class isn't.
   * Called from the constructor.
   */
  void setFromEnv();

  ZgyWriterArgsV3& copyFromV1(const ZgyWriterArgs& old);
  ZgyWriterArgsV3& copyFromV2(const ZgyWriterArgsV2& old);

public:

  /**
   * \brief Declare when to compute low resolution bricks.
   * \param mode
   *        - Default:    Let the library choose. Currently "Early".
   *        - Early:      Generates as soon as possible.
   *        - Early1:     Generates as soon as possible, only lod1.
   *        - Never:      Don't generate. Might leave stale lowres behind.
   *        - Rebuild:    Do a full rebuild when the file is finalized.
   * \details
   *
   * If a file is going to be reopened and appended to at a later
   * time, it needs to be opened with "Never" on each open except the
   * last one that should specify "Force". Otherwise the "Default" mode,
   * which is the default, might work best both performance wise and
   * to help the application to make a correct progress bar.
   *
   * This is a hard requirement if the file is compressed. Currently
   * it applies to uncompressed files as well. Lifting the restriction
   * for uncompressed files might be too much work to be feasible.
   *
   * Calling lodmode() is a convenience for setting the four lower
   * level flags lod{Incr|Last}{Force|Level}. The rules for what
   * those should be set to are rather obscure.
   *
   * Most lodmode() will set both "Levels" to -1 meaning that there is no
   * additional limitation on how many lowres levels can be generated now.
   *
   * lodmode(Early) uses Force = 1 (only generate bricks where all
   * input bricks are ready) while writing, and 2 (generate for all
   * bricks that are still dirty) when finalizing. Early1 is similar
   * but defers lod2 and above to the finalize step. At that point the
   * histogram will be complete. Which could give better results for
   * the weighted average algorithm used in lod2 and above.
   *
   * lodmode(Never/Force) uses Force = 0 (i.e. never) during writes
   * and 0 or 3 respectively during finalize.
   *
   * See enableIncrementalLOD() where the settings are used.
   */
  ZgyWriterArgsV3& lodmode(LodMode);

  /**
   * Interleaved with writing, decide whether to compute low resolution data.
   *  * 0: Not.
   *  * 1: Only for lod bricks where all input is ready.
   *  * 2: Forbidden: All bricks that are still dirty.
   *  * 3: Forbidden: All bricks. Lowres gets build ridiculously many times.
   */
  ZgyWriterArgsV3& lodIncrForce(int);

  /**
   * Interleaved with writing, decide how many low resolution levels are
   * allowed to be generated. This limitation is in addition to the Force
   * setting. Levels up to and including the supplied level can be
   * generated. "-1" means all (the default) and "0" means none.
   * Setting "1" might concievably give higher quality lowres data
   * at the cost of more work to be done during finalize.
   */
  ZgyWriterArgsV3& lodIncrLevel(int);

  /**
   * When finalizing, decide whether to compute low resolution data.
   *  * 0: Not.
   *  * 1: Only for lod bricks where all input is ready.
   *       Forbidden until (if at all) the code allows to re-open
   *       a file and continue writing incremental lods.
   *  * 2: All bricks that are still dirty.
   *  * 3: All bricks. Do not trust dirty state.
   */
  ZgyWriterArgsV3& lodLastForce(int);

  /**
   * When finalizing, decide how many low resolution levels are
   * allowed to be generated. This limitation is in addition to the
   * Force setting. Levels up to and including the supplied level can
   * be generated. "-1" means all (the default) and "0" means none.
   * Other settings don't make much sense during finalize except for
   * some obscure performance measurements.
   */
  ZgyWriterArgsV3& lodLastLevel(int);

  /**
   * Force a specific histogram range instead of letting the library
   * figure it out by itself. In the latter case the range typically
   * ends up a little too wide. I.e. with the first few and last few
   * buckets being empty.
   *
   * Currently this setting is ignored for integral types.
   * It is used only when a float file is initially created, or when
   * the entire survey is set to a constant value. The latter might
   * happen both on create and during an update.
   */
  ZgyWriterArgsV3& historange(float lo, float hi);

  /**
   * Choose the decimation algorithm to be used. The first entry
   * in the vector is used for LOD 1, the second for LOD 2, etc.
   * Levels higher than the vector length use the last entry.
   */
  ZgyWriterArgsV3& decimation(const std::vector<DecimationType>);

  /**
   * Copy guids from an existing file.
   *
   * This should ONLY be done if the target file has, or will have,
   * the exact same contents.
   *
   * See doc/keepguids.md for a very long documentation.
   */
  ZgyWriterArgsV3& guidsfrom(const std::shared_ptr<OpenZGY::IZgyReader>&);

  /**
   * Output the argument package as text, for debugging only.
   */
  void dump(std::ostream& out) const;

  /**
   * \brief Set file to open.
   * \details If starting with sd:// this opens a file on seismic store.
   */
  ZgyWriterArgsV3& filename(const std::string& value);

  /**
   * \brief Set credentials and other configuration.
   * \details This depends on the back-end.
   * For local files you normally don't need to pass anything here.
   */
  ZgyWriterArgsV3& iocontext(const IOContext *value);

  /**
   * \brief Set functor for compressing full resolution data.
   */
  ZgyWriterArgsV3& compressor(const compressor_t& value);

  /**
   * \brief Set functor for compressing full resolution data.
   * \details This overload uses a factory to look up the compressor.
   */
  ZgyWriterArgsV3& compressor(const std::string& name, const std::vector<std::string>& args);

  /**
   * \brief Set functor for compressing full resolution data.
   * \details This overload uses a factory to look up the ZGY compressor.
   * It is just a convenience that is shorter to type.
   */
  ZgyWriterArgsV3& zfp_compressor(float snr);

  /**
   * \brief Set functor for compressing low resolution data.
   */
  ZgyWriterArgsV3& lodcompressor(const compressor_t& value);

  /**
   * \brief Set functor for compressing low resolution data.
   * \details This overload uses a factory to look up the compressor.
   */
  ZgyWriterArgsV3& lodcompressor(const std::string& name, const std::vector<std::string>& args);

  /**
   * \brief Set functor for compressing low resolution data.
   * \details This overload uses a factory to look up the ZGY compressor.
   * It is just a convenience that is shorter to type.
   */
  ZgyWriterArgsV3& zfp_lodcompressor(float snr);

  /**
   * \brief Set size of the survey.
   * \param ni number of inlines (slowest varying index).
   * \param nj number of crosslines.
   * \param nk number of samples per trace (fastest).
   */
  ZgyWriterArgsV3& size(std::int64_t ni, std::int64_t nj, std::int64_t nk);

  /**
   * \brief Set size of one brick.
   * \details Almost always (64,64,64). Change at your own peril.
   */
  ZgyWriterArgsV3& bricksize(std::int64_t ni, std::int64_t nj, std::int64_t nk);

  /**
   * \brief Set type of samples in each brick.
   */
  ZgyWriterArgsV3& datatype(SampleDataType value);

  /**
   * \brief Set scaling factors.
   * \details For integral storage this specifies the two floating
   * point numbers that correspond to the lowest and highest
   * representable integer value. So for int8 files, -128 will be
   * converted to "lo" and +127 will be converted to "hi".
   *
   * ZGY files require that min<max. This is not enforced here,
   * but is checked when the file is actually created.
   */
  ZgyWriterArgsV3& datarange(float lo, float hi);

  /**
   * \brief Set vertical unit.
   * \param dimension time or depth (a.k.a. length).
   * \param name for annotation only.
   * \param factor multiply by this factor to convert from storage units to SI units.
   */
  ZgyWriterArgsV3& zunit(UnitDimension dimension, const std::string& name, double factor);

  /**
   * \brief Set horizontal unit.
   * \param dimension cartesian length or (unsupported) polat coordinates,
   * \param name for annotation only.
   * \param factor multiply by this factor to convert from storage units to SI units.
   */
  ZgyWriterArgsV3& hunit(UnitDimension dimension, const std::string& name, double factor);

  /// \brief Set first (ordinal 0) inline number.
  /// \details For maximum portability the inline and crossline start
  /// and increment should be integral numbers. Some applications
  /// might choose to convert them to int.
  ZgyWriterArgsV3& ilstart(float value);

  /// \brief Set inline number increment between two adjacent ordinal values.
  /// \copydetails ilstart
  ZgyWriterArgsV3& ilinc(float value);

  /// \brief Set first (ordinal 0) crossline number.
  /// \copydetails ilstart
  ZgyWriterArgsV3& xlstart(float value);

  /// \brief Set crossline number increment between two adjacent ordinal values.
  /// \copydetails ilstart
  ZgyWriterArgsV3& xlinc(float value);

  /// \brief Set first time/depth.
  /// \details Vertical annotation is generally safe to have non-integral.
  ZgyWriterArgsV3& zstart(float value);

  /// \brief Set increment (distance between samples) in vertical direction.
  /// \copydetails zstart
  ZgyWriterArgsV3& zinc(float value);

  /// \brief Set survey corner points in world coordinates.
  /// \details The corners are ordered origin, last inline (i.e. i=last, j=0),
  /// last crossline, diagonal.
  ZgyWriterArgsV3& corners(const corners_t& value);

  /**
   * \brief Copy metadata from existing file.
   *
   * Set most of the metadata from an open file descriptor.
   * Typically used when the file to be created is to be computed from
   * an existing file. Typically this will be called first so the
   * settings don't inadverently shadow something set explicitly.
   */
  ZgyWriterArgsV3& metafrom(const std::shared_ptr<OpenZGY::IZgyReader>&);

  /**
   * \brief Copy metadata from another ZgyWriterArgs.
   *
   * Copy only those settings that have been explicitly changed in the
   * supplied "other" ZgyWriterArgs into *this. If other.metafrom()
   * has been called it is unspecified what gets copied.
   *
   * Compared to the version in ZgyWriterArgsV2, this one is more
   * consistent. Including filename, iocontext, and the two compressors.
   * Neither of those are likely to be useful, but it is safer to do
   * everything.
   */
  ZgyWriterArgsV3& merge(const ZgyWriterArgsV3&);
};

} // namespace
