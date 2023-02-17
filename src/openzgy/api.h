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

#include <cstdint>
#include <vector>
#include <array>
#include <ostream>
#include <iostream>
#include <functional>
#include <memory>
#include <ostream>
#include <iostream>
#include <sstream>
#include <string>
#include <mutex>

#include "declspec.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY {
  /**
   * \mainpage
   *
   * The %OpenZGY C++ API allows read/write access to files
   * stored in the ZGY format. The main part of the API is here:
   *
   * \li IZgyReader and its IZgyMeta base class.
   * \li IZgyWriter and its ZgyWriterArgs argument package.
   * \li IZgyUtils for anything not read or write.
   * \li \ref exceptions
   *  that you might want to catch.
   * \if SSTORE
   * \li SeismicStoreIOContext for cloud credentials.
   * \endif
   * \li ProgressWithDots example of progress reporting.
   * \li \ref Example Example application.
   *
   * If you are reading this document from doxygen/native/apidoc.pdf
   * in the source tree then please see doxygen/README.md for an
   * explanation of why the documentation produced by the build might
   * be better.
   *
   * \if IMPL
   * If you are viewing the full Doxygen documentation then this
   * covers both the API and most of the implementation. So if you
   * look at the list of classes and methods this might seem a bit
   * daunting. All you really need to use the API should be in the
   * above list. Excluding trivial structs that will be cross
   * referenced as needed. So you don't need to go looking for them.
   * Of course, if you want to work om %OpenZGY itself then you
   * probably need everything.
   *
   * See also the following related pages:
   *
   * \li \ref physicalformat
   * \li \ref implementation
   * \li \ref lowres
   * \li \ref migration
   *
   * \endif
   *
   * \internal
   *
   * The public version of the documentation will include both the
   * IZgyReader and the ZgyReader class, with identical descriptions.
   * Same for IZgyWriter et cetera. I don't want ZgyReader and friends
   * to show up but the documentation is attached to those classes and
   * the interface just has a copydoc. This is done deliberately to
   * keep the documentation close to the inplementation so it is easy
   * to keep in sync.
   *
   * Doxygen's cond statement doesn't do what I want. That will exclude
   * the docmentation completely. Whay I want is to keep the information
   * for use by \copydoc but ignore it for any other purpose.
   *
   * \endinternal
   *
   * \namespace OpenZGY
   *
   * \brief The entire public API is in this namespace.
   *
   * The main interface is IZgyReader and IZgyWriter.
   *
   * \namespace ::OpenZGY::Errors
   *
   * \brief Exceptions that can be thrown by OpenZGY.
   *
   * \namespace ::OpenZGY::Impl
   *
   * \brief Implementation of the abstract interfaces in OpenZGY
   *
   * \namespace ::OpenZGY::Formatters
   *
   * \brief operator<< for readable output of enums etc.
   *
   * \if IMPL
   * \namespace InternalZGY   *
   * \brief Implementation not visible to clients.
   * \endif
   *
   * \page Example
   * \include example.h
   */
}

namespace Test {
  class ZgyWriterMock;
}

namespace OpenZGY {
#if 0
}
#endif

class IOContext;
class IZgyReader;

namespace Impl {
  class ZgyMeta;
  class ZgyReader;
  class ZgyWriter;
  class EnumMapper;
}

/** \file api.h
 *
 * \brief IZgyReader, IZgyWriter, and other user visible classes.
 *
 * This file contains the public OpenZGY API.
 *
 * The API is modeled roughly after the API exposed by the Python wrapper
 * around the existing C++ ZGY-Public API. This is probably just as good
 * a starting point than anything else. And it makes testing simpler for
 * those tests that compare the old and new behavior.
 *
 * OpenZGY::IZgyMeta
 *
 * \li Has a private reference to a impl.meta.ZgyInternalMeta instance.
 * \li Contains a large number of properties exposing meta data,
 *     most of which will present information from the ZgyInternalMeta
 *     in a way that is simpler to use and doesn't depend on the
 *     file version.
 * \li End users will access methods from this class. Most likely
 *     via the derived ZgyReader and ZgyWriter classes. The end users
 *     might not care that there is a separate base class.
 *
 * OpenZGY::ZgyMetaAndTools extends IZgyMeta
 *
 * \li Add coordinate conversion routines to a ZgyMeta instance.
 *
 * OpenZGY::ZgyReader extends ZgyMetaAndTools with read, readconst, close
 * OpenZGY::ZgyWriter extends ZgyMetaAndTools with write, writeconst, finalize, close
 *
 * \li Has a private reference to a impl.meta.ZgyInternalBulk instance.
 * \li Add bulk read and write functionality, forwarding the requests
 *     to the internal bulk instance.
 * \li These classes with their own and inherited methods and properties
 *     comprise the public OpenZGY API.
 * \li Currently the ZgyWriter does not expose the bulk read() function
 *     but it does allow accessing all the metadata. Allowing read()
 *     might be added later if it appears to be useful.
 *     In practice this just means to let ZgyWriter inherit ZgyReader.
 *
 * Note that enums and exceptions are problematic when it comes to
 * encapsulation. Enums might be:
 *
 * \li    Only visible in the public api. Often end up being mapped
 *        to a corresponding internal enum. No problem, but the mapping
 *        can be a bit tedious to maintain.
 *
 * \li    Only visible internally. Possibly representing integer values
 *        written to file. Which makes them an implementation detail,
 *        which must absolutely not be exposed publically.
 *
 * \li    Defined by the public api but passed unchanged from the api
 *        layer to the implementation layer. So there will be an
 *        include in some impl/xxx.h file to a part of the public API.
 *        This is frowned upon.
 *
 * \li    Defined by the internal api but may need to be used from
 *        applications, i.e. also visible in the public api. There will
 *        be an include of e.g. "impl/ugly.h" from one of the public
 *        header files and/or application code.
 *        This is strongly discouraged except for testing.
 *
 * Exceptions have similar issues. It is possible to catch all exceptions
 * raised in the impl layer and convert those to exceptions owned by the
 * api layer. But re-throwing an exception makes debuggig harder.
 */

/**
 * \brief Sample data type used in the public API.
 *
 * Corresponds to RawDataType used in the ZGY file format.
 * Data may be read and written using this type or 32-bit float.
 *
 * Thread safety: Enums do not have race conditions.
 */
enum class OPENZGY_API SampleDataType
{
  unknown = 1000,
  int8    = 1001,
  int16   = 1002,
  float32 = 1003,               // Note, cannot call it just "float".
};

namespace Formatters {
  extern OPENZGY_API std::ostream& operator<<(std::ostream& os, SampleDataType value);
}

/**
 * \brief Horizontal or vertical dimension as used in the public API.
 *
 * Horizontal dimension may be length or arc angle, although most
 * applications only support length. Vertical dimension may be time or
 * length. Vertical length is of course the same as depth. Arguably
 * there should have been separate enums for horizontal and vertical
 * dimension since the allowed values differ.
 *
 * Thread safety: Enums do not have race conditions.
 */
enum class OPENZGY_API UnitDimension
{
  unknown  = 2000,
  time     = 2001,
  length   = 2002,
  arcangle = 2003,
};

namespace Formatters {
  extern OPENZGY_API std::ostream& operator<<(std::ostream& os, UnitDimension value);
}

/**
 * \brief Possible algorithms to generate LOD bricks.
 *
 * We might trim this list later to what is actually in use.
 * The "classic" ZGY only uses the first two.
 *
 * Maps to LodAlgorithm in the implementation layer.
 *
 * Thread safety: Enums do not have race conditions.
 */
enum class OPENZGY_API DecimationType
{
  LowPass          = 100, ///< \brief Lowpass Z / decimate XY
  WeightedAverage  = 101, ///< \brief Weighted averaging (depends on global stats)
  Average          = 102, ///< \brief Simple averaging
  Median           = 103, ///< \brief Somewhat more expensive averaging
  Minimum          = 104, ///< \brief Minimum value
  Maximum          = 105, ///< \brief Maximum value
  //MinMax         = 106, ///< \brief Checkerboard of minimum and maximum values
  Decimate         = 107, ///< \brief Simple decimation, use first sample
  DecimateSkipNaN  = 108, ///< \brief Use first sample that is not NaN
  //DecimateRandom = 109, ///< \brief Random decimation using a fixed seed
  AllZero          = 110, ///< \brief Just fill the LOD brick with zeroes
  //WhiteNoise     = 111, ///< \brief Fill with white noise
  MostFrequent     = 112, ///< \brief The value that occurs most frequently
  MostFrequentNon0 = 113, ///< \brief The non-zero value that occurs most frequently
  AverageNon0      = 114, ///< \brief Average value, but treat 0 as NaN.
};

namespace Formatters {
  extern OPENZGY_API std::ostream& operator<<(std::ostream& os, DecimationType value);
}

/**
 * \brief Argument to _close() and _decimate()
 *
 * Currently only used in private methods of ZgyWriter.
 * The enum might end up being used in a public method later.
 *
 * Thread safety: Enums do not have race conditions.
 */
enum class OPENZGY_API FinalizeAction
{
  /**
   * Remove any existing information. If force=true this is unconditinal.
   * Otherwise it changes to "Keep" if the information already exists
   * and is not stale. With force=off, Delete and Keep actually do the same.
   */
  Delete = 3001,

  /**
   * Keep any existing information. If force=true this will even keep the
   * information if it is looks stale. Otherwise the mode changes to "Delete"
   * if OpenZGY cannot guarantee that the information is still correct.
   */
  Keep = 3002,

  /**
   * Do the minimum amount of work needed to bring the derived data up to
   * date. May cause the statistics and histogram to be less accurate.
   * Changes to "Keep" if the information is already up to date. changes
   * to "BuildFull" if information does not exist already or if the
   * file is compressed or if incremental cannot work for other reasons.
   * Setting force=true might do a little more work, to be defimed.
   */
  BuildIncremental = 3003,

  /**
   * Delete and re-create derived information from scratch. If force=true
   * this is unconditional. That may be useful if you want to change the
   * decimation algorithm. Otherwise change to "Keep" if the information
   * is already up to date. Never changes to BuildIncremental.
   */
  BuildFull = 3004,

  /**
   * As BuildFull, but do not collect or store histogram and statistics.
   * This might be useful to speed up writes that we know will only be
   * read by Petrel. Or other apps that don't strictly need them.
   * CAVEAT: This is incompatible with DecimationType::WeightedAverage.
   * So, indirectly this option will lower the quality of LoD data,
   */
  BuildNoHistogram = 3005,

  /**
   * Value used as default argument in finalize, inside close, etc.
   */
  BuildDefault = BuildFull,
};

namespace Formatters {
  extern OPENZGY_API std::ostream& operator<<(std::ostream& os, FinalizeAction value);
}

/**
 * \brief Statistics of all sample values on the file.
 *
 * \details Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_API SampleStatistics
{
public:
  std::int64_t cnt;  /**< \brief Number of added samples. */
  double sum;        /**< \brief Sum of added samples. */
  double ssq;        /**< \brief Sum-of-squares of added samples. */
  double min;        /**< \brief Minimum added sample value. */
  double max;        /**< \brief Maximum added sample value. */
  SampleStatistics()
    : cnt(0), sum(0), ssq(0), min(0), max(0)
  {
  }

  /** \brief Create a new instance with all contents filled in. */
  SampleStatistics(std::int64_t cnt_in,
                 double sum_in, double ssq_in,
                 double min_in, double max_in)
    : cnt(cnt_in), sum(sum_in), ssq(ssq_in), min(min_in), max(max_in)
  {
  }
};

/**
 * \brief Histogram of all sample values on the file.
 *
 * \internal
 * KEEP THIS TEXT IN SYNC WITH InternalZGY::HistogramData
 * \endinternal
 *
 * The histogram is described by the fixed total number of bins, the
 * center value of the samples in the first bin, and the center value
 * of the samples in the last bin.
 *
 * The width of each bin is given by (max - min) / (nbins - 1).
 * Bin 0 holds samples with values min_ +/-binwidth/2.
 *
 * This means that the total range of samples that can be represented in the
 * histogram is actually min-binwidth/2 to max+binwidth/2, which is slightly
 * larger than just min..max _unless_ each of the bins can only hold a single
 * value (e.g. data converted from 8-bit storage, in 256 bins).
 *
 * This definition has some subtle effects, best illustrated by a few examples.
 *
 * If the source data is ui8_t, the logical range is 0..255. Assume nbins=256.
 * This gives binwidth=1, and the true range of the histogram -0.5..+255.5.
 * But since the input is integral, the actual range is just 0..255 inclusive.
 * Try to fill the histogram with evenly distrubuted random data and you end
 * up with each bin having roughly the same number of elements.
 *
 * Now consider ui16_t, range 0..65535 and nbins is still 256. This gives
 * binwidth=257, not 256. The true range of the histogram is -128.5..+65663.5.
 * Try to fill the histogram with evenly distrubuted random data and you end
 * up with the first and the last bin having approximately half as many
 * elements as all the others. This is not really a problem, but may seem
 * a bit surprising.
 *
 * The range should have the zero-centric property. Zero, if present in the
 * range, should map to the center of a bin. Otherwise the histogram close to
 * zero might look odd.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_API SampleHistogram
{
public:
  std::int64_t samplecount; /**< \brief Sum of counts of all bins */
  double minvalue;          /**< \brief Center value of data in first bin */
  double maxvalue;          /**< \brief Center value of data in last bin */
  std::vector<std::int64_t> bins; /**< \brief array of counts by bin */
  SampleHistogram()
    : samplecount(0), minvalue(0), maxvalue(0), bins()
  {
  }

  /** \brief Create a new instance with all contents filled in. */
  SampleHistogram(std::int64_t samplecount_in,
                double minvalue_in, double maxvalue_in,
                const std::vector<std::int64_t>& bins_in)
    : samplecount(samplecount_in)
    , minvalue(minvalue_in)
    , maxvalue(maxvalue_in)
    , bins(bins_in)
  {
  }
};

/**
 * \brief Return value from filestats().
 * \details
 *
 * Meta information about the file that might be reported to a curious
 * end user. Such as the compression factor that was achieved. Or whether
 * this file was written by a legacy application that chose to store
 * alpha tiles in the file. The information is not needed to consume
 * the data stored in the file.
 *
 * This is a concrete class, not an interface. This is done to make
 * the class a POD and ensure it is copyable. Note this means you
 * should not assign the filestats() result to a referece.
 *
 * Note that the reported compression is about file sizes, not the
 * requested signal to noise factor. The latter would need to be
 * saved explicitly in the file. Currently it isn't.
 *
 * Normally the only wasted size in a ZGY file is the padding between
 * the header area and the first brick. For a file with aplha tiles
 * (cannot currently be cretated by OpenZGY) there might be more data
 * lost due to alignment. And files written to cloud storage in a way
 * that updates the same brick more than once may end up with holes.
 *
 * In the statistics, the padding after the header area is captured
 * in padding_size and all other holes are captured in wasted_size.
 * Caveat: Leftover space after the end of a compressed block will
 * currently not be detected.
 */
class OPENZGY_API FileStatistics
{
 // The code that fills in the contents do so by accessing the data fields
 // directly. There isn't much point in the additional boilerplace code
 // to implement write accessors.
 friend class Impl::ZgyMeta;
 friend class Impl::ZgyReader;
 friend class Impl::ZgyWriter;

private:
  std::int64_t _file_version;
  std::int64_t _file_size;
  std::int64_t _header_size;
  std::int64_t _data_start;
  std::vector<std::int64_t> _segment_sizes;
  //std::int64_t _padding_size;
  //std::int64_t _wasted_size;
  std::int64_t _alpha_normal_count;
  std::int64_t _alpha_normal_size_per_entry;
  std::int64_t _alpha_compressed_count;
  std::int64_t _alpha_compressed_size;
  std::int64_t _alpha_missing_count;
  std::int64_t _alpha_constant_count;
  std::int64_t _brick_normal_count;
  std::int64_t _brick_normal_size_per_entry;
  std::int64_t _brick_compressed_count;
  std::int64_t _brick_compressed_size;
  std::int64_t _brick_missing_count;
  std::int64_t _brick_constant_count;
  // Derived information
  std::int64_t _used_size;
  std::int64_t _used_if_uncompressed;
  double _compression_factor;
  bool _is_compressed;

public:
FileStatistics()
    : _file_version(0)
    , _file_size(0)
    , _header_size(0)
    , _data_start(-1)
    , _segment_sizes()
    //, _padding_size(0)
    //, _wasted_size(0)
    , _alpha_normal_count(0)
    , _alpha_normal_size_per_entry(0)
    , _alpha_compressed_count(0)
    , _alpha_compressed_size(0)
    , _alpha_missing_count(0)
    , _alpha_constant_count(0)
    , _brick_normal_count(0)
    , _brick_normal_size_per_entry(0)
    , _brick_compressed_count(0)
    , _brick_compressed_size(0)
    , _brick_missing_count(0)
    , _brick_constant_count(0)
    , _used_size(0)
    , _used_if_uncompressed(0)
    , _compression_factor(1.0)
    , _is_compressed(false)
  {
  }

  /// Version number from the main file header.
  std::int64_t fileVersion() const { return _file_version; }
  /// Total size of the file on disk or cloud.
  std::int64_t fileSize() const { return _file_size; }
  /// Size of all headers.
  std::int64_t headerSize() const { return _header_size; }
  /// Lowest address of any brick or tile, or -1 if there are none.
  std::int64_t dataStart() const { return _data_start; }
  /// Used for cloud storage only.
 const std::vector<std::int64_t>& segmentSizes() const {return _segment_sizes;}
  // Wasted due to first brick alignment.
  //std::int64_t paddingSize() const { return _padding_size; }
  // Wasted due to other reasons.
  //std::int64_t wastedSize() const { return _wasted_size; }
  /// Number of uncompressed tiles.
  std::int64_t alphaNormalCount() const { return _alpha_normal_count; }
  /// Size used by one uncompressed tile.
  std::int64_t alphaNormalSizePerEntry() const { return _alpha_normal_size_per_entry; }
  /// Number of compressed tiles.
  std::int64_t alphaCompressedCount() const { return _alpha_compressed_count; }
  /// Total size used by compressed tiles.
  std::int64_t alphaCcompressedSize() const { return _alpha_compressed_size; }
  /// Number of compressed tiles.
  std::int64_t alphaMissingCount() const { return _alpha_missing_count; }
  /// Number of constant value tiles.
  std::int64_t alphaConstantCount() const { return _alpha_constant_count; }
  /// Number of uncompressed bricks.
  std::int64_t brickNormalCount() const { return _brick_normal_count; }
  /// Size used by one uncompressed brick.
  std::int64_t brickNormalSizePerEntry() const { return _brick_normal_size_per_entry; }
  /// Number of compressed bricks.
  std::int64_t brickCompressedCount() const { return _brick_compressed_count; }
  /// Total size used by compressed bricks.
  std::int64_t brickCompressedSize() const { return _brick_compressed_size; }
  /// Number of compressed bricks.
  std::int64_t brickMissingCount() const { return _brick_missing_count; }
  /// Number of constant value bricks.
  std::int64_t brickConstantCount() const { return _brick_constant_count; }
  // Derived information
  /**
   * \brief Space used by headers and data bricks.
   * \details
   * File size not including padding and holes, combining all LOD
   * levels and the main headers. The padding between the header area
   * and the first brick will not be included. Nor will holes between
   * uncompressed bricks contribute. Currently any holes between
   * compressed bricks are not detected which means that they will be
   * counted as used. This can be derived from the other information.
   */
  std::int64_t usedSize() const { return _used_size; }
  /**
   * \brief Space needed if the file is/was uncompressed.
   * \details
   * As used_size if the file is/was uncompressed. This can be derived
   * from the other information.
   */
  std::int64_t usedIfUncompressed() const { return _used_if_uncompressed; }
  /**
   * \brief Measure how successful the compression was.
   * \details
   * Estimate the relative size of this possibly compressed file
   * compared to the same file if uncompressed. Will be 1.0 if file is
   * already uncompressed but a value of 1.0 doesn't technically imply
   * that the file is not compressed. Padding is ignored so the result
   * will not match precisely what you get by uncompressing the file
   * and storing it on disk. Also not taken into account is that the
   * padding after the header area might differ between the compressed
   * and uncompressed formats.
   */
  double compressionFactor() const { return _compression_factor; }
  /**
   * True if at least one brick is flagged as compressed, even in the
   * unlikely case where the compression didn't actually reduce the
   * file size. This can be derived from the other information.
   */
  bool isCompressed() const { return _is_compressed; }
  /**
   * For debugging. Output most of the information to the supplied ostream.
   */
  void dump(std::ostream& out, const std::string& prefix = "") const {
    std::stringstream segs;
    for (std::int64_t it : _segment_sizes)
      segs << " " << it;
    out << prefix << "ZGY version " << _file_version
        << " file compressed to "
        << int(100.0 * _compression_factor) << "% of original\n"
        << prefix << "Size:  "
        << _file_size << " bytes of which "
        << _header_size << " are in headers and "
        << _file_size - _used_size << " wasted\n"
        << prefix << "Segments:" << segs.str() << ", "
        << "Data area starts at: " << _data_start << "\n"
        << prefix << "Alpha: "
        << _alpha_missing_count    << " missing, "
        << _alpha_constant_count   << " constant, "
        << _alpha_normal_count     << " normal ("
        << _alpha_normal_count * _alpha_normal_size_per_entry << " bytes), "
        << _alpha_compressed_count << " compressed ("
        << _alpha_compressed_size  << " bytes)\n"
        << prefix << "Brick: "
        << _brick_missing_count    << " missing, "
        << _brick_constant_count   << " constant, "
        << _brick_normal_count     << " normal ("
        << _brick_normal_count * _brick_normal_size_per_entry << " bytes), "
        << _brick_compressed_count << " compressed ("
        << _brick_compressed_size  << " bytes)\n";
  }
};

/**
 * \brief Argument package for creating a ZGY file.
 *
 * This is a short lived helper class for passing arguments to the functions
 * that create a ZGY file. Needed because there are a lot of arguments
 * and C++ doesn't allow keyword arguments. Do NOT use this class
 * for holding on to the information.
 *
 * Information not set explicitly will be defaulted. The following
 * two statements give the same result:
 * \code
 *   ZgyWriterArgs default_args = ZgyWriterArgs();
 *
 *   ZgyWriterArgs default_args = ZgyWriterArgs()
 *     .filename("")
 *     .iocontext(nullptr)
 *     .compressor(nullptr)
 *     .lodcompressor(nullptr)
 *     .size(0, 0, 0)
 *     .bricksize(64, 64, 64)
 *     .datatype(SampleDataType::float32)
 *     .datarange(0, -1)
 *     .zunit(UnitDimension::unknown, "", 1.0)
 *     .hunit(UnitDimension::unknown, "", 1.0)
 *     .ilstart(0)
 *     .ilinc(0)
 *     .xlstart(0)
 *     .xlinc(0)
 *     .zstart(0)
 *     .zinc(0)
 *     .corners(ZgyWriterArgs::corners_t{0,0,0,0,0,0,0,0});
 * \endcode
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_API ZgyWriterArgs
{
public:
  typedef double float64_t;
  typedef std::array<std::array<float64_t,2>,4> corners_t;
  // Type aliases duplicated from the InternalZGY layer.
  typedef std::pair<std::shared_ptr<const void>, std::int64_t> rawdata_t;
  typedef std::function<rawdata_t(const rawdata_t&,const std::array<int64_t,3>&)> compressor_t;

private:
  friend class Impl::ZgyWriter; // Processing the high level parts.
  friend class ::Test::ZgyWriterMock;
  friend class Impl::EnumMapper; // Mapping the rest to the internal struct.
  friend class IZgyWriter; // For IZgyWriter::reopen()
  std::string _filename;
  std::shared_ptr<const IOContext> _iocontext; // Short lived but I'll use a smart pointer just in case.
  compressor_t _compressor;
  compressor_t _lodcompressor;
  std::array<std::int64_t,3> _size;
  std::array<std::int64_t,3> _bricksize;
  SampleDataType _datatype;
  std::array<float,2> _datarange;
  UnitDimension _zunitdim, _hunitdim;
  std::string _zunitname, _hunitname;
  double _zunitfactor, _hunitfactor;
  float _zstart, _zinc;
  std::array<float,2> _annotstart, _annotinc;
  corners_t _corners;
  // Keeps track of what information has been changed from the default.
  // Currently metafrom() will set all to true. That might change.
  // The information is only used by merge(). Transient information
  // such as iocontext is ignored.
  bool _have_size;
  bool _have_bricksize;
  bool _have_datatype;
  bool _have_datarange;
  bool _have_zunit;
  bool _have_hunit;
  bool _have_ilstart;
  bool _have_ilinc;
  bool _have_xlstart;
  bool _have_xlinc;
  bool _have_zstart;
  bool _have_zinc;
  bool _have_corners;

public:
  ZgyWriterArgs()
    : _filename("")
    , _iocontext()
    , _compressor()
    , _lodcompressor()
    , _size{0,0,0}
    , _bricksize{64,64,64}
    , _datatype(SampleDataType::float32)
    , _datarange{0, -1}
    , _zunitdim(UnitDimension::unknown)
    , _hunitdim(UnitDimension::unknown)
    , _zunitname("")
    , _hunitname("")
    , _zunitfactor(1.0)
    , _hunitfactor(1.0)
    , _zstart(0)
    , _zinc(0)
    , _annotstart{0,0}
    , _annotinc{0,0}
      , _corners{ {{0,0}} }
    , _have_size(false)
    , _have_bricksize(false)
    , _have_datatype(false)
    , _have_datarange(false)
    , _have_zunit(false)
    , _have_hunit(false)
    , _have_ilstart(false)
    , _have_ilinc(false)
    , _have_xlstart(false)
    , _have_xlinc(false)
    , _have_zstart(false)
    , _have_zinc(false)
    , _have_corners(false)
  {
  }

  /**
   * \brief Output in human readable form for debugging.
   */
  void dump(std::ostream& out) const
  {
    out << "ZgyWriterArgs\n"
        << "  filename:    \"" << _filename << "\"\n"
        << "  iocontext:   " << (_iocontext ? "*" : "(null)") << "\n"
        << "  compressor:  " << (_compressor ? "*" : "(null)") << "\n"
        << "  lodcompress: " << (_lodcompressor ? "*" : "(null)") << "\n"
        << "  " << (_have_size?"*":"") << "size:        (" << _size[0] << "," << _size[1] << "," << _size[2] << ")\n"
        << "  " << (_have_bricksize?"*":"") << "bricksize:   (" << _bricksize[0] << "," << _bricksize[1] << "," << _bricksize[2] << ")\n"
        << "  " << (_have_datatype?"*":"") << "datatype:    " << int(_datatype) << "\n"
        << "  " << (_have_datarange?"*":"") << "datarange:   " << _datarange[0] << " to " << _datarange[1] << "\n"
        << "  " << (_have_zunit?"*":"") << "zunit:       " << int(_zunitdim) << " \"" << _zunitname << "\" " << _zunitfactor << "\n"
        << "  " << (_have_hunit?"*":"") << "hunit:       " << int(_hunitdim) << " \"" << _hunitname << "\" " << _hunitfactor << "\n"
        << "  " << (_have_ilstart||_have_ilinc?"*":"") << "ilstart/inc: " << _annotstart[0] << " / " << _annotinc[0] << "\n"
        << "  " << (_have_xlstart||_have_xlinc?"*":"") << "xlstart/inc: " << _annotstart[1] << " / " << _annotinc[1] << "\n"
        << "  " << (_have_zstart||_have_zinc?"*":"") << "zstart/inc:  "  << _zstart << " / " << _zinc << "\n"
        << "  " << (_have_corners?"*":"") << "corner0:     " << _corners[0][0] << ", " << _corners[0][1] << "\n"
        << "  " << (_have_corners?"*":"") << "corner1:     " << _corners[1][0] << ", " << _corners[1][1] << "\n"
        << "  " << (_have_corners?"*":"") << "corner2:     " << _corners[2][0] << ", " << _corners[2][1] << "\n"
        << "  " << (_have_corners?"*":"") << "corner3:     " << _corners[3][0] << ", " << _corners[3][1] << "\n";
      }

  /**
   * \brief Set file to open.
   * \if SSTORE
   * \details If starting with sd:// this opens a file on seismic store.
   * \endif
   */
  ZgyWriterArgs& filename(const std::string& value) { _filename = value; return *this; }
  /**
   * \brief Set credentials and other configuration.
   * \details This depends on the back-end.
   * For local files you normally don't need to pass anything here.
   */
  ZgyWriterArgs& iocontext(const IOContext *value); // _iocontext = value->clone()
  /**
   * \brief Set functor for compressing full resolution data.
   */
  ZgyWriterArgs& compressor(const compressor_t& value) { _compressor = value; return *this; }
  /**
   * \brief Set functor for compressing full resolution data.
   * \details This overload uses a factory to look up the compressor.
   */
  ZgyWriterArgs& compressor(const std::string& name, const std::vector<std::string>& args);
  /**
   * \brief Set functor for compressing full resolution data.
   * \details This overload uses a factory to look up the ZGY compressor.
   * It is just a convenience that is shorter to type.
   */
  ZgyWriterArgs& zfp_compressor(float snr);
  /**
   * \brief Set functor for compressing low resolution data.
   */
  ZgyWriterArgs& lodcompressor(const compressor_t& value) { _lodcompressor = value; return *this; }
  /**
   * \brief Set functor for compressing low resolution data.
   * \details This overload uses a factory to look up the compressor.
   */
  ZgyWriterArgs& lodcompressor(const std::string& name, const std::vector<std::string>& args);
  /**
   * \brief Set functor for compressing low resolution data.
   * \details This overload uses a factory to look up the ZGY compressor.
   * It is just a convenience that is shorter to type.
   */
  ZgyWriterArgs& zfp_lodcompressor(float snr);
  /**
   * \brief Set size of the survey.
   * \param ni number of inlines (slowest varying index).
   * \param nj number of crosslines.
   * \param nk number of samples per trace (fastest).
   */
  ZgyWriterArgs& size(std::int64_t ni, std::int64_t nj, std::int64_t nk) { _size[0]=ni; _size[1]=nj; _size[2]=nk; _have_size = true; return *this; }
  /**
   * \brief Set size of one brick.
   * \details Almost always (64,64,64). Change at your own peril.
   */
  ZgyWriterArgs& bricksize(std::int64_t ni, std::int64_t nj, std::int64_t nk) { _bricksize[0]=ni; _bricksize[1]=nj; _bricksize[2]=nk; _have_bricksize = true; return *this; }
  /**
   * \brief Set type of samples in each brick.
   */
  ZgyWriterArgs& datatype(SampleDataType value) { _datatype = value; _have_datatype = true; return *this; }
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
ZgyWriterArgs& datarange(float lo, float hi) { _datarange[0] = lo; _datarange[1] = hi; _have_datarange = true; return *this; }
  /**
   * \brief Set vertical unit.
   * \param dimension time or depth (a.k.a. length).
   * \param name for annotation only.
   * \param factor multiply by this factor to convert from storage units to SI units.
   */
  ZgyWriterArgs& zunit(UnitDimension dimension, const std::string& name, double factor) {
    _zunitdim = dimension;
    _zunitname = name;
    _zunitfactor = factor;
    _have_zunit = true;
    return *this;
  }
  /**
   * \brief Set horizontal unit.
   * \param dimension cartesian length or (unsupported) polat coordinates,
   * \param name for annotation only.
   * \param factor multiply by this factor to convert from storage units to SI units.
   */
  ZgyWriterArgs& hunit(UnitDimension dimension, const std::string& name, double factor) {
    _hunitdim = dimension;
    _hunitname = name;
    _hunitfactor = factor;
    _have_hunit = true;
    return *this;
  }
  /// \brief Set first (ordinal 0) inline number.
  /// \details For maximum portability the inline and crossline start
  /// and increment should be integral numbers. Some applications
  /// might choose to convert them to int.
  ZgyWriterArgs& ilstart(float value) { _annotstart[0] = value; _have_ilstart = true; return *this; }
  /// \brief Set inline number increment between two adjacent ordinal values.
  /// \copydetails ilstart
  ZgyWriterArgs& ilinc(float value) { _annotinc[0] = value; _have_ilinc = true; return *this; }
  /// \brief Set first (ordinal 0) crossline number.
  /// \copydetails ilstart
  ZgyWriterArgs& xlstart(float value) { _annotstart[1] = value; _have_xlstart = true; return *this; }
  /// \brief Set crossline number increment between two adjacent ordinal values.
  /// \copydetails ilstart
  ZgyWriterArgs& xlinc(float value) { _annotinc[1] = value; _have_xlinc = true; return *this; }
  /// \brief Set first time/depth.
  /// \details Vertical annotation is generally safe to have non-integral.
  ZgyWriterArgs& zstart(float value) { _zstart = value; _have_zstart = true; return *this; }
  /// \brief Set increment (distance between samples) in vertical direction.
  /// \copydetails zstart
  ZgyWriterArgs& zinc(float value) { _zinc = value; _have_zinc = true; return *this; }
  /// \brief Set survey corner points in world coordinates.
  /// \details The corners are ordered origin, last inline (i.e. i=last, j=0),
  /// last crossline, diagonal.
  ZgyWriterArgs& corners(const corners_t& value) { _corners = value; _have_corners = true; return *this; }
  /**
   * \brief Copy metadata from existing file.
   *
   * Set most of the metadata from an open file descriptor.
   * Typically used when the file to be created is to be computed from
   * an existing file. Typically this will be called first so the
   * settings don't inadverently shadow something set explicitly.
   */
  ZgyWriterArgs& metafrom(const std::shared_ptr<OpenZGY::IZgyReader>&);

  /**
   * \brief Copy metadata from another ZgyWriterArgs.
   *
   * Copy only those settings that have been explicitly changed in the
   * supplied "other" ZgyWriterArgs into *this. If other.metafrom()
   * has been called it is unspecified what gets copied.
   */
  ZgyWriterArgs& merge(const ZgyWriterArgs&);

  // TODO-Low: Add accessors as well. But these are only for internal
  // use so the ugliness of accessing data members isn't that bad.
};

/**
 * \brief Base class of IZgyReader and IZgyWriter.
 * \details Thread safety: Interfaces do not have race conditions.
 */
class OPENZGY_API IZgyMeta
{
public:
  typedef std::int8_t int8_t;
  typedef std::int16_t int16_t;
  typedef std::int32_t int32_t;
  typedef std::int64_t int64_t;
  typedef float float32_t;
  typedef double float64_t;
  typedef std::array<int64_t,3> size3i_t;
  typedef std::array<std::array<float64_t,2>,4> corners_t;
  // Type aliases duplicated from the InternalZGY layer.
  typedef std::pair<std::shared_ptr<const void>, std::int64_t> rawdata_t;
  typedef std::function<rawdata_t(const rawdata_t&,const std::array<int64_t,3>&)> compressor_t;

public:
  virtual ~IZgyMeta();
  virtual size3i_t       size()           const = 0; /**< \brief Size in inline, crossline, vertical directions. */
  virtual SampleDataType datatype()       const = 0; /**< \brief Type of samples in each brick. */
  virtual std::array<float32_t,2> datarange() const = 0; /**< \brief Used for float to int scaling. */
  virtual std::array<float32_t,2> raw_datarange() const = 0; /**< \brief datarange before adjustment. */
  virtual UnitDimension  zunitdim()       const = 0; /**< \brief Vertical dimension. */
  virtual UnitDimension  hunitdim()       const = 0; /**< \brief Horizontal dimension. */
  virtual std::string    zunitname()      const = 0; /**< \brief For annotation only. Use hunitfactor, not the name, to convert to or from SI. */
  virtual std::string    hunitname()      const = 0; /**< \brief For annotation only. Use hunitfactor, not the name, to convert to or from SI. */
  virtual float64_t      zunitfactor()    const = 0; /**< \brief Multiply by this factor to convert from storage units to SI units. */
  virtual float64_t      hunitfactor()    const = 0; /**< \brief Multiply by this factor to convert from storage units to SI units. */
  virtual float32_t      zstart()         const = 0; /**< \brief First time/depth. */
  virtual float32_t      zinc()           const = 0; /**< \brief Increment in vertical direction. */
  virtual std::array<float32_t,2> annotstart() const = 0; /**< \brief First inline, crossline. */
  virtual std::array<float32_t,2> annotinc()   const = 0; /**< \brief Increment in inline, crossline directions. */
  virtual const corners_t corners()      const = 0; /**< \brief Survey corner points in world coordinates. */
  virtual const corners_t indexcorners() const = 0; /**< \brief Survey corner points in ordinal (i,j) coordinates. */
  virtual const corners_t annotcorners() const = 0; /**< \brief Survey corner points in inline, crossline coordinates. */
  virtual size3i_t       bricksize()      const = 0; /**< \brief Size of one brick. Almost always (64,64,64), change at your own peril. */
  virtual std::vector<size3i_t> brickcount() const = 0; /**< \brief Number of bricks at each resolution (LOD) level. */
  virtual int32_t        nlods()          const = 0; /**< \brief Number of resolution (LOD) levels. */
  // Only expose the guids we think the application will need.
  // Note that OpenZGY doesn't really support updates,
  // so dataid() and previd() are not very useful.
  //virtual std::string dataid() const = 0;    /**< GUID set on file creation. */
  virtual std::string verid()  const = 0;    /**< GUID set each time the file is changed. */
  //virtual std::string previd() const = 0;    /**< GUID before last change. */

  // The Python version has meta() as a dict holding all the meta data,
  // this isn't really useful in C++ and just makes it harder to see which
  // metadata is being used. [set_]numthreads is N/A in this accessor.
  // If the code allows enabling multi threading then this would be configured
  // in IOContext or as arguments to the compress plugin.
  //virtual void           meta()           const = 0;
  //virtual int32_t        numthreads()     const = 0;
  //virtual void           set_numthreads(int32_t) = 0;
  virtual void           dump(std::ostream&) const = 0; /**< \brief Output in human readable form for debugging. */
  virtual SampleStatistics statistics()     const = 0; /**< \brief Statistics of all sample values on the file. */
  virtual SampleHistogram  histogram()      const = 0; /**< \brief Histogram of all sample values on the file. */
  virtual std::shared_ptr<const FileStatistics> filestats() const = 0; /**< \brief For display purposes only. */
};

/**
 * \brief Base class of IZgyReader and IZgyWriter.
 * \details Thread safety: Interfaces do not have race conditions.
 */
class OPENZGY_API IZgyTools : virtual public IZgyMeta
{
public:
  virtual ~IZgyTools();
  /**
   * \brief General coordinate conversion of an array of points.
   * (NOT IMPLEMENTED YET)
   *
   * \param from: control points in the current coordinate system.
   * \param to:   control points in the desired coordinate system.
   *
   * \details Convert coordinates in place, with the conversion
   * defined by giving the values of 3 arbitrary control points in
   * both the "from" and "to" coordinate system. A common choice of
   * arbitrary points is to use three of the lattice corners.
   * As a convenience the "from" and "to" parameters are declared as
   * corners_t so the caller can pass corners(), annotcorners(), or
   * indexcorners() directly.
   *
   * \internal TODO-Low I haven't implemented the array version yet,
   * only the one that converts a single point.
   */
  virtual void transform(const corners_t& from, const corners_t& to, std::vector<std::array<float64_t,2>>&) const = 0;
  /**
   * \brief General coordinate conversion of a single coordinate pair.
   * \copydetails IZgyTools::transform()
   */
  virtual std::array<float64_t,2> transform1(const corners_t& from, const corners_t& to, const std::array<float64_t,2>&) const = 0;
  virtual std::array<float64_t,2> annotToIndex(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
  virtual std::array<float64_t,2> annotToWorld(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
  virtual std::array<float64_t,2> indexToAnnot(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
  virtual std::array<float64_t,2> indexToWorld(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
  virtual std::array<float64_t,2> worldToAnnot(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
  virtual std::array<float64_t,2> worldToIndex(const std::array<float64_t,2>&) const = 0; /**< \brief Convert a single coordinate pair. */
};

/**
 * \brief Main API for reading ZGY files.
 *
 * Obtain a concrete instance by calling the factory method IZgyReader::open().
 * You can then use the instance to read both meta data and bulk data.
 * It is recommended to explicitly close the file when done with it.
 *
 * Note: Yes, I understand that the open() factory method should not have
 * been lexically scoped inside the ostensibly pure %IZgyReader interface.
 * But this reduces the number of classes a library user needs to relate to.
 *
 * Thread safety: Interfaces do not have race conditions.
 */
class OPENZGY_API IZgyReader : virtual public IZgyTools
{
public:
  virtual ~IZgyReader();
  /// \copydoc Impl::ZgyReader::read(const size3i_t&,const size3i_t&,float*,int) const
  virtual void read(const size3i_t& start, const size3i_t& size, float* data, int lod = 0) const = 0;
  /// \copydoc Impl::ZgyReader::read(const size3i_t&,const size3i_t&,std::int16_t*,int) const
  virtual void read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod = 0) const = 0;
  /// \copydoc Impl::ZgyReader::read(const size3i_t&,const size3i_t&,std::int8_t*,int) const
  virtual void read(const size3i_t& start, const size3i_t& size, std::int8_t* data, int lod = 0) const = 0;
  /// \copydoc Impl::ZgyReader::readconst()
  virtual std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, int lod = 0, bool as_float = true) const = 0;
  /// \copydoc Impl::ZgyReader::close()
  virtual void close() = 0;
  /// \brief Open a ZGY file for reading.
  static std::shared_ptr<IZgyReader> open(const std::string& filename, const IOContext* iocontext = nullptr);
};

/**
 * \brief Main API for creating ZGY files.
 *
 * Obtain a concrete instance by calling the factory method IZgyWriter::open().
 * All meta data is specified in the call to open(), so meta data will appear
 * to be read only. You can use the instance to write bulk data. The file
 * becomes read only once the instance is closed.
 *
 * It is recommended to call finalize() after all bulk has been written.
 * but if you forget this will be called from close() with default arguments.
 * It is required to explicitly close() the file when done with it,
 * Technically the destructor will call close() but being a destructor
 * it needs to swallow any exceptions.
 *
 * Note: Yes, I understand that the open() factory method should not have
 * been lexically scoped inside the ostensibly pure %IZgyWriter interface.
 * But this reduces the number of classes a library user needs to relate to.
 *
 * Thread safety: Interfaces do not have race conditions.
 */
class OPENZGY_API IZgyWriter : virtual public IZgyTools
{
public:
  virtual ~IZgyWriter();
  /// \copydoc Impl::ZgyWriter::read(const size3i_t&,const size3i_t&,float*) const
  virtual void read(const size3i_t& start, const size3i_t& size, float* data) const = 0;
  /// \copydoc Impl::ZgyWriter::read(const size3i_t&,const size3i_t&,std::int16_t*) const
  virtual void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const = 0;
  /// \copydoc Impl::ZgyWriter::read(const size3i_t&,const size3i_t&,std::int8_t*) const
  virtual void read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const = 0;
  /// \copydoc Impl::ZgyReader::readconst()
  virtual std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float = true) const = 0;
  /// \copydoc Impl::ZgyReader::close()
  virtual void write(const size3i_t& start, const size3i_t& size, const float* data) = 0; /**< \copydoc Impl::ZgyWriter::write(const size3i_t&,const size3i_t&,const float*) */
  virtual void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data) = 0; /**< \copydoc Impl::ZgyWriter::write(const size3i_t&,const size3i_t&,const std::int16_t*) */
  virtual void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data) = 0; /**< \copydoc Impl::ZgyWriter::write(const size3i_t&,const size3i_t&,const std::int8_t*) */
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const float* data) = 0; /**< \copydoc Impl::ZgyWriter::writeconst(const size3i_t&,const size3i_t&,const float*) */
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) = 0; /**< \copydoc Impl::ZgyWriter::writeconst(const size3i_t&,const size3i_t&,const std::int16_t*) */
  virtual void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data) = 0; /**< \copydoc Impl::ZgyWriter::writeconst(const size3i_t&,const size3i_t&,const std::int8_t*) */
  /// \copydoc Impl::ZgyWriter::finalize()
  virtual void finalize(
               const std::vector<DecimationType>& decimation = std::vector<DecimationType>(),
               const std::function<bool(std::int64_t,std::int64_t)>& progress = nullptr,
               FinalizeAction action = FinalizeAction::BuildDefault,
               bool force = false) = 0;
  /// \copydoc Impl::ZgyWriter::close_incomplete()
  virtual void close_incomplete() = 0;
  /// \copydoc Impl::ZgyWriter::close()
  virtual void close() = 0;
  /// \cond IMPL
  virtual bool errorflag() const = 0;
  virtual void set_errorflag(bool) = 0;
  /// \endcond
  /// \brief Create a ZGY file and open it for writing.
  static std::shared_ptr<IZgyWriter> open(const ZgyWriterArgs& args);
  /// \brief Open an existing ZGY file for writing.
  static std::shared_ptr<IZgyWriter> reopen(const ZgyWriterArgs& args);
};

/**
 * \brief Operations other than read and write.
 *
 * Any operations that don't fit into IZgyReader or IZgyWriter go here.
 * Such as deleting a file. Or any other operation that does not need
 * the file to be open first.
 *
 * Thread safety: Interfaces do not have race conditions.
 */
class OPENZGY_API IZgyUtils
{
public:
  virtual ~IZgyUtils();
  /**
   * \brief Delete a file. Works both for local and cloud files.
   *
   * Note that the instance must be of the correct (local or cloud) type.
   */
  virtual void deletefile(const std::string& filename, bool missing_ok=true) = 0;
  /**
   * \brief Return a url that might be opened faster.
   *
   * Using this method might help performance in a distributed system
   * where a large number of processing nodes need to open the same
   * file for read. One node would call this method which costs about
   * the same as opening the file locally. All other nodes can then
   * use the alternate url which opens a lot faster due to the meta
   * data having already been retrieved and serialized in the alturl.
   *
   * The alternate url is valid for a limited time only. Typically it
   * expires when the current access token expires but it might expire
   * earlier.
   *
   * Note that the alternate url might contain secrets, so it needs to
   * be protected just as you would protect an access token.
   *
   * The methid is not useful for local files where the file name will
   * just be returned unmodified. The method is unlikely to be useful
   * in a single-node application. In that case the application can
   * simply keep the file open instead of keeping the alturl
   * reference.
   */
  virtual std::string alturl(const std::string& filename) = 0;
  /*
   * Return an IDToken or the equivalent for the current back end.
   * This might be useful when sharing credentials between files.
   */
  virtual std::string idtoken() = 0;
  /**
   * \brief Create a new concrete instance of IZgyUtils.
   * \param prefix File name or file name prefix.
   * \param iocontext Credentials and other configuration.
   *
   * The reason you need to supply a file name or a file name prefix is that
   * you need to provide enough information to identify the back-end that
   * \if SSTORE
   * this instance will be bound to. So both "sd://some/bogus/file.zgy" and
   * just "sd://" will produce an instance that works for the seismic store.
   * \else
   * this instance will be bound to. So if you have registered a back-end
   * named "xx", both "xx://some/bogus/file.zgy" and just "xx://" will
   * produce an instance that works for your XX backend,
   * \endif
   *
   * For performance reasons you should consider caching one IZgyUtils
   * instance for each back end you will be using. Instead of just creating
   * a new one each time you want to invoke a method. Just remember that
   * most operations need an instance created with the same prefix.
   */
  static std::shared_ptr<IZgyUtils> utils(const std::string& prefix, const IOContext* iocontext);
};

/**
 * \brief Simple progress bar.
 *
 * Progress bar that writes dots (51 by default) to standard output.
 * This can be user as-is for simple command line apps, or you can use
 * the source code as an example on how to write your own.
 *
 * The default of 51 dots will print one dot at startup and then one
 * additional dot for each 2% work done.
 *
 * If you are using this to write to the cloud a file that is smaller
 * than ~10 GB then the progress bar will probably move in larger
 * jumps. Because writing to a cloud back-end uses very large buffers.
 * Most cloud back-ends cannot report progress inside a "write block".
 *
 * When passing a progress reporter to a function, make sure you do not
 * pass the class itself. You need to create an instance of it.
 *
 * When passed as the second parameter to IZgyWriter::finalize() you will
 * need to pass it as std::ref(p) since the instance is noncopyable.
 * finalize() expects a reference to a function. The copy is inserted by
 * the compiler when converting ProgressWithDots::operator() into an
 * std::function. If OpenZGY were to add the std::ref itself the API would
 * get messier. Adding some ProgressWithDots::make_functor_with_ref()
 * method would work but hardly be any better.
 *
 * Thread safety:
 * Protected by a mutex.
 *
 * TODO-Low: Add a pimpl pattern to avoid needing the <mutex> header.
 * I really want to keep the include count low in this main header file.
 */
class OPENZGY_API ProgressWithDots
{
private:
  int _dots_printed;
  int _length;
  std::ostream& _outfile;
  std::mutex _mutex;

private:
  ProgressWithDots(const ProgressWithDots&) = delete;
  ProgressWithDots& operator=(const ProgressWithDots&) = delete;

public:
  /**
   * \param length Size of progress bar, default 51 dots.
   * \param outfile Stream to write output, default std::cerr
   *
   * \copydoc OpenZGY::ProgressWithDots
   */
  ProgressWithDots(int length=51, std::ostream& outfile = std::cerr);
  /**
   * \brief Callback invoked to report progress.
   *
   * \param done Number of work units done.
   * \param total Number of work units in total.
   *
   * The callback will normally get called exactly once with done==0,
   * before processing starts but after "total" is known.
   * And exactly once with done==total signifying that the work is
   * finished and the function being monitored should soon return.
   *
   * This particular callback will always return true,
   * meaning that the operation is not to be aborted.
   */
  bool operator()(std::int64_t done, std::int64_t total);
};

} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
