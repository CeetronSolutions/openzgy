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

#ifdef HAVE_ZFP // Entire file

#include "compression.h"
#include "fancy_timers.h"
#include "environment.h"
#include "exception.h"

#include <zfp.h>

#include <string>
#include <cstring>
#include <vector>
#include <memory>
#include <functional>
#include <algorithm>
#include <cmath>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \brief Implement ZFP compression.
 * \details The methods in this class will called from a factory. \details
 *
 * If you are comparing this file with the Python version you may see
 * they are rather different. The Python code was used to evaluate
 * ZFP performance and to experiment with different compression settings.
 * So it contains code for telemetry and there are remnants of code
 * allowing to use ZFP in different ways.
 *
 * On the other hand, the Python bindings for ZFP are way simpler to
 * use than the C++ api.
 *
 * Thread safety: The plug-ins and the compression algorithms they call
 * are all supposed to be thread safe. Registering the factory currently
 * is not.
 */
class ZfpCompressPlugin
{
public:

  /** \copydoc NullCompressPlugin::getCompressor */
  static compressor_t getCompressor(const std::vector<std::string>& args)
  {
    if (args.size() != 2 || args[0] != "snr" || args[1].empty())
      throw OpenZGY::Errors::ZgyUserError("Bad arguments to ZFP compress.");
    std::size_t idx;
    float snr = std::stof(args[1], &idx);
    if (idx != args[1].size())
      throw OpenZGY::Errors::ZgyUserError("ZFP \"snr\" argument not a number.");
    return [snr](const rawdata_t& data, const index3_t& shape) {
             return compress(data, shape, snr);
           };
  }

  /**
   * \copydoc NullCompressPlugin::compress
   * \param snr Desired signal to noise ratio.
   */
  static rawdata_t compress(const rawdata_t& data, const index3_t& shape, float snr)
  {
    static SummaryPrintingTimerEx timer("ZFP-compress");
    SimpleTimerEx tt(timer);
    return _try_compress(data, shape, snr);
  }

  /** \copydoc NullCompressPlugin::decompress */
  static rawdata_t decompress(const rawdata_t& cdata,
                              const BrickStatus& status,
                              const index3_t& shape)
  {
    // ZFP encodes size in its own header so we ignore what is passed
    // by the caller.
    if (status != BrickStatus::Compressed ||
        cdata.second < 4 ||
        0 != std::strncmp((const char*)cdata.first.get(), "zfp", 3))
      return rawdata_t{nullptr, 0}; // Not ours.
    static SummaryPrintingTimerEx timer("ZFP-decompress");
    SimpleTimerEx tt(timer);
    // TODO-Worry: Is ZFP really handling byte swapping?
    return _try_decompress(cdata, shape);
  }

  /** \copydoc NullCompressPlugin::Register */
  class Register
  {
  public:
    Register()
    {
      CompressFactoryImpl::registerCompressor
        ("ZFP", &ZfpCompressPlugin::getCompressor);
      CompressFactoryImpl::registerDecompressor
        ("Zfp", &ZfpCompressPlugin::decompress);
    }
  };

private:
  static int _zfp_precision_from_snr(rawdata_t data, float want_snr)
  {
    if (want_snr < 10 || want_snr > 70)
      return 0;

    // TODO-Performance, do I really need to check?
    // Maybe the file could have a global "may contain non-finite" ?
    const float *ptr = static_cast<const float*>(data.first.get());
    const std::int64_t size = data.second / sizeof(float);
    if (size != std::count_if(ptr, ptr + size, [](float x){return std::isfinite(x);}))
      return 0; // found an Inf or a NaN. Try lossless later.

    // ZFP precision has a valid range of 1 to 64 bit planes.
    // If compression is done simply by reducing the number of bits:
    // If signal uses 8 bits (-128..+127) the average magnitude of the
    // signal would be 64 and the acerage quantization noise 0.5.
    // So SNR with the current metric would be 20*log10(128) ~= 42 dB.
    // Assume ZFP is better than that if asked to use N bit planes:
    // 48 dB at 8 bits, or simply 6 dB per bit.
    // A heuristic taken from hca_filt.zgy gets a somewhat different
    // result, but is quite accurate for that data set.

    return (int(want_snr) + 23) / 5;
  };

  /**
   * \brief ZFP compression.
   *
   * Uses "precision" more if want_snr>0, or lossless otherwise.
   */
  static rawdata_t _zfp_try_one_compress(const rawdata_t& idata, const index3_t& shape, float want_snr)
  {
    if (idata.second != shape[0]*shape[1]*shape[2]*(std::int64_t)sizeof(float))
      throw OpenZGY::Errors::ZgyInternalError("ZFP got wrong size hint");

    const int precision = _zfp_precision_from_snr(idata, want_snr);
    const zfp_type type = zfp_type_float; // array scalar type

    // TODO-Low Am I expcted to transpose the data?
    // It might not make any practical difference as long as it is
    // done consistently. More testing would be nice though.
    // Allocate meta data for the 3D uncompressed array a[nz][ny][nx]
    // Note that the FIRST size argument is the fastest varying one.
    // TODO-Low: Check at runtime that shape[n] fits in an integer.
    zfp_field* field = zfp_field_3d(
        const_cast<void*>(idata.first.get()),
        type,
        static_cast<unsigned>(shape[2]), // nx
        static_cast<unsigned>(shape[1]),
        static_cast<unsigned>(shape[0]));

    // Allocate meta data for a compressed stream. Bitstream attached later.
    zfp_stream* zfp = zfp_stream_open(nullptr);

    // Set "precision" compression mode with chosen quality, or leave lossless.
    if (precision > 0)
      zfp_stream_set_precision(zfp, precision);

    // Allocate buffer for compressed data. Plus fudge for header.
    const size_t cdata_size = 100+zfp_stream_maximum_size(zfp, field);
    if (cdata_size <= 100) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      return rawdata_t{nullptr,0};
    }
    // The explicit destructor is there to prevent valgrind from
    // complaining about delete used where delete[] was indicated. In
    // C++17 or thereabouts just use std::shared_ptr<char[]>. Older
    // compilers allow that trick with unique_ptr but not shared_ptr.
    auto cdata_buffer = std::shared_ptr<char>
      (new char[cdata_size], std::default_delete<char[]>());

    // Associate a new bit stream for compressed data with allocated buffer.
    bitstream* bits = stream_open(cdata_buffer.get(), cdata_size);
    zfp_stream_set_bit_stream(zfp, bits);
    zfp_stream_rewind(zfp);

    // Compress array and output compressed stream.
    // Field has an internal pointer to the uncompressed "array".
    // zfp is associated with stream which in turn knows about cdata_buffer.
    // Note that hdrsize is in bits, zfpsize is in bytes including header.
    const size_t hdrsize = zfp_write_header(zfp, field, ZFP_HEADER_FULL);
    const size_t zfpsize = zfp_compress(zfp, field); // return actual size

    /* clean up */
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(bits);

    // zfpsize might be smaller than allocated data in cdata_buffer, that is ok.
    // But if zero this means the compression failed.
    // zfpsize does include hdrsize (zfp_compress returns current offset).
    return hdrsize && zfpsize ? rawdata_t{cdata_buffer, zfpsize} : rawdata_t{nullptr,0};
  }

  /**
   * \brief Undo ZFP compression in either precision mode or lossless mode.
   *
   * The expected result size is passed in and will be checked against
   * the header of the compressed data. The expected type is "float";
   * this will also be checked.
   */
  static rawdata_t _zfp_try_one_decompress(const rawdata_t& cdata, const index3_t& size)
  {
    rawdata_t result{nullptr,0};
    const zfp_type type = zfp_type_float; // array scalar type
    zfp_field* field = zfp_field_alloc();
    // TODO-Low Do I need to transpose or alternatively set strides?
    // It might not make any practical difference as long as it is
    // done consistently. More testing would be nice though.

    // Allocate meta data for a compressed stream. Bitstream attached later.
    zfp_stream* zfp = zfp_stream_open(nullptr);

    // Associate a new bit stream and associate it with the compressed data.
    // ZFP does not mind if the buffer contains more data than needed.
    bitstream* bits = stream_open(const_cast<void*>(cdata.first.get()), cdata.second);

    try {
      zfp_stream_set_bit_stream(zfp, bits);
      zfp_stream_rewind(zfp);

      const size_t header_size = zfp_read_header(zfp, field, ZFP_HEADER_FULL);
      if (header_size == 0)
        throw OpenZGY::Errors::ZgyFormatError("ZFP corrupted header");

      // Technically I don't need the ZFP_HEADER_META part of the header,
      // but it makes things easier for the Python version. Since I have
      // the data I might as well do a consistency check.
      size_t sizearray[4]{0}; // size[0] is nx, i.e. fastest varying index.
      zfp_field_size(field, &sizearray[0]);
      if (zfp_field_type(field) != type ||
          zfp_field_dimensionality(field) != 3 ||
          sizearray[2] != size[0] || // nz -> slowest varying
          sizearray[1] != size[1] || // ny
          sizearray[0] != size[2]) // nx
      {
        //std::cout << ".type " << (int)zfp_field_type(field)
        //        << ".dims " << zfp_field_dimensionality(field)
        //        << ".size " << sizearray[0]
        //        << ", " << sizearray[1]
        //        << ", " << sizearray[2]
        //        << std::endl;
        throw OpenZGY::Errors::ZgyFormatError("ZFP type or size miamatch");
      }

      const std::int64_t idata_size = size[0]*size[1]*size[2]*sizeof(float);
      auto idata_buffer = std::shared_ptr<char>
        (new char[idata_size], std::default_delete<char[]>());
      zfp_field_set_pointer(field, idata_buffer.get());

      // The function returns how many compressed bytes that it consumed.
      // Not important here, except that a zero return signals an error.
      // Note: Should I raise an exception on error? This isn't simply a
      // matter of "somebody else might be able to decompress this".
      const size_t cdata_size = zfp_decompress(zfp, field);
      if (cdata_size == 0)
        throw OpenZGY::Errors::ZgyFormatError("ZFP corrupted body");
      result.first = idata_buffer;
      result.second = idata_size;
    }
    catch (...) {
      zfp_field_free(field);
      zfp_stream_close(zfp);
      stream_close(bits);
      throw;
    }

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(bits);

    return result;
  }

  /**
   * \brief Try ZFP compression in precision mode and lossless mode.
   *
   * The method needs to make sure the compressed size is not larger than
   * the input. Not only is this inefficient but there is code that
   * assumes the compressed version will be smaller.
   *
   * Note that the maximum size is the size this data will take up on
   * disk if left uncompressed. Currently we enforce that compressed
   * files must be declared as float and that data fed to the
   * compressor is also float. As long as those rules hold there is
   * no need to worry.
   */
  static rawdata_t _try_compress(const rawdata_t& idata, const index3_t& shape, float want_snr)
  {
    if (want_snr <= 0) // or zfpy dynamic load failed?
      return rawdata_t{nullptr,0};
    rawdata_t r = _zfp_try_one_compress(idata, shape, want_snr);
    if (r.first && r.second > 0 && r.second < idata.second * 0.9)
      return r;
    r = _zfp_try_one_compress(idata, shape, 0);
    if (r.first && r.second > 0 && r.second < idata.second * 0.9)
      return r;
    return rawdata_t{nullptr,0};
  }

  /**
   * \brief Undo ZFP compression in either precision mode or lossless mode.
   *
   * The expected result size is passed in and will be checked against
   * the header of the compressed data. The expected type is "float";
   * this will also be checked.
   */
  static rawdata_t _try_decompress(const rawdata_t& cdata, const index3_t& size)
  {
    return _zfp_try_one_decompress(cdata, size);
  }
};

} // namespace

/**
 * Register this compressor in a static initializer so it will always
 * be available. Note that ~Register() needs to be empty when we use
 * the class this way.
 */
namespace {
  InternalZGY::ZfpCompressPlugin::Register static_initializer;
}

#endif // Entire file disabled if !HAVE_ZFP
