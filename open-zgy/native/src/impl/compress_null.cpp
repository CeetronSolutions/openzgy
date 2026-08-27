// Copyright 2017-2020, Schlumberger
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

#include "compression.h"
#include "../exception.h"

/**
 * \file compress_null.cpp
 * \brief Example code for the compression plug-in mechanism.
 * \details A real implementation will typically use Doxygen's copydoc
 * statements insteat of copyint the overly verbose explanations in this file.
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * \brief Example compression plug-in that always fails to compress.
 * \details This class is here just as an example and for documentation.
 * Real compression classes may use "copydoc" instead of repeating
 * what is explained here.
 *
 * Thread safety: The plug-ins and the compression algorithms they call
 * are all supposed to be thread safe. Registering the factory currently
 * is not.
 */
class NullCompressPlugin
{
public:
  /**
   * \internal This documentation will be copied to each of the
   * compression classes. So keep it general please. \endinternal
   *
   * \brief Get a functor for compressing data.
   *
   * Parse the stringly typed arguments that the factory sent us.
   * Do a consistency check on the arguments. Capture the now
   * strongly typed compressor arguments in a lambda that can be
   * used to compress one block of data.
   *
   * The returned lambda function is supposed to be thread safe.
   *
   * The return is allowed to be an empty function. The caller will
   * treat this as if no compression had been requested. This is not
   * precisely the same as returning a factor that always returns "I
   * give up". The latter case would limit the file to float32 type
   * and might turn off alignments in the file.
   *
   * If the compressor needs mutable state, e.g. to accumulate
   * perfromance measurements, it can be declared as an instance
   * method allowing the lambda to capture "this". The class would
   * normally be noncopyable. The returned lambda can be copied at
   * will; the copies will share the same state. Remember to use
   * locks to protect the shared state in "this".
   *
   * Another option for mutable state is to have this function
   * allocate some kind of "state" instance and capture it in
   * the lambda. How to access that state later, e.g. to extract
   * any timing measurements, is left as an exercise ;-).
   */
  static compressor_t getCompressor(const std::vector<std::string>&)
  {
    return [](const rawdata_t& data, const index3_t& shape) {
             return compress(data, shape);
           };
  }

  /**
   * \internal This documentation will be copied to each of the
   * compression classes. So keep it general please. \endinternal
   *
   * \brief Compress using the chosen algorithm.
   * \param data Smart pointer to data, and size in bytes of data.
   * \param shape hint for 3d layout of data in number of samples.
   * \returns Data in the same format as the input.
   *
   * \details
   *
   * The number and type of parameters, except for the first one,
   * varies depending on the algorithm being used. Some algorithms
   * might not need additional parameters at all.
   *
   * Current assumptions that apply to all compression algorithms.
   * Except for the one about thread safety these are subject to change.
   *
   * \li The method must be thread safe.
   *
   * \li The algorithm is responsible any for big / little endian
   *     conversion.
   *
   * \li The compressed data stream must never be larger than the
   *     uncompressed data. If it looks like this is going to happen
   *     then just return {nullptr, 0} telling the caller that
   *     compression is not possible.
   *
   * \li If the algorithm decides that compression is inadvisable for
   *     other reasons, such as NaN or Inf in the data, it is also
   *     acceptable to return {nullptr, 0}. Never return the input
   *     data unchanged because this will then be flagged as compressed.
   *     And the uncompress algorithm will fail miserably.
   *
   * \li The compressed data stream should start with a magic
   *     number that the decompressor can recognize.
   *
   * \internal
   * TODO-Worry I am not sure ZFP handles byte swapping correctly.
   * See the ZFP documentation. A special compilation flag is needed
   * on big endian machines. Also I suspect that the ZFP optional
   * header (which this code uses) might need explicit byte swapping.
   * \endinternal
   */
  static rawdata_t compress(const rawdata_t& data, const index3_t& shape)
  {
    return rawdata_t{nullptr, 0};
  }

  /**
   * \internal This documentation will be copied to each of the
   * compression classes. So keep it general please. \endinternal
   *
   * \brief Decompress data if it was compressed by our algorithm.
   * \param cdata Compressed data, possibly with trailing garbage.
   * \param status Might be used to identfy which algorithm was used.
   * \param shape rank and size in samples of the uncompressed result.
   *
   * The algorithm should check both the status argument and any magic
   * number at the start of the buffer to decide whether this is data
   * we know how to decompress. If not then return {nullptr, 0} which
   * should make the caller try something else.
   *
   * If the algorithm is recognized but the data cannot be decompressed
   * because it is corrupt then the method may throw an exception.
   * In that case no further algorithms will be tried.
   *
   * Older versions of this function also received the value type of the
   * data that was compressed and the value type we want it returned in
   * This is now redundant since we only allow compressing float data.
   * Anything else is simply not useful.
   *
   * Passing an uncompressed brick to this function is an error. I.e.
   * status will never be Normal. We don't have enough context to
   * handle uncompressed bricks that might require byteswapping and
   * fix for legacy quirks. Also cannot handle constant bricks,
   * missing bricks, etc.
   *
   * Current assumptions that apply to all decompression algorithms.
   * Except for the one about thread safety these are subject to change.
   *
   * \li The method must be thread safe.
   *
   * \li The compressed data stream is allowed to have trailing
   *     garbage; this must be silently ignored by the decompressor.
   *
   * \li The algorith may assume that the compressed data stream will
   *     never be longer than the uncompressed data. This is enforced
   *     by the compressor.
   *
   * \li The reason for the two assumptions above is an implementation
   *     detail. The reported size of a compressed brick isn't completely
   *     reliable. This might change in the next version of the file format.
   *
   * \li The decompressor must verify that both the status (usually
   *     set to BrickStatus.Compressed) and the magic number at the
   *     start of the compressed stream matches is correct. Relying on
   *     only status is currently not sufficient because all
   *     compressed bricks get the same status.
   *
   * \internal
   * ZFP already has a magic number. Compressors that don't have this
   * will be more awkward to support. Maybe the compressor /
   * decompressor for such types could be modified to add an extra
   * header with the compressed size and a magic number. Or we might
   * change the code to add a (size, algorithm number) header to every
   * compressed block to relieve the specific compressor /
   * decompressor from worrying about this. Or the brick status could
   * be used to encode which algorithm was used, picked up from the
   * MSB of the lup entry. Which would also require the compressor to
   * return both the actual compressed data and the code to identify
   * the decompressor. That is the main reason we are also passed the
   * "status" arg. Caveat: If adding an extra header, keep in mind
   * that this header must be included when checking that the
   * compressed stream is not larger than the input.
   * \endinternal
   */
  static rawdata_t decompress(const rawdata_t& cdata,
                              const BrickStatus& status,
                              const index3_t& shape)
  {
    return rawdata_t{nullptr, 0};
  }

  /**
   * \brief %Register the compress and decompress functions in the factory.
   *
   * Compressors use a factory so the application code can write
   * something like
   *
   * \code
   * compressor = ZgyCompressFactory("ZFP", "snr", "30");
   * \endcode
   *
   * without getting a link time dependency on the "ZFP" code.
   *
   * Decompressors use a factory because when reading blocks from the
   * file it is not possible to hard code a specific decompressor.
   *
   * It is safe to construct an instance from a static initializer.
   * Do note that this is not true for the destructor. If %Register
   * will be instanciated from a static initializer then the
   * destructor should be empty or omitted altogether.
   *
   * The reason for this is that the destruction is done in an
   * unspecified order at exit. Also, cleaning up the registry on exit
   * would be pointless anyway.
   */
  class Register
  {
  public:
    /**
     * \brief %Register this plug-in.
     *
     * It is safe to construct an instance from a static initializer.
     * Do note that this is not true for the destructor. If %Register
     * will be instanciated from a static initializer then the
     * destructor should be empty or omitted altogether.
     */
    Register()
    {
      CompressFactoryImpl::registerCompressor
        ("Null", &NullCompressPlugin::getCompressor);
      CompressFactoryImpl::registerDecompressor
        ("Null", &NullCompressPlugin::decompress);
    }
    /**
     * \brief Unregister this plug-in.
     *
     * Do NOTHING if %Register will be used as a static initializer.
     * We don't want any cleanup to happen on application exit. That
     * is pointless and it is also tricky to get destructors called
     * in the right order.
     *
     * If %Register is going to be instanciated only when needed
     * then a destructor is in order. See e.g. MockCompressPlugin
     * and MockCompressInstaller in the unit tests.
     */
    ~Register()
    {
      //CompressFactoryImpl::registerCompressor("Null", nullptr);
      //CompressFactoryImpl::registerDecompressor("Null", nullptr);
    }
  };
};

} // namespace

/**
 * %Register this compressor in a static initializer so it will always
 * be available. Note that ~Register() needs to be empty when we use
 * the class this way.
 */
namespace {
  InternalZGY::NullCompressPlugin::Register static_initializer;
}
