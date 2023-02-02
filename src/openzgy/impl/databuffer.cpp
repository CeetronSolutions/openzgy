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

#include "databuffer.h"
#include "roundandclip.h"
#include "exception.h"
#include "environment.h"
#include "fancy_timers.h"
#include "minmaxscan.h"

#include <cstdint>
#include <array>
#include <algorithm>
#include <vector>
#include <memory.h>
#include <typeinfo>
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <atomic>

#define XXX_ENABLE_WARMUP_LOGGING 0 /* beware, will slow down execution */

namespace InternalZGY {
#if 0
}
#endif

namespace {
#if 0
}
#endif

namespace {
  /**
   * \brief Instrumenting DataBufferNd for performance measurements.
   *
   * The following lists functions that are public, potentially
   * expensive, and are not called from other members of DataBuffer
   * which might cause us to count the same task twice. Beware that
   * I did not check this list very thoroughly.
   *
   * Declaring all timers in one spot is not really needed, at least
   * not when using static instances that report when the app exits.
   * If I need a report each time a DataBuffer goes out of scope
   * it would be different.
   *
   * For slightly more advanced usage, a global cleartimers() can
   * be called to output all collected details at that point.
   *
   * \details Thread safety: Safe because SummaryPrintingTimerEx is
   * thread safe when used correctly.
   */
  class AdHocTimers
  {
  public:
    SummaryPrintingTimerEx ctor;
    SummaryPrintingTimerEx allsame;
    SummaryPrintingTimerEx fill;
    SummaryPrintingTimerEx range;
    SummaryPrintingTimerEx clone;
    SummaryPrintingTimerEx scaletof;
    SummaryPrintingTimerEx scaletos;
    SummaryPrintingTimerEx copyscalar;
    SummaryPrintingTimerEx copysubset;
#if XXX_ENABLE_WARMUP_LOGGING
    SummaryPrintingTimerEx warm;
#endif
    AdHocTimers()
      : ctor      ("DBNd.ctor")
      , allsame   ("DBNd.allsame")
      , fill      ("DBNd.fill")
      , range     ("DBNd.range")
      , clone     ("DBNd.clone")
      , scaletof  ("DBNd.scaletof")
      , scaletos  ("DBNd.scaletos")
      , copyscalar("DBNd.copyscalar")
      , copysubset("DBNd.copysubset")
#if XXX_ENABLE_WARMUP_LOGGING
      , warm      ("DBNd.warm")
#endif
    {}
    static AdHocTimers& instance()
    {
      static AdHocTimers instance_;
      return instance_;
    }
    void cleartimers(bool show) {
      if (show) {
#if XXX_ENABLE_WARMUP_LOGGING
        warm.print();
#endif
        copysubset.print();
        copyscalar.print();
        scaletos.print();
        scaletof.print();
        clone.print();
        range.print();
        fill.print();
        allsame.print();
        ctor.print();
      }
      else {
#if XXX_ENABLE_WARMUP_LOGGING
        warm.reset();
#endif
        copysubset.reset();
        copyscalar.reset();
        scaletos.reset();
        scaletof.reset();
        clone.reset();
        range.reset();
        fill.reset();
        allsame.reset();
        ctor.reset();
      }
    }
  };
}

#if 0
/**
 * The default for OPENZGY_COPYSUBSET_SHORTCUT is 0 i.e. off, because
 * it hasn't been shown that it gives a measurable speedup and it
 * hasn't been sufficiantly tested.
 */
static int
copysubset_shortcut()
{
  static int enable = Environment::getNumericEnv("OPENZGY_COPYSUBSET_SHORTCUT", 0);
  return enable;
}
#endif

/**
 * Copy subset of one array into a subset of another.
 *
 * In typical usage, either source or destination will be a single
 * brick read from or being written to the file. "orig" is the
 * survey-relative position corresponding to the first sample in the
 * memory buffer. The other buffer is typically something the user
 * supplied. With "orig" being the first sample the user wants to read
 * or write.
 *
 * @param ndim      Highest dimension to consider.
 *
 * @param srcorig   Survey coordinates of the first sample in the destination.
 * @param srcsize   Size of the source array.
 * @param srcstride The stride to use, per dimension, for the source array.
 * @param srcbuf    Pointer to buffer holding source array.
 *
 * @param dstorig   Survey coordinates of the first sample in the destination.
 * @param dstsize   Size of the destination array.
 * @param dststride The stride to use, per dimension, for the destination array.
 * @param dstbuf    Pointer to buffer holding destination array.
 *
 * @param cpyorig   if non-null: Don't copy data before this, even if possible.
 * @param cpysize   If non-null: Don't copy past cpyorig + cpysize, even if
 *                  possible. If setting either of cpyorig or copysize then
 *                  both should be set.
 *
 * Based on CopySubset() in ArrayConvert.cpp in the old accessor.
 * This can also replace Permute(), as that one is
 * apparently just a special case with
 *
 *   srcorig = dstorig = cpyorig = (0,0,0);
 *   srcsize = dstsize = cpysize = size;
 *
 * As for the other code in ArrayConvert: RescaleConvert may be too general
 * for my use because OpenZGY doesn't allow conversion between arbitrary
 * representations such as between int8/int16 or int16/uint16. Big-endian
 * architectures are currently not supported so this makes many of the other
 * ArrayConvert methods moot as well.
 *
 * CopySubset needs to be templated on item type only in order to efficiently
 * handle the case where the fastest varying dimension isn't 1 in both source
 * and target. If this was not a requirement I could have used void* and an
 * explicit elementsize instead.
 *
 * TODO-Performance: The code is supposed to handle any buffer strides,
 * but is a LOT more efficient if the last dimension has step 1. This is
 * normally the case in OpenZGY but somebody might want to re-use it in
 * pther contexts.
 */
template <typename T>
static void
CopySubset(std::int32_t ndim,
           const std::int64_t* srcorig,
           const std::int64_t* srcsize,
           const std::int64_t* srcstride,
           const T* srcbuf,
           const std::int64_t* dstorig,
           const std::int64_t* dstsize,
           const std::int64_t* dststride,
           T* dstbuf,
           const std::int64_t* cpyorig = 0,
           const std::int64_t* cpysize = 0)
{
  // Sanity check.
  if (ndim < 1)
    return;

  if (cpyorig == nullptr || cpysize == nullptr) {
    // This basically ignores the "copy" limits.
    cpyorig = srcorig;
    cpysize = srcsize;
  }

#if 0
  static std::atomic<std::int64_t> totalcall{0}, fastcall{0};
  ++totalcall;
  if (ndim==3 && copysubset_shortcut() > 1) {
    static auto three= [](const std::int64_t *a) -> std::string {
                         if (!a)
                           return std::string("(nullptr)");
                         std::stringstream ss;
                         ss << "(" << a[0] << "," << a[1] << "," << a[2] << ")";
                         return ss.str();
                       };
    if (totalcall > 0 && (totalcall % 1000) == 0 && copysubset_shortcut() > 1) {
      std::cerr << "Fast CopySubset " << fastcall << "/" << totalcall
                << " " << three(srcsize) << " * " << sizeof(T)
                << " lockfree " << totalcall.is_lock_free()
                << std::endl;
    }
    if (copysubset_shortcut() > 2)
      std::cerr << "CopySubset shortcut " << copysubset_shortcut()
                << " srcorig "   << three(srcorig)
                << " srcsize "   << three(srcsize)
                << " srcstride " << three(srcstride)
                << " dstorig "   << three(dstorig)
                << " dstsize "   << three(dstsize)
                << " dststride " << three(dststride)
                << " cpyorig "   << three(cpyorig)
                << " cpysize "   << three(cpysize)
                << std::endl;
  }

  // Do a short cut if the entire copy operation can be made
  // using a single memcpy. The test is very specific and there
  // are other cases where a single memcpy would have worked.
  // But inside OpenZGY I believe this is the only one that
  // matters. The code will be hit when reading or writing a
  // single brick. NOTE that it might be better to test for that
  // case at a higher level. In which case the test here can
  // be removed. Initially the test needs to be here to debug
  // a particular issue.
  if (ndim == 3 &&
      copysubset_shortcut() &&
      // Source is C-Contiguous?
      srcstride[2] == 1 &&
      srcstride[1] == srcsize[0] &&
      srcstride[0] == srcsize[0] * srcsize[1] &&
      // Neither source nor destination has any offset?
      srcorig[0] == 0 &&
      srcorig[1] == 0 &&
      srcorig[2] == 0 &&
      dstorig[0] == 0 &&
      dstorig[1] == 0 &&
      dstorig[2] == 0 &&
      // Same size?
      srcsize[0] == dstsize[0] &&
      srcsize[1] == dstsize[1] &&
      srcsize[2] == dstsize[2] &&
      // Same stride?
      srcstride[1] == dststride[1] &&
      srcstride[2] == dststride[2] &&
      srcstride[0] == dststride[0] &&
      // At this point, target must also be C-Contiguous.
      // Not clipped by survey?
      cpyorig[0] <= srcorig[0] && cpysize[0] >= srcsize[0] + cpyorig[0] &&
      cpyorig[1] <= srcorig[1] && cpysize[1] >= srcsize[1] + cpyorig[1] &&
      cpyorig[2] <= srcorig[2] && cpysize[2] >= srcsize[2] + cpyorig[2])
  {
    ++fastcall;
    memcpy(dstbuf, srcbuf, dstsize[0] * dstsize[1] * dstsize[2] * sizeof(T));
    return;
  }
#endif

  // Calculate overlap range.
  const std::int64_t beg = std::max(std::max(dstorig[0], srcorig[0]), cpyorig[0]);
  const std::int64_t end = std::min(std::min(dstorig[0] + dstsize[0], srcorig[0] + srcsize[0]), cpyorig[0] + cpysize[0]);

  if (beg >= end) {
    // degenerate case: empty overlap. Nothing to do.
  }
  else if (ndim == 1) {
    // Degenerate case: 1-dimensional.
    if (srcstride[0] == 1 && dststride[0] == 1) {
      memcpy(dstbuf + (beg - dstorig[0]),
             srcbuf + (beg - srcorig[0]),
             (end - beg)*sizeof(T));
    }
    else {
      // TODO-Test remember to test this case with a nontrivial region.
      // Copy one element at a time. Using memcpy is not as expensive
      // as it looks. As long as the compiler can see the constant
      // size it is free to replace memcpy with an assignment.
      // memcpy() is more robust wrt. misaligned data. Also, we probably
      // won't get here very often because the stride of the last dim is
      // normally 1.
      T* dst       = &dstbuf[(beg - dstorig[0])*dststride[0]];
      const T* src = &srcbuf[(beg - srcorig[0])*srcstride[0]];
      for (std::int64_t ii = beg; ii < end; ++ii,
             dst += dststride[0], src += srcstride[0])
        memcpy(dst, src, sizeof(T));
    }
  }
  else {
    // general case: recursive
    for (std::int64_t ii = beg; ii < end; ++ii) {
      CopySubset(ndim - 1,
                 srcorig + 1, srcsize + 1, srcstride + 1, srcbuf + (ii - srcorig[0])*srcstride[0],
                 dstorig + 1, dstsize + 1, dststride + 1, dstbuf + (ii - dstorig[0])*dststride[0],
                 cpyorig + 1, cpysize + 1);
    }
  }
}

/**
 * Copy a constant value into a subset of an array.
 * Apart from the source being a constant this is identical to CopySubset.
 */
template <typename T>
static void
CopyScalar(std::int32_t ndim,
           const std::int64_t* srcorig,
           const std::int64_t* srcsize,
           const T scalar,
           const std::int64_t* dstorig,
           const std::int64_t* dstsize,
           const std::int64_t* dststride,
           T* dstbuf,
           const std::int64_t* cpyorig = 0,
           const std::int64_t* cpysize = 0)
{
  // Sanity check.
  if (ndim < 1)
    return;

  if (cpyorig == nullptr || cpysize == nullptr) {
    // This basically ignores the "copy" limits.
    cpyorig = dstorig;
    cpysize = dstsize;
  }

  // Calculate overlap range.
  const std::int64_t beg = std::max(std::max(dstorig[0], srcorig[0]), cpyorig[0]);
  const std::int64_t end = std::min(std::min(dstorig[0] + dstsize[0], srcorig[0] + srcsize[0]), cpyorig[0] + cpysize[0]);

  if (beg >= end) {
    // degenerate case: empty overlap. Nothing to do.
  }
  else if (ndim == 1) {
    std::int64_t step = dststride[0];
    T* it_dst = &dstbuf[(beg - dstorig[0]) * step];
    T* it_end = &dstbuf[(end - dstorig[0]) * step];
    if (step == 1) {
      std::fill(it_dst, it_end, scalar);
    }
    else {
      for (; it_dst < it_end; it_dst += step)
          *it_dst = scalar;
    }
  }
  else {
    // general case: recursive
    for (std::int64_t ii = beg; ii < end; ++ii) {
      CopyScalar(ndim - 1,
                 srcorig + 1, srcsize + 1, scalar,
                 dstorig + 1, dstsize + 1, dststride + 1, dstbuf + (ii - dstorig[0])*dststride[0],
                 cpyorig + 1, cpysize + 1);
    }
  }
}

/**
 * Convert a voidptr + size to a DataBuffer. See makeDataBuffer3d.
 *
 * To avoid too much copy/paste, this private method is able to create
 * several different buffer types by calling different DataBufferNd
 * constructors. The code to convert from a RawDataType enum into a
 * C++ template is messy. I don't want to repeat that that for all 4
 * cases.
 *
 * If voidptr is null this creates an uninitialized 3d DataBuffer with
 * the given size and type and isScalar() == false. nbytes should be 0.
 *
 * If voidptr is not null then this creates a scalar DataBuffer if
 * nbytes indicates there is room for just one sample. If there are
 * enough bytes it creates a regular DataBuffer pointing to the voidptr
 * that was passed in.
 *
 * The valid values for nbytes are:
 *
 *    0              - If and only if voidptr is null. Allocate new data.
 *    sizeof(T)      - (T*)voidptr contains the scalar to use.
 *    sizeof(double) - (double*)voidptr contains the scalar to use,
 *    size[0]*size[1]*size[2]*sizeof(T) - (T*)voidptr is used as-is.
 *
 * Any other value will raise an exception.
 *
 * The function needs to be explicitly told whether the caller wants a
 * scalar buffer or not. Otherwise there is an ambiguity when
 * size0*size1*size2 is either 1 or sizeof(double)/sizeof(T).
 */
template<typename T, int NDim>
static std::shared_ptr<DataBuffer>
_makeDataBufferT(const std::shared_ptr<void>& raw, std::int64_t rawsize,
                 const std::array<std::int64_t,NDim> bricksize, bool make_scalar)
{
  const std::int64_t nsamples = bricksize[0] * bricksize[1] * bricksize[2];
  if ((reinterpret_cast<std::intptr_t>(raw.get()) % sizeof(T)) != 0)
    throw OpenZGY::Errors::ZgyInternalError("Misaligned buffer.");
  std::shared_ptr<T> rawdata = std::static_pointer_cast<T>(raw);
  if (make_scalar) {
    if (rawdata && rawsize == sizeof(T)) {
      T value;
      memcpy(&value, raw.get(), sizeof(T));
      return std::make_shared<DataBufferNd<T,NDim>>(value, bricksize);
    }
    else if (rawdata && rawsize == (std::int64_t)sizeof(double)) {
      double dblvalue;
      memcpy(&dblvalue, rawdata.get(), sizeof(double));
      T value = static_cast<T>(dblvalue);
      return std::make_shared<DataBufferNd<T,NDim>>(value, bricksize);
    }
    else {
      throw OpenZGY::Errors::ZgyInternalError("Creating a scalar buffer with wrong input size.");
    }
  }
  else if (!rawdata && rawsize == 0) {
    return std::shared_ptr<DataBufferNd<T,NDim>>
      (new DataBufferNd<T,NDim>(bricksize));
  }
  else if (rawsize == (std::int64_t)sizeof(T) * nsamples) {
    return std::shared_ptr<DataBufferNd<T,NDim>>
        (new DataBufferNd<T,NDim>(rawdata, bricksize));
  }
  else {
    std::cerr << "@@ scalar? " << make_scalar
              << " rawdata " << (rawdata ? "ok" : "nullptr")
              << " rawsize " << rawsize
              << " bricksize ("
              << bricksize[0] << ","
              << bricksize[1] << ","
              << bricksize[2] << ")"
              << " vtsize " << sizeof(T)
              << std::endl;
    throw OpenZGY::Errors::ZgyInternalError("Invalid arguments to makeDataBuffer.");
  }
}

std::shared_ptr<DataBuffer>
_makeDataBuffer3d(
     const std::shared_ptr<void>& raw, std::int64_t nbytes,
     const std::array<std::int64_t,3>& size,
     RawDataType dtype, bool make_scalar)
{
  switch (dtype) {
  case RawDataType::SignedInt8:    return _makeDataBufferT<std::int8_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::UnsignedInt8:  return _makeDataBufferT<std::uint8_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::SignedInt16:   return _makeDataBufferT<std::int16_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::UnsignedInt16: return _makeDataBufferT<std::uint16_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::SignedInt32:   return _makeDataBufferT<std::int32_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::UnsignedInt32: return _makeDataBufferT<std::uint32_t,3>(raw, nbytes, size, make_scalar);
  case RawDataType::Float32:       return _makeDataBufferT<float,3>(raw, nbytes, size, make_scalar);
  case RawDataType::IbmFloat32:
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

} // end anonymous namespace

DataBuffer::~DataBuffer()
{
}

/**
 * Convert voidptr + nbytes + size3d to a DataBuffer.
 *
 * The method is allowed to make a copy of the data, so the caller
 * must not assume that the pointer it passed to the constructor
 * can still be used to access the buffer contents.
 *
 * The input data needs to be correctly aligned for its type. This is
 * checked in _makeDataBufferT(). If it turns out that unaligned input
 * should be supported then it is trivial to just make a new buffer
 * and copy the data to it.
 *
 * What is less trivial is to leave the data misaligned and assume
 * that all consumers such as CopySubset() will be using memcpy()
 * which doesn't have problems with alignment. It might actually work
 * but the risk of introducing odd bugs is too high.
 */
std::shared_ptr<DataBuffer>
DataBuffer::makeDataBuffer3d(
     const std::shared_ptr<void>& raw, std::int64_t nbytes,
     const std::array<std::int64_t,3>& size,
     RawDataType dtype)
{
  // Hide the feature that a nullptr allocates a new buffer.
  // Force callers to explicitly use makeNewBuffer3d for that.
  if (!raw || nbytes == 0)
    throw OpenZGY::Errors::ZgyInternalError("makeDataBuffer3d called without data.");
  std::shared_ptr<DataBuffer> result =
    _makeDataBuffer3d(raw, nbytes, size, dtype, /*make_scalar=*/false);
  if (result->datatype() != dtype || result->size3d() != size)
    throw OpenZGY::Errors::ZgyInternalError("makeDataBuffer3d made a mistake.");
  return result;
}

/**
 * \brief Create a DataBuffer with a newly allocated, not intialized buffer.
 *
 * \details
 * Caller is responsible for either initializing the buffer to a known
 * default value or making sure the entire buffer has valid data
 * written to it.
 */
std::shared_ptr<DataBuffer>
DataBuffer::makeNewBuffer3d(
     const std::array<std::int64_t,3>& size,
     RawDataType dtype)
{
  std::shared_ptr<DataBuffer> result =
    _makeDataBuffer3d(nullptr, 0, size, dtype, /*make_scalar=*/false);
  if (result->datatype() != dtype || result->size3d() != size)
    throw OpenZGY::Errors::ZgyInternalError("makeNewBuffer3d made a mistake.");
  return result;
}

/**
 * \brief Create a DataBuffer representing a scalar, passed in as a double.
 *
 * \details
 * If you want to pass the scalar value as the buffer's actual data
 * type, as an array length 1 if T, an overload can be created for
 * that but currently there is no need for it.
 */
std::shared_ptr<DataBuffer>
DataBuffer::makeScalarBuffer3d(
     double value,
     const std::array<std::int64_t,3>& size,
     RawDataType dtype)
{
  std::shared_ptr<DataBuffer> result =
    _makeDataBuffer3d(std::make_shared<double>(value), sizeof(double), size, dtype, /*make_scalar=*/true);
  if (result->datatype() != dtype || result->size3d() != size)
    throw OpenZGY::Errors::ZgyInternalError("makeScalarBuffer3d made a mistake.");
  return result;
}

template <typename T, int NDim>
DataBufferNd<T,NDim>::DataBufferNd(T scalar, const std::array<std::int64_t,NDim>& size)
  : _size(size)
  , _stride({0})
  , _is_scalar(true)
  , _data(std::shared_ptr<T>(new T[1]{scalar}, std::default_delete<T[]>()))
{
  for (int dim=0; dim<NDim; ++dim)
    if (size[dim] < 1)
      throw OpenZGY::Errors::ZgyInternalError("DataBuffer size must be positive");
}

template <typename T, int NDim>
DataBufferNd<T,NDim>::DataBufferNd(const std::array<std::int64_t,NDim>& size)
  : _size(size)
  , _stride(cstride(size))
  , _is_scalar(false)
  , _data(std::shared_ptr<T>(new T[make_product(size)], std::default_delete<T[]>()))
{
  for (int dim=0; dim<NDim; ++dim)
    if (size[dim] < 1)
      throw OpenZGY::Errors::ZgyInternalError("DataBuffer size must be positive");
  // TODO-Worry: Leaving the buffer not initialized might turn an easily
  // reproducible bug into a heissenbug. But zeroing the buffers here
  // is too wasteful. If you must, then do a fill() just after constructing
  // a new buffer,
  //memset(_data.get(), 0, make_product(size) * sizeof(T));
}

template <typename T, int NDim>
DataBufferNd<T,NDim>::DataBufferNd(const std::shared_ptr<T>& data, const std::array<std::int64_t,NDim>& size)
  : _size(size)
  , _stride(cstride(size))
  , _is_scalar(false)
  , _data(data)
{
  for (int dim=0; dim<NDim; ++dim)
    if (size[dim] < 1)
      throw OpenZGY::Errors::ZgyInternalError("DataBuffer size must be positive");
  if (!data)
    throw OpenZGY::Errors::ZgyInternalError("DataBuffer storage cannot be null");
}

/**
 * Set up a default stride with the last dimension varying fastest
 * and with no holes in the data.
 */
template <typename T, int NDim>
typename DataBufferNd<T,NDim>::ndsize_t
DataBufferNd<T,NDim>::cstride(const ndsize_t& size)
{
  ndsize_t result;
  result[NDim-1] = 1;
  for (int dim = NDim-2; dim >= 0; --dim)
    result[dim] = result[dim+1] * size[dim+1];
  return result;
}

/**
 * Convert a raw pointer to a fixed-length array. There is an unsafe
 * assumption thay the array passed by the caller contains the correct
 * number of elements.
 */
template <typename T, int NDim>
typename DataBufferNd<T,NDim>::ndsize_t
DataBufferNd<T,NDim>::make_array(const std::int64_t *in)
{
  ndsize_t result;
  for (int ii=0; ii<NDim; ++ii)
    result[ii] = in[ii];
  return result;
}

template <typename T, int NDim>
std::int64_t
DataBufferNd<T,NDim>::make_product(const ndsize_t& size)
{
  std::int64_t result = 1;
  for (int dim = 0; dim < NDim; ++dim)
    result *= size[dim];
  return result;
}

template <typename T, int NDim>
std::string
DataBufferNd<T,NDim>::toString() const
{
  std::stringstream ss;
  ss << "type " << typeid(T).name() << " size (";
  for (int dim = 0; dim < NDim; ++dim)
    ss << _size[dim] << (dim == NDim-1 ? ")" : ", ");
  ss << " stride (";
  for (int dim = 0; dim < NDim; ++dim)
    ss << _stride[dim] << (dim == NDim-1 ? ")" : ", ");
  if (isScalar())
    ss << " constant " << (double)scalarValue();
  else
    ss << " ptr 0x" << std::hex << (intptr_t)voidData().get() << std::dec;
  return ss.str();
}

template <typename T, int NDim>
std::int64_t
DataBufferNd<T,NDim>::allocsize() const
{
  return isScalar() ? 1 : totalsize();
}

template <typename T, int NDim>
std::int64_t
DataBufferNd<T,NDim>::totalsize() const
{
  return make_product(_size);
}

template <typename T, int NDim>
std::int64_t
DataBufferNd<T,NDim>::itemsize() const
{
  return sizeof(T);
}

/**
 * Return the size of the buffer, assuming it is 3d. Which it almost always is.
 * As a code smell this suggests that DataBuffer maybe shouldn't have been
 * templated on NDim in the first place. Or have NDim as a regular varable.
 * If fewer than 3 dimensions exist then the first index or indices are
 * treated as size=1, stride=0. If more then 3 dimensions only the last 3
 * will be returned. Which is probably not very useful.
 */
template <typename T, int NDim>
std::array<std::int64_t,3>
DataBufferNd<T,NDim>::size3d() const
{
  std::array<std::int64_t,3> result;
  result[2] = (NDim >= 1) ? _size[NDim-1] : 0;
  result[1] = (NDim >= 2) ? _size[NDim-2] : 1;
  result[0] = (NDim >= 3) ? _size[NDim-3] : 1;
  return result;
}

/**
 * Return the layout of the buffer, assuming it is 3d. See size3d().
 * Layout is always c-contiguous so this can be trivially computed
 * from the size.
 */
template <typename T, int NDim>
std::array<std::int64_t,3>
DataBufferNd<T,NDim>::stride3d() const
{
  std::array<std::int64_t,3> result;
  result[2] = (NDim >= 1) ? _stride[NDim-1] : 0;
  result[1] = (NDim >= 2) ? _stride[NDim-2] : 0;
  result[0] = (NDim >= 3) ? _stride[NDim-3] : 0;
  return result;
}

/**
 * \brief Set all samples to the same value.
 *
 * The value is passed as a float and will be cast to the correct type.
 * There is NO check for overflow.
 *
 * Some optimizing is attempted because gprof claims that (a) calling
 * this function amounts to 50% of the read overhead, and (b) memset,
 * if possible, would reduce that time significantly. I suspect that
 * both of those numbers are false. An ad-hoc Timer reports 10% in both
 * cases. But just in case I am wrong I have added the optimization.
 *
 * Using memset is possible if the value to be set is zero (the common
 * case for float cubes) and if the sample size is 1 (all int8 cubes).
 * Int16 cubes are less likely to benefit.
 */
template <typename T, int NDim>
void
DataBufferNd<T,NDim>::fill(double value)
{
  const T val = static_cast<T>(value);
#if XXX_ENABLE_WARMUP_LOGGING
  if (!isScalar()) {
    SimpleTimerEx ww(AdHocTimers::instance().warm);
    if (value == 0)
      memset(_data.get(), 0, allocsize() * sizeof(T));
    else
      std::fill(_data.get(), _data.get() + allocsize(), val);
    AdHocTimers::instance().warm.addBytesWritten(allocsize() * sizeof(T));
  }
#endif
  SimpleTimerEx tt(AdHocTimers::instance().fill);
  if (isScalar())
    _data.get()[0] = val;
  else {
    if (true && value == 0)
      memset(_data.get(), 0, allocsize() * sizeof(T));
    else if (false && sizeof(T) == 1)
      memset(_data.get(), static_cast<int>(val), allocsize() * sizeof(T));
    else
      std::fill(_data.get(), _data.get() + allocsize(), val);
  }
}

/**
 * Make the buffer empty. Size and stride are unchanged but the data will be
 * nullptr. This will cause an access violation if trying to read or write
 * data. The reason this is useful is that _data might point to something
 * that is not reference counted by our smart pointer. The pointer's deleter
 * is then a no-op. clear() should be called if we can no longer trust
 * that the pointer is valid. A reproducible null pointer exception is way
 * better than random crashes.
 */
template <typename T, int NDim>
void
DataBufferNd<T,NDim>::clear()
{
  _data.reset();
}

/**
 * \brief Find lowest and highest sample excluding +/- Inf and NaN.
 *
 * \details The OpemMP directive is commented out because it might
 * unfortunately turn out to run slower. The inner loop is very
 * inexpensive which means the code is sensitive to how much overhead
 * OpenMP will add. I have made real-world measurements but only in
 * one hardware / OS / compiler environment and only for float input.
 * I observed that 1 million elements was the point where
 * multi-threading barely started to have a benefit. In a typical
 * scenario I don't expect to see more then 5 times that number. That
 * is a typical size of one brick-column. With 5 million floats the
 * task ran 3.5 times faster using 16 threads. With a less capable
 * OpenMP installation it might easily have run slower. It could
 * potentially run a *lot* slower. In my test, range() accounted for
 * 5% of total elapsed time. So I am, weighing a potential 5%
 * performance improvement against the risk of slowing down. I am
 * chicken. I hope I can instead parallelize this at a higher level.
 * Making this issue moot.
 */
template <typename T, int NDim>
std::pair<double, double>
DataBufferNd<T,NDim>::range(const std::int64_t *used_in) const
{
  SimpleTimerEx tt(AdHocTimers::instance().range);
  if (used_in)
    for (int dim=0; dim<NDim; ++dim)
      if (used_in[dim] > _size[dim])
        throw OpenZGY::Errors::ZgyInternalError("Used size > buffer");
  std::array<std::int64_t,NDim> used;
  for (int dim=0; dim<NDim; ++dim)
    used[dim] = used_in ? used_in[dim] : _size[dim];
  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::lowest();
  if (isScalar()) {
    if (IsFiniteT(scalarValue()))
      min = max = scalarValue();
  }
  else if (used == _size) {
    // Buffer is always c-contiguous. Treat as 1d.
    const T * const ptr = data();
    const std::int64_t totalsize = allocsize();
    if (std::is_same<T, float>::value) {
        auto float_ptr = reinterpret_cast<const float * const>(data());
        float fMin, fMax;
        MinMaxScan::scanArray(float_ptr, allocsize(), 1, fMin, fMax);
        // The cast is just to keep the compiler happy; we already know that T is float.
        min = (T)fMin;
        max = (T)fMax;
    }
    else {
        //#pragma omp parallel for reduction(min: min) reduction(max: max) if(totalsize >= 1*1024*1024)
        for (std::int64_t ii = 0; ii < totalsize; ++ii) {
            const T value = ptr[ii];
            if (!IsFiniteT(value)) continue;
            if (min > value) min = value;
            if (max < value) max = value;
        }
    }
  }
  else {
    std::array<std::int64_t,NDim> loop{0};
    while (loop[0] < used[0]) {
      const T* ptr = data();
      for (int dim=0; dim<NDim; ++dim) {
        ptr += loop[dim] * _stride[dim];
      }
      // Roll out the inner loop; also assume last stride is 1.
      for (const T* end = ptr + used[NDim-1]; ptr < end; ++ptr) {
        const T value = *ptr;
        if (!IsFiniteT(value)) continue;
        if (min > value) min = value;
        if (max < value) max = value;
      }
      // Increment iterator, but not final dim since that is rolled out.
      // The slowest dim falls off the end so the loop can terminate;
      // the intermediate dims will wrap and increment the next slower dim.
      for (int dim = NDim-2; dim >= 0; --dim) {
        if (++loop[dim] < used[dim])
          break;
        if (dim != 0)
          loop[dim] = 0;
      }
    }
  }
  return std::make_pair(static_cast<double>(min), static_cast<double>(max));
}

template <typename T, int NDim>
bool
DataBufferNd<T,NDim>::isAllSame(const std::int64_t *used_in) const
{
  if (isScalar())
    return true;
  SimpleTimerEx tt(AdHocTimers::instance().allsame);
  // TODO-Low yet another place I wish I hadn't templeted on NDim.
  std::array<std::int64_t,NDim> used;
  for (int dim=0; dim<NDim; ++dim)
    used[dim] = used_in ? used_in[dim] : _size[dim];

  if (used[0]<=0 || used[1] <= 0 || used[2] <= 0 || isScalar()) {
    // Has zero or one sample, so clearly all samples are the same.
    return true;
  }
  else if (used[0] > _size[0] || used[1] > _size[1] || used[2] > _size[2]) {
    // Caller made an error. Should maybe assert.
    return false;
  }
  else if (used == _size) {
    // Short cut: No holes and no unused area. Treat as 1d.
    const T *ptr = data();
    const T *end = ptr + allocsize();
    const T first = *ptr++;
    if (IsNanT(first)) {
      for (; ptr < end; ++ptr)
        if (!IsNanT(*ptr))
          return false;
    }
    else {
      for (; ptr < end; ++ptr)
        if (*ptr != first)
          return false;
    }
  }
  else {
    const T *ptr = data();
    const T first = *ptr;
    if (IsNanT(first)) {
      for (int ii=0; ii<used[0]; ++ii)
        for (int jj=0; jj<used[1]; ++jj)
          for (int kk=0; kk<used[2]; ++kk)
            if (!IsNanT(ptr[ii*_stride[0]+jj*_stride[1]+kk*_stride[2]]))
              return false;
    }
    else{
      for (int ii=0; ii<used[0]; ++ii)
        for (int jj=0; jj<used[1]; ++jj)
          for (int kk=0; kk<used[2]; ++kk)
            if (ptr[ii*_stride[0]+jj*_stride[1]+kk*_stride[2]] != first)
              return false;
    }
  }
  return true;
}

template <typename T, int NDim>
void
DataBufferNd<T,NDim>::copySubset(const ndsize_t& srcorig, const self_type& src,
                                 const ndsize_t& dstorig,       self_type& dst,
                                 const ndsize_t& cpyorig, const ndsize_t& cpysize)
{
  if (dst.isScalar())
    throw OpenZGY::Errors::ZgyInternalError("Attempted copy to scalar buffer.");

  if (src.isScalar())
    CopyScalar(NDim,
               srcorig.data(), src._size.data(), src.scalarValue(),
               dstorig.data(), dst._size.data(), dst._stride.data(), dst.data(),
               cpyorig.data(), cpysize.data());
  else
    CopySubset(NDim,
               srcorig.data(), src._size.data(), src._stride.data(), src.data(),
               dstorig.data(), dst._size.data(), dst._stride.data(), dst.data(),
               cpyorig.data(), cpysize.data());
}

template <typename T, int NDim>
void
DataBufferNd<T,NDim>::copySubset(const ndsize_t& srcorig, const self_type& src,
                                 const ndsize_t& dstorig,       self_type& dst)
{
  if (dst.isScalar())
    throw OpenZGY::Errors::ZgyInternalError("Attempted copy to scalar buffer.");

  if (src.isScalar())
    CopyScalar(NDim,
               srcorig.data(), src._size.data(), src.scalarValue(),
               dstorig.data(), dst._size.data(), dst._stride.data(), dst.data(),
               nullptr, nullptr);
  else
    CopySubset(NDim,
               srcorig.data(), src._size.data(), src._stride.data(), src.data(),
               dstorig.data(), dst._size.data(), dst._stride.data(), dst.data(),
               nullptr, nullptr);
}

template <typename T, int NDim>
void
DataBufferNd<T,NDim>::copyFrom(const DataBuffer* src,
                               const std::int64_t *srcorig,
                               const std::int64_t *dstorig,
                               const std::int64_t *cpyorig,
                               const std::int64_t *cpysize)
{
  SimpleTimerEx tt(src && src->isScalar() ?
                   AdHocTimers::instance().copyscalar :
                   AdHocTimers::instance().copysubset);

  const self_type* source = dynamic_cast<const self_type*>(src);
  if (!source) {
    if (src)
      throw OpenZGY::Errors::ZgyInternalError("Attempted copy from wrong buffer type.");
    else
      throw OpenZGY::Errors::ZgyInternalError("Attempted copy from null buffer.");
  }
  if (cpyorig && cpysize)
    copySubset(make_array(srcorig), *source,
               make_array(dstorig), *this,
               make_array(cpyorig), make_array(cpysize));
  else
    copySubset(make_array(srcorig), *source,
               make_array(dstorig), *this);
}

/**
 * Copy the input DataBuffer into a newly allocated buffer of float.
 * The input is assumed to have been read from a ZGY file,
 * so this function will also apply the storage to float transform.
 *
 * If the buffer is already floating point then it needs neither
 * type conversion nor scaling and the function will return an
 * empty pointer. Note that it does *not* return its input argument,
 * because the input is a plain pointer and I don't want to mess with
 * std::enable_shared_from_this.
 *
 * Note that both s_scaleToFloat and s_scaleFromFloat are private static
 * methods that are templated on the integral type. This means that
 * s_scaleToFloat is templated on the input type and s_scaleFromFloat
 * is templated on the desired target type. The public interface is
 * in the virtual scaleToFloat() and scaleToStorage() methods.
 *
 * If the buffer is non contiguous the contents of the padding area is
 * unspecified. It might be garbage or values converted from the input.
 * This violates the principle of least surprise but is likely harmless.
 *
 * TODO-Low: Should I try to combine the type conversion and scaling
 * done here with the functionality of copySubset? This might allow
 * skipping one temporary buffer but I don't know how easy it will be
 * for the bulk access to make use of it.
 */
template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::s_scaleToFloat(const DataBuffer* in,
                                     const std::array<double,2>& factors)
{
  typedef self_type src_type;
  typedef DataBufferNd<float,NDim> dst_type;
  // Factors are used as-is when converting from int storage to float user.
  const double a = factors[0]; // slope
  const double b = factors[1]; // intercept
  const self_type *src = dynamic_cast<const src_type*>(in);
  if (src == nullptr) {
    if (in == nullptr)
      return std::shared_ptr<dst_type>();
    else
      throw OpenZGY::Errors::ZgyInternalError("Mixing types in DataBufferNd::scaleToFloat");
  }
  else if (!std::numeric_limits<T>::is_integer) {
    return std::shared_ptr<dst_type>();
  }
  else if (src->isScalar()) {
      return std::shared_ptr<dst_type>(
          new dst_type(static_cast<float>(src->scalarValue() * a + b), src->safesize()));
  }
  else {
    auto dst = std::shared_ptr<dst_type>(new dst_type(src->safesize()));
    // TODO-Worry: Beware of edge bricks.
    // It might not be a good idea to convert the padding area.
    // Note: dst_type::value_type is always float, but just in
    // case this changes e.g. to double I will spell it out.
    // src_type::value_type is the same as T.
    typename dst_type::value_type *dst_ptr = dst->data();
    const src_type::value_type *src_ptr = src->data();
    const src_type::value_type *src_end = src_ptr + src->allocsize();
#if 0
    std::cerr << "@@ Converting" << std::setprecision(14)
              << " " << *src_ptr
              << " to " << static_cast<typename dst_type::value_type>(*src_ptr * a + b)
              << " using factors " << a << " " << b
              << std::endl;
#endif
    while (src_ptr < src_end)
      *dst_ptr++ = static_cast<typename dst_type::value_type>(*src_ptr++ * a + b);
    return dst;
  }
}

/**
 * See s_scaleToFloat() for details.
 */
template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::s_scaleFromFloat(const DataBuffer* in,
                                       const std::array<double,2>& factors)
{
  typedef DataBufferNd<float,NDim> src_type;
  typedef self_type dst_type;
  // Factors are inverted when converting from float user to int storage.
  const double a = 1.0 / factors[0]; // slope
  const double b = -factors[1] * a;  // intercept
  // TODO-Low, should I check for NaN and replace with defaultvalue,
  // or is that just too expensive?
  T defaultstorage = RoundAndClip<T>(0 * a + b);
  //T defaultvalue = RoundAndClip<T>(-factors[1]/factors[0]) * factors[0] + factors[1];
  auto src = dynamic_cast<const src_type*>(in);
  if (!src) {
    if (in == nullptr)
      return std::shared_ptr<dst_type>();
    else
      throw OpenZGY::Errors::ZgyInternalError("Wrong type in DataBufferNd::scaleFromFloat");
  }
  else if (!std::numeric_limits<T>::is_integer) {
    return std::shared_ptr<dst_type>();
  }
  else if (src->isScalar()) {
    T val = RoundAndClip<T>(src->scalarValue() * a + b, defaultstorage);
    return std::shared_ptr<dst_type>(new dst_type(val, src->safesize()));
  }
  else {
    auto dst = std::shared_ptr<dst_type>(new dst_type(src->safesize()));
    // TODO-Worry: Beware of edge bricks.
    // It might not be a good idea to convert the padding area.
    //
    // TODO-Performance: If caller has already computed the value
    // range, it might be able to tell us whether clipping or testing
    // for infinite or both is needed.
    dst_type::value_type *dst_ptr = dst->data();
    const typename src_type::value_type *src_ptr = src->data();
    const std::int64_t totalsize = src->allocsize();
    // TODO-Performance: There is a risk of the OpemMP overhead being larger
    // then the speedup gained by multiple threads. I ran tests only in one
    // environment. It seemed safe by a wide margin if there was at least one
    // standard-sized brick (256 KB to 1 MB) being processed. I am still
    // worrying, because technically there is no upper limit to how much
    // overhead OpenMP might add. While doing serial processing has a fixed
    // and not that dramatic cost.
#pragma omp parallel for if(totalsize >=256*1024)
    for (std::int64_t ii=0; ii<totalsize; ++ii) {
      dst_ptr[ii] = RoundAndClip<dst_type::value_type>(src_ptr[ii] * a + b, defaultstorage);
    }
    return dst;
  }
}

template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::scaleToFloat(const std::array<double,2>& factors)
{
  SimpleTimerEx tt(AdHocTimers::instance().scaletof);
  return s_scaleToFloat(this, factors);
}

template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::scaleToStorage(const std::array<double,2>& factors,
                                     RawDataType dt)
{
  SimpleTimerEx tt(AdHocTimers::instance().scaletos);
  // A better way of converting an enum to a template argument would be welcome.
  switch (dt) {
  case RawDataType::SignedInt8:    return DataBufferNd<std::int8_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::UnsignedInt8:  return DataBufferNd<std::uint8_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::SignedInt16:   return DataBufferNd<std::int16_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::UnsignedInt16: return DataBufferNd<std::uint16_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::SignedInt32:   return DataBufferNd<std::int32_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::UnsignedInt32: return DataBufferNd<std::uint32_t,NDim>::s_scaleFromFloat(this, factors);
  case RawDataType::Float32:       return std::shared_ptr<DataBuffer>();
  case RawDataType::IbmFloat32:    return std::shared_ptr<DataBuffer>();
  default: throw OpenZGY::Errors::ZgyInternalError("Unrecognized type enum");
  }
}

/**
 * Make a deep copy of the buffer. Works both for scalars and normal buffers.
 */
template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::clone() const
{
  SimpleTimerEx tt(AdHocTimers::instance().clone);
  std::shared_ptr<DataBufferNd<T,NDim>> result;
  if (this->isScalar()) {
    result.reset(new DataBufferNd<T,NDim>(this->scalarValue(), this->safesize()));
  }
  else {
    result.reset(new DataBufferNd<T,NDim>(this->safesize()));
    memcpy(result->data(), this->data(), this->allocsize() * this->itemsize());
  }
  return result;
}

/**
 * \brief Make a shallow view over the buffer, showing only part of the input.
 *
 * \detailed
 * Only the slowest varying dimension may be sliced. This means the output
 * will still be c-contiguous.
 *
 * If size[0] is 1 then it is allowed but not required to slice on dim=1,
 * and if size[1]==1 then it is allowed to slice on dim=2 as well.
 *
 * Currently this is only called from lodalgo.cpp
 *  createLodMT -> createLodPart -> slice1
 */
template <typename T, int NDim>
std::shared_ptr<DataBuffer>
DataBufferNd<T,NDim>::slice1(int dim, std::int64_t start, std::int64_t size) const
{
  std::shared_ptr<DataBufferNd<T,NDim>> result;
  ndsize_t newsize = this->_size;
  newsize[dim] = size;
  if (this->isScalar()) {
    result.reset(new DataBufferNd<T,NDim>(this->scalarValue(), newsize));
  }
  else {
    for (int dd = dim-1; dd >= 0; --dd)
      if (_size[dd] != 1)
        throw OpenZGY::Errors::ZgyInternalError
          ("DataBuffer::slice1() cannot slice on dim=" + std::to_string(dim)
           + " when size[" + std::to_string(dd) + "]="
           + std::to_string(_size[dd]) + ".");
    // Adjust the bulk ptr to start at the part we are interested in.
    T* ptr = _data.get() + start * _stride[dim];
    // Alias this to the original smart pointer.
    std::shared_ptr<T> smartptr(_data, ptr);
    result.reset(new DataBufferNd<T,NDim>(smartptr, newsize));
  }
  return result;
}

/**
 * For ad-hoc performance measurements.
 */
void
DataBuffer::cleartimers(bool show)
{
  AdHocTimers::instance().cleartimers(show);
}

// Explicit template instanciation.
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::int8_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::uint8_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::int16_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::uint16_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::int32_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<std::uint32_t,3>)
OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(DataBufferNd<float,3>)

}
// namespace
