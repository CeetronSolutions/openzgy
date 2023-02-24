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

#include "../declspec.h"
#include "enum.h"
#include "types.h"

#include <cstdint>
#include <array>
#include <memory>
#include <string>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file databuffer.h
 * \brief Each DataBuffer instance represents some in memory data.
 */

/**
 * \brief Each DataBuffer instance represents some in memory data.
 *
 * Warning, creeping features. The class starter out as a trivial
 * pointer + 3d size struct but is growing dangerously fat.
 *
 * This data type represents an unsafe buffer pointer, plus:
 *
 *   - (template parameter) value type
 *   - (template parameter) number of dimensions
 *   - Size in each dimension and implied total buffer size.
 *   - Stride in each dimension to convert 1d/Nd, always c-contiguous.
 *   - Option to flag the entire buffer as constant-value and store
 *     just the constant.
 *
 * There is no knowledge of where in the survey this data belongs.
 *
 * Support for views and arbitrary strides have been removed. Buffer
 * layout is now always c-contiguous, i.e. the last index is varying
 * fastest and there are no holes.
 *
 * About the previous version of this class: The original design and
 * reasoning behind allowing non-contiguous buffers was that it would
 * allow skipping a buffer copy needed e.g. for reshaping the user's
 * buffer into brick-sized chunks and use a view instead. Much like
 * Python't numpy.ndarray.
 *
 * This was intended to improve performance. But the performance
 * effect is unknown and might even be detrimental because methods
 * that operate on the buffer contents would often need to take
 * strides into account.
 *
 * Worse, the current OpenZGY/C++ implementation has a large number of
 * places that just assume the buffer is c-contiguous in order to save
 * time. It is very unlikely that anybody will ever bother to change
 * this. And if somebody does, bringing back the fancy DataBuffer
 * features from git will be the least of the problems.
 *
 * All this is about the C++ implementation. In Python the slicing is
 * built into the language and is fully supported by numpy. So in
 * Python the OpenZGY code does make use of slices.
 *
 * DataBuffer is a non-templated base class to allow passing buffers
 * around without needing to have the calling code templated as well.
 *
 * copyFrom() is somewhat unsafe. It can verify (using dynamic_cast)
 * that "src" is of the same type as "this". But srcorig etc. are
 * passed as dumb pointers so the code cannot know whether they are
 * of the correct size. Maybe I am being too general here; maybe the
 * number of dimensions should always be 3?
 *
 * The actual bulk data is stored as a std::shared_ptr<T> which means
 * it can be passed around independantly of the DataBuffer instance.
 * The voidData() method returns the actual smart pointer as a
 * std::shared_ptr<void> which can be downcast to the correct type.
 * the data() method in the templated leaf type returns the raw
 * pointer. It is possible to also have methods returning void*
 * and std::shared_ptr<T> but don't add those unless actually used.
 *
 * Constructors that expect an external buffer are required to provide
 * it as a smart pointer. The DataBuffer will share ownership with the
 * caller. There is nothing that prevents the callers from passing a
 * shared_ptr with a no-op deleter, effectively removing reference
 * counting. But that will of course void their warranty. And make the
 * callers 100% responsible for the data being valid long enough.
 *
 * Unsafe external buffers are currently used for read() and write()
 * requests from the application code and for delivery from the file
 * back end. Reference counted external buffers are used when the
 * buffer was built from decompressed data.
 *
 * Each of those scenarios should be ok since no buffers should be
 * held onto after the functions return. TODO-Worry: Could async
 * read requests get delivered after the call to read() got aborted
 * via an exception? TODO-Worry: Might somebody decide to hold on to
 * a buffer in order to implement delayed write? In 99% of cases the
 * application's data buffer needs to be copied into a properly
 * reference counted buffer anyway. Maybe make that 100% by not allowing
 * short cuts when the application writes exactly one brick. Another
 * problem is if the writer decides to copy out brick at a time into
 * a reusable one brick buffer. This could also lead to the user's
 * buffer being held on to longer. Maybe make sure all bricks are
 * copied up front before any writing starts. This might also simplify
 * parallelized compression.
 *
 * Changes I have considered:
 *
 *   * Move all data members to the parent. _data will be declared as
 *     std::shared_ptr<void>, and data() will need to do a cast. Which
 *     is safe with respect to alignmemt as long as the cast is to the
 *     concrete type that the buffer initially had. Several of the
 *     methods in the DataBuffer base class can now be made non virtual
 *     if desired.
 *
 *   * TODO-Low: Remove the template on number of dimensions. This would
 *     make the API simpler. Current usage
 *     in OpenZGY is always NDim=3. If the Alpha tile feature is resurrected
 *     then there will be a need for 2d arrays, but can probably be handled
 *     as 2d arrays with nk=1.
 *     Not quite as elegant as the current template but might still be
 *     a good idea overall.
 *
 * Thread safety: Modification may lead to a data race. The users of
 * this class are responsible for any serializing. Normally this
 * should not be an issue, as users typically thread these instances
 * as they would treat a POD array.
 */
class OPENZGY_TEST_API DataBuffer
{
public:
  //** No copy or assign in the abstract base class.
  DataBuffer(const DataBuffer&) = delete;
  //** No copy or assign in the abstract base class.
  DataBuffer& operator=(const DataBuffer&) = delete;
  DataBuffer() = default;
  virtual ~DataBuffer();
  virtual std::string  toString()  const = 0;
  virtual std::shared_ptr<void>        voidData()        = 0;
  virtual std::shared_ptr<const void>  voidData()  const = 0;
  //** Number of elements in internal buffer. 1 if scalar, else totalsize().
  virtual std::int64_t allocsize() const = 0;
  //** Number of elements. Computed as the product of size() in each dimension.
  virtual std::int64_t totalsize() const = 0;
  //** Size in bytes of one element.
  virtual std::int64_t itemsize()  const = 0;
  //** Number of elements in each of the (assumed) 3 dimensions
  virtual std::array<std::int64_t,3> size3d() const = 0;
  //** Buffer layout in each of the (assumed) 3 dimensions
  virtual std::array<std::int64_t,3> stride3d() const = 0;
  //** True if buffer represents a constant value. data() will be length 1.
  virtual bool isScalar() const = 0;
  //** True if buffer is not scalar but samples are all the same.
  virtual bool isAllSame(const std::int64_t *used_in) const = 0;
  //** The constant value, if any, as a double. Unspecified if !isScalar().
  virtual double scalarAsDouble() const = 0;
  //** Convert the templated Type argument of the leaf class to an enum.
  virtual RawDataType datatype() const = 0;
  //** Set all samples to the same value.
  virtual void fill(double value) = 0;
  //** Make the instance unusable. May be needed if data is not ref counted.
  virtual void clear() = 0;
  //** Value range excluding infinites and nan.
  virtual std::pair<double, double> range(const std::int64_t *used_in = 0) const = 0;
  /** Corresponds to openzgy.impl._partialCopy(). To make templates work
   ** better the C++ version is an instance method with 'this' as destination.
   **/
  virtual void copyFrom(
            const DataBuffer* src,
            const std::int64_t *srcorig,const std::int64_t *dstorig,
            const std::int64_t *cpyorig,const std::int64_t *cpysize) = 0;
  /** Deep copy of a DataBuffer */
  virtual std::shared_ptr<DataBuffer> clone() const = 0;
  /** Convert DataBufferNd<T,N> to DataBufferNd<float,N>, hiding leaf types. */
  virtual std::shared_ptr<DataBuffer> scaleToFloat(const std::array<double,2>&) = 0;
  /** Convert DataBufferNd<float,N> to DataBufferNd<T,N>, hiding leaf types. */
  virtual std::shared_ptr<DataBuffer> scaleToStorage(const std::array<double,2>&, RawDataType) = 0;
  virtual std::shared_ptr<DataBuffer> slice1(int dim, std::int64_t start, std::int64_t size) const = 0;
  static std::shared_ptr<DataBuffer> makeDataBuffer3d(const std::shared_ptr<void>& raw, std::int64_t nbytes, const std::array<std::int64_t,3>& size, RawDataType dtype);
  static std::shared_ptr<DataBuffer> makeNewBuffer3d(const std::array<std::int64_t,3>& size, RawDataType dtype);
  static std::shared_ptr<DataBuffer> makeScalarBuffer3d(double scalar, const std::array<std::int64_t,3>& size, RawDataType dtype);
  static void cleartimers(bool show);
};

/**
 * See the non-template base class DataBuffer for details.
 *
 * Thread safety: Modification may lead to a data race. The users of
 * this class are responsible for any serializing. Normally this
 * should not be an issue, as users typically thread these instances
 * as they would treat a POD array.
 */
template <typename T, int NDim>
class DataBufferNd : public DataBuffer
{
public:
  typedef T value_type;
  typedef DataBufferNd<T, NDim> self_type;
  typedef std::array<std::int64_t,NDim> ndsize_t;
  enum { ndim = NDim };

  // TODO-Low: Group functions better.
public: // These are only available in the templated class.
  //** Default copy/assign would work, but don't offer it if I don't need it.
  DataBufferNd(const self_type&) = delete;
  //** Default copy/assign would work, but don't offer it if I don't need it.
  self_type& operator=(const self_type&) = delete;
  //** Buffer represents a constant value. No allocated data.
  DataBufferNd(T scalar, const std::array<std::int64_t,NDim>& size);
  //** Allocate a new buffer and take ownership of it.
  DataBufferNd(const std::array<std::int64_t,NDim>& size);
  //** External buffer, shared ownership. Beware we might get a "fake" smartptr.
  DataBufferNd(const std::shared_ptr<T>& data, const std::array<std::int64_t,NDim>& size);
  //** Plain, strongly typed pointer to data.
  T* data() { return _data.get(); }
  //** Plain, strongly typed pointer to data.
  const T* data() const { return _data.get(); }
  //** The constant value, if any. Return first value if !isScalar().
  T scalarValue() const { return data()[0]; }
  //** True if buffer represents a constant value.
  bool isScalar() const override { return _is_scalar; }
  //** True if buffer is not scalar but samples are all the same.
  bool isAllSame(const std::int64_t *used_in) const override;
  //** The constant value, if any, as a double.
  double scalarAsDouble() const override { return (double)scalarValue(); }
  RawDataType datatype() const override { return RawDataTypeTraits<T>::datatype; }
  void fill(double value) override;
  void clear() override;
  std::pair<double, double> range(const std::int64_t *used_in = 0) const override;

private:
  ndsize_t _size;
  ndsize_t _stride;
  bool _is_scalar;
  std::shared_ptr<T> _data;

  static ndsize_t cstride(const ndsize_t& size);
  static ndsize_t make_array(const std::int64_t *in);
  static std::int64_t make_product(const ndsize_t& size);

public:
  std::string  toString()  const override;
  std::shared_ptr<void> voidData()        override { return _data; }
  std::shared_ptr<const void> voidData()  const override { return _data; }
  std::int64_t allocsize() const override;
  std::int64_t totalsize() const override;
  std::int64_t itemsize() const override;
  std::array<std::int64_t,3> size3d() const override;
  std::array<std::int64_t,3> stride3d() const override;
  const ndsize_t& safesize() const {return _size;} // TODO: is nonvirtual deliberate?
  const ndsize_t& safestride() const {return _stride;} // TODO: is nonvirtual deliberate?
  std::shared_ptr<DataBuffer> clone() const override;
  std::shared_ptr<DataBuffer> scaleToFloat(const std::array<double,2>&) override;
  std::shared_ptr<DataBuffer> scaleToStorage(const std::array<double,2>&, RawDataType) override;
  std::shared_ptr<DataBuffer> slice1(int dim, std::int64_t start, std::int64_t size) const override;

  void copyFrom(const DataBuffer* src,
                const std::int64_t *srcorig,const std::int64_t *dstorig,
                const std::int64_t *cpyorig,const std::int64_t *cpysize) override;

  static void copySubset(const ndsize_t& srcorig, const self_type& src,
                         const ndsize_t& dstorig,       self_type& dst,
                         const ndsize_t& cpyorig, const ndsize_t& cpysize);

  static void copySubset(const ndsize_t& srcorig, const self_type& src,
                         const ndsize_t& dstorig,       self_type& dst);

public:
  // Logically private, but needs to be available from a DteBufferNd with different template arguments. Easier to make them public.
  static std::shared_ptr<DataBuffer> s_scaleToFloat(const DataBuffer*, const std::array<double,2>&);
  static std::shared_ptr<DataBuffer> s_scaleFromFloat(const DataBuffer*, const std::array<double,2>&);
};

// Explicit template instanciation.
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::int8_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::uint8_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::int16_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::uint16_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::int32_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<std::uint32_t, 3>)
OPENZGY_DECLARE_EXPLICIT_TEMPLATE(DataBufferNd<float, 3>)

} // namespace
