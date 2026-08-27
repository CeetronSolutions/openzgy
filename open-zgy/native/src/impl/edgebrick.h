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

#pragma once

#include "../exception.h"

#include <iterator>
#include <cstdint>
#include <array>

namespace InternalZGY {
#if 0
}
#endif

/**
 * This is a very simplified STL style iterator for EdgeBrick.
 * It only supports a forward iterator over its contents. It is probably
 * missing several details that even a forward iterator is expected to have.
 * Case in point: A templatized allocator, and a const version. I might need
 * to add the latter.
 *
 * Iterate over a region of interest in a brick. Storage has dim 2 a.k.a.
 * vertical fastest varying. The ROI runs from {0,0,0} up to but not including
 * the provided "roi" argument.
 *
 * The expected usage is collecting histogram and statistics for a brick that
 * has padding past the survey edge that must be excluded.
 */
template<typename T>
class BrickIterator
{
public:
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;
  typedef T value_type;
  typedef T* pointer;
  typedef T& reference;

  BrickIterator()
    : data_(nullptr)
    , data_ptr_(nullptr)
    , size_{ 0,0,0 }
    , roi_{ 0,0,0 }
    , pos_{ 0,0,0 }
  {
  }

  BrickIterator(pointer data,
    const std::array<std::int64_t, 3>& size,
    const std::array<std::int64_t, 3>& roi,
    bool atend)
    : data_(data)
    , data_ptr_(data + (atend ? size[0] * size[1] * size[2] : 0))
    , size_(size)
    , roi_(roi)
    , pos_{ (atend ? size[0] : 0), 0, 0 }
  {
    for (int dim = 0; dim < 3; ++dim)
      if (roi[dim] < 0 || size[dim] < 0 || roi[dim] > size[dim])
        throw OpenZGY::Errors::ZgyInternalError("Bad size in BrickIterator");
  }

  reference operator*() const { return *data_ptr_; }
  pointer operator->() const { return data_ptr_; }
  BrickIterator& operator++() { increment(); return *this; }
  BrickIterator operator++(int) { BrickIterator tmp = *this; ++(*this); return tmp; }
  friend bool operator== (const BrickIterator& a, const BrickIterator& b) { return a.data_ptr_ == b.data_ptr_; };
  friend bool operator!= (const BrickIterator& a, const BrickIterator& b) { return a.data_ptr_ != b.data_ptr_; };

private:
  pointer data_;
  pointer data_ptr_;
  std::array<std::int64_t, 3> size_;
  std::array<std::int64_t, 3> roi_;
  std::array<std::int64_t, 3> pos_;

  void increment()
  {
    ++data_ptr_;
    if (++(pos_[2]) >= roi_[2]) {
      pos_[2] = 0;
      data_ptr_ += (size_[2] - roi_[2]);
      if (++(pos_[1]) >= roi_[1]) {
        pos_[1] = 0;
        if (++(pos_[0]) >= roi_[0]) {
          pos_[0] = size_[0]; // roi end -> brick end.
        }
        data_ptr_ = data_ + pos_[0] * size_[1] * size_[2];
      }
    }
  }
};

/**
 * This is a very simplified STL container for a 3d buffer.
 * It only supports a forward iterator over its contents, limited to an ROI.
 * It is probably missing several details that a container is expected to have.
 * Case in point: A templatized allocator, and the possibility of returning
 * a const_itertor. I might need to add the latter. I an not sure an
 * EdgeBrick<const float>::begin() would work like EdgeBrick<float>::cbegin()
 *
 * Storage has dim 2 a.k.a. vertical fastest varying. The ROI runs from
 * {0,0,0} up to but not including the provided "roi" argument.
 *
 * The expected usage is collecting histogram and statistics for a brick that
 * has padding past the survey edge that must be excluded.
 */
template<typename T>
class EdgeBrick
{
public:
  typedef BrickIterator<T> iterator;

private:
  iterator begin_;
  iterator end_;

public:
  EdgeBrick(T* data,
    const std::array<std::int64_t, 3>& size,
    const std::array<std::int64_t, 3>& roi)
    : begin_(data, size, roi, false)
    , end_(data, size, roi, true)
  {
  }

  const iterator& begin() const { return begin_; }
  const iterator& end() const { return end_; }
};

} // namespace
