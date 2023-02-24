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

#include "readwritemirror.h"
#include "readwriterelay.h"
#include "../impl/guid.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY { namespace Tools {
#if 0
}}
#endif

/**
 * Shamelessly duplicated from meta.cpp
 *
 * Compute the size of the file in bricks, for all LOD levels
 * and all 3 dimensions. Indirectly also compute the number
 * of LOD levels, since by definition the last LOD level
 * is the one that has size (1, 1, 1) i.e. the entire cube
 * is decimated enought to fit inside a single brick.
 */
static std::vector<std::array<std::int64_t,3>>
calcLodSizes(const std::array<std::int64_t,3>& size_in,
	     const std::array<std::int64_t,3>& bricksize)
{
  std::array<std::int64_t,3> size = size_in;
  std::vector<std::array<std::int64_t,3>> result;
  if (size[0] < 1 || size[1] < 1 || size[2] < 1) {
    // Empty file is the only case where the highest lod isn't one brick.
    // There will be one lod and its size will be zero bricks.
    result.push_back(std::array<std::int64_t,3>{0,0,0});
  }
  else {
    // number of full and partial bricks:
    for (int ii=0; ii<3; ++ii)
      size[ii] = (size[ii] + bricksize[ii] - 1) / bricksize[ii];
    result.push_back(size);
    // Continue to add new levels until we get a level with just one brick.
    while (size[0] > 1 || size[1] > 1 || size[2] > 1) {
      for (int ii=0; ii<3; ++ii)
        size[ii] = (size[ii] + 1) / 2;
      result.push_back(size);
    }
  }
  return result;
}

ZgyReaderMirror::ZgyReaderMirror(
     const std::shared_ptr<IZgyReader>& relay_in,
     const std::array<int,3>& mirrors)
  : ZgyReaderRelay(relay_in)
{
  const size3i_t bs = relay().bricksize();
  const size3i_t orig_size = relay().size();
  this->orig_cropped_size_ = size3i_t
    {
     (orig_size[0] / bs[0]) * bs[0],
     (orig_size[1] / bs[1]) * bs[1],
     (orig_size[2] / bs[2]) * bs[2]
    };

  this->inflated_size_ = size3i_t{
    this->orig_cropped_size_[0] * mirrors[0],
    this->orig_cropped_size_[1] * mirrors[1],
    this->orig_cropped_size_[2] * mirrors[2]};

  this->indexcorners_ = corners_t
    {{{                        0.0,                         0.0},
      {(double)inflated_size_[0]-1,                         0.0},
      {                        0.0, (double)inflated_size_[1]-1},
      {(double)inflated_size_[0]-1, (double)inflated_size_[1]-1}}};
  for (int ii=0; ii<4; ++ii)
    this->annotcorners_[ii] = relay().indexToAnnot(this->indexcorners_[ii]);
  for (int ii=0; ii<4; ++ii)
    this->worldcorners_[ii] = relay().indexToWorld(this->indexcorners_[ii]);

  // brickcount a.k.a. lodsizes, and nlods. Re-calculate from size.
  this->brickcount_ = calcLodSizes(inflated_size_, relay().bricksize());
  this->nlods_ = static_cast<std::int32_t>(this->brickcount_.size());

  // Several less important metadata, included because it is cheap.
  this->verid_ = InternalZGY::GUID::makeGUID();

  const int factor = mirrors[0] * mirrors[1] * mirrors[2];
  statistics_ = relay().statistics();
  statistics_.cnt *= factor;
  statistics_.sum *= factor;
  statistics_.ssq *= factor;

  histogram_ = relay().histogram();
  histogram_.samplecount *= factor;
  for (auto& it : histogram_.bins)
    it *= factor;

  filestats_ = relay().filestats();
  // I should have multiplied all the brick counts by factor.
  // But it is very unlikely that anybody will notice, or care.
  // Other members such as file size might not make much sense.
}

ZgyReaderMirror::~ZgyReaderMirror()
{
}

template<typename T> static void
flip(const T* inbuf, T* result,
     const std::array<std::int64_t,3>& size,
     const std::array<int,3>& replicant)
{
  if (((replicant[0] | replicant[1] | replicant[2]) & 1) == 0) {
    memcpy(result, inbuf, size[0]*size[1]*size[2]*sizeof(T));
  }
  else {
    const std::int64_t ibeg = (replicant[0] & 1) ? size[0] - 1 : 0;
    const std::int64_t jbeg = (replicant[1] & 1) ? size[1] - 1 : 0;
    const std::int64_t kbeg = (replicant[2] & 1) ? size[2] - 1 : 0;

    const std::int64_t iend = (replicant[0] & 1) ? -1 : size[0];
    const std::int64_t jend = (replicant[1] & 1) ? -1 : size[1];
    const std::int64_t kend = (replicant[2] & 1) ? -1 : size[2];

    const std::int64_t istep = (replicant[0] & 1) ? -1 : 1;
    const std::int64_t jstep = (replicant[1] & 1) ? -1 : 1;
    const std::int64_t kstep = (replicant[2] & 1) ? -1 : 1;

    for (std::int64_t iin = ibeg; iin != iend; iin += istep)
      for (std::int64_t jin = jbeg; jin != jend; jin += jstep)
        for (std::int64_t kin = kbeg; kin != kend; kin += kstep)
          *result++ = inbuf[iin * size[1] * size[2] + jin * size[2] + kin];
  }
}

template<typename T> static void flipI(T* buf, std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  if (buf != NULL)
    for (std::int64_t i=0; i<ni/2; ++i)
      for (std::int64_t j=0; j<nj; ++j)
        for (std::int64_t k=0; k<nk; ++k)
          std::swap(buf[i*nj*nk + j*nk + k], buf[(ni-i-1)*nj*nk + j*nk + k]);
}

template<typename T> static void flipJ(T* buf, std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  if (buf != NULL)
    for (std::int64_t j=0; j<nj/2; ++j)
      for (std::int64_t i=0; i<ni; ++i)
        for (std::int64_t k=0; k<nk; ++k)
          std::swap(buf[i*nj*nk + j*nk + k], buf[i*nj*nk + (nj-j-1)*nk + k]);
}

template<typename T> static void flipK(T* buf, std::int64_t ni, std::int64_t nj, std::int64_t nk)
{
  if (buf != NULL)
    for (std::int64_t i=0; i<ni; ++i)
      for (std::int64_t j=0; j<nj; ++j)
        for (std::int64_t k=0; k<nk/2; ++k)
          std::swap(buf[i*nj*nk + j*nk + k], buf[i*nj*nk + j*nk + (nk-k-1)]);
}

/**
 * Slightly slower flipping, but does not require a separate buffer.
 */
template<typename T> static void
flipIJK(T* slice,
     const std::array<std::int64_t,3>& size,
     const std::array<int,3>& replicant)
{
  if (replicant[0] & 1) flipI(slice, size[0], size[1], size[2]);
  if (replicant[1] & 1) flipJ(slice, size[0], size[1], size[2]);
  if (replicant[2] & 1) flipK(slice, size[0], size[1], size[2]);
}

std::array<int,3>
ZgyReaderMirror::getReplicant(const size3i_t& start, const size3i_t& size) const
{
  std::array<int,3> result{0,0,0};
  for (int ii=0; ii<3; ++ii) {
    const std::int64_t orig_size = this->orig_cropped_size_[ii];
    if (start[ii] < 0 || size[ii] <= 0)
      throw std::runtime_error("ZgyReaderMirror: Bad start or size");
    if ((start[ii] % size[ii]) != 0)
      throw std::runtime_error("ZgyReaderMirror: Misaligned read");
    const int lo = static_cast<int>(start[ii] / orig_size);
    const int hi = static_cast<int>((start[ii] + size[ii] - 1) / orig_size);
    if (lo != hi)
      throw std::runtime_error("ZgyReaderMirror: Crossing boundary");
    result[ii] = lo;
  }
  return result;
}

std::array<std::int64_t,3>
ZgyReaderMirror::getStart(
     const size3i_t& start,
     const size3i_t& size,
     const std::array<int,3>& replicant) const
{
  std::array<std::int64_t,3> result {
    start[0] % orig_cropped_size_[0],
    start[1] % orig_cropped_size_[1],
    start[2] % orig_cropped_size_[2]
  };
  for (int dim = 0; dim < 3; ++dim) {
    if (replicant[dim] & 1) {
      result[dim] = this->orig_cropped_size_[dim] - result[dim] - size[dim];
      // Earlier checks should have caugth this.
      if (result[dim] < 0 || result[dim] + size[dim] > orig_cropped_size_[dim])
        throw std::runtime_error("ZgyReaderMirror::getStart Internal error");
    }
  }
  return result;
}

void
ZgyReaderMirror::read(
     const size3i_t& start, const size3i_t& size,
     float* data, int lod) const
{
  std::array<int,3> rep = getReplicant(start, size);
  relay().read(getStart(start, size, rep), size, data, lod);
  flipIJK(data, size, rep);
}

void
ZgyReaderMirror::read(
     const size3i_t& start, const size3i_t& size,
     std::int16_t* data, int lod) const
{
  std::array<int,3> rep = getReplicant(start, size);
  relay().read(getStart(start, size, rep), size, data, lod);
  flipIJK(data, size, rep);
}

void
ZgyReaderMirror::read(
     const size3i_t& start, const size3i_t& size,
     std::int8_t* data, int lod) const
{
  std::array<int,3> rep = getReplicant(start, size);
  relay().read(getStart(start, size, rep), size, data, lod);
  flipIJK(data, size, rep);
}

std::pair<bool,double>
ZgyReaderMirror::readconst(
     const size3i_t& start, const size3i_t& size,
     int lod, bool as_float) const
{
  std::array<int,3> rep = getReplicant(start, size);
  return relay().readconst(getStart(start, size, rep), size, lod, as_float);
}

ZgyWriterMirror::ZgyWriterMirror(
     const std::shared_ptr<IZgyWriter>& relay_in,
     const std::array<int,3>& mirrors)
  : ZgyWriterRelay(relay_in)
  , mirrors_(mirrors)
{
  // ZgyWriterMirror already forces the original size
  // to align with the bricksize. If adding more code
  // e.g. to allow mirroring to succeed as long as
  // it is aligned to chunk size and maybe not brick
  // size, then this might not work. But why on earth
  // (or below it, in this case) would anybody do that?

  // The real writer already shows the inflated size,
  // need to trim it down to only show the non-mirrored
  // part. Because that is the only part the application
  // needs to write.
  const size3i_t orig_size = relay().size();
  this->orig_cropped_size_ = size3i_t
    {
      orig_size[0] / mirrors[0],
      orig_size[1] / mirrors[1],
      orig_size[2] / mirrors[2]
    };
}

ZgyWriterMirror::~ZgyWriterMirror()
{
}

void ZgyWriterMirror::missing(const std::string& name)
{
  throw std::runtime_error("ZgyWriterMirror::" + name + " is not implemented yet");
}

static void
writeit(
     const std::array<std::int64_t,3>& read_start,
     const std::array<std::int64_t,3>& read_size,
     const std::array<std::int64_t,3>& orig_size,
     const std::array<int,3>& mirror,
     std::function<void(const std::array<std::int64_t,3>& start, const std::array<int,3>& rep)> fn)
{
  std::array<int,3> rep{};
  for (rep[0] = 0; rep[0] < mirror[0]; ++rep[0]) {
    for (rep[1] = 0; rep[1] < mirror[1]; ++rep[1]) {
      for (rep[2] = 0; rep[2] < mirror[2]; ++rep[2]) {
        std::array<std::int64_t,3> this_start{};
        for (int dim = 0; dim < 3; ++dim) {
          if (rep[dim] & 1) {
            this_start[dim] =
              ((rep[dim]+1) * orig_size[dim]) -
              (read_start[dim] + read_size[dim]);
          }
          else {
            this_start[dim] =
              rep[dim] * orig_size[dim] + read_start[dim];
          }
        }
        fn(this_start, rep);
      }
    }
  }
}

void ZgyWriterMirror::write(const size3i_t& start, const size3i_t& size, const float* data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>& rep)
          {
            // Modifying the buffer is ok as long as we flip it back.
            flipIJK(const_cast<float*>(data), size, rep);
            relay().write(this_start, size, data);
            flipIJK(const_cast<float*>(data), size, rep);
          });
}

void ZgyWriterMirror::write(const size3i_t& start, const size3i_t& size, const std::int16_t *data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>& rep)
          {
            flipIJK(const_cast<std::int16_t*>(data), size, rep);
            relay().write(this_start, size, data);
            flipIJK(const_cast<std::int16_t*>(data), size, rep);
          });
  }

void ZgyWriterMirror::write(const size3i_t& start, const size3i_t& size, const std::int8_t* data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>& rep)
          {
            flipIJK(const_cast<std::int8_t*>(data), size, rep);
            relay().write(this_start, size, data);
            flipIJK(const_cast<std::int8_t*>(data), size, rep);
          });
}

void ZgyWriterMirror::writeconst(const size3i_t& start, const size3i_t& size, const float* data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>&)
          {
            relay().writeconst(this_start, size, data);
          });
}

void ZgyWriterMirror::writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>&)
          {
            relay().writeconst(this_start, size, data);
          });
}

void ZgyWriterMirror::writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data)
{
  writeit(start, size, this->orig_cropped_size_, this->mirrors_,
          [&](const size3i_t& this_start, const std::array<int,3>&)
          {
            relay().writeconst(this_start, size, data);
          });
}

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
