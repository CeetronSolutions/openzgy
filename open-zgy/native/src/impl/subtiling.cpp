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

#include "subtiling.h"

#include <cstdint>
#include <array>
#include <string.h>

namespace InternalZGY {
#if 0
}
#endif

/** \file: subtiling.cpp

\brief Handle 8x8 subtiling.

Version 1 of the uncompressed ZGY format reorganized the contents of alpha tiles
and data bricks before writing. Instead of keeping the 64x64(x64) layout, they
were changed to consisting of 8x8 subtiles/bricks each holding 8x8(x64) samples.
This was done to provide faster access to inline/crossline sections, as the
reader would then not need to load full 64x64(x64) tiles/bricks, but could
concentrate only on the much smaller (1/64th) 8x8(x64) subtiles/bricks that
intersected the section of interest.

Starting with version 2, such subtiling is no longer done.

When reading and writing files of version 1, the 8x8 subtiling is handled by
treating the data as being of higher dimension and using a general permutation
function. The rationale is as follows.

Consider a brick of size N[3]. In the two first dimensions it is subtiled
with tiles of size n[2]. The offset of sample i[3] is then:

\verbatim
  offset(i) = (i[0] DIV n[0]) * (n[0]*N[1]*N[2])
            + (i[1] DIV n[1]) * (n[0]*n[1]*N[2])
            + (i[0] MOD n[0]) * (n[1]*N[2])
            + (i[1] MOD n[1]) * (N[2])
            + (i[2])          * (1)
\endverbatim

where the two first lines correspond to indexing a specific subtile, and the
last three lines correspond to indexing within the subtile.

For the non-subtiled variant of the brick, the offset of sample i[3] is
simply:

\verbatim
  offset(i) = i[0]*N[1]*N[2] + i[1]*N[2] + i[2]
\endverbatim

This can be rearranged into

\verbatim
  offset(i) = (i[0] DIV n[0]) * (n[0]*N[1]*N[2])
            + (i[1] DIV n[1]) * (n[1]*N[2])
            + (i[0] MOD n[0]) * (N[1]*N[2])
            + (i[1] MOD n[1]) * (N[2])
            + (i[2])          * (1)
\endverbatim

which takes the same form as the offset for the subtiled brick.

Both cases can now be viewed as finding the offset of sample x[5] in a 5-
dimensional dataset of size { N[0]/n[0], N[1]/n[1], n[0], n[1], N[2] } with
stride S[5]:

\verbatim
  offset(x) = x[0]*S[0] + x[1]*S[1] + x[2]*S[2] + x[3]*S[3] + x[4]*S[4]
\endverbatim

where the indices in both cases are defined by

\verbatim
  x = { i[0] DIV n[0], i[1] DIV n[1], i[0] MOD n[0], i[1] MOD n[1], i[2] }
\endverbatim

and the strides are

\verbatim
  S = { n[0]*N[1]*N[2], n[0]*n[1]*N[2], n[1]*N[2], N[2], 1 }
  S = { n[0]*N[1]*N[2], n[1]*N[2], N[1]*N[2], N[2], 1 }
\endverbatim

for the subtiled and non-subtiled bricks, respectively.

With only strides differing (in zero-based element 1 and 2), we can then convert
between the two by applying a general 5-dimensional permutation. Since the last
two strides are equal, we can in practice achieve the correct permutation by
using only 4-dimensions. This corresponds to treating each vertical plane of
the subtiles as a 1-dimensional array.
*/

static void
permute(int ndim, const std::int64_t* size,
        const std::int64_t* dststride,       void* dstbuf_in,
        const std::int64_t* srcstride, const void* srcbuf_in,
        std::int64_t itemsize)
{
  // Enable pointer arithmetic
  char* dstbuf = static_cast<char*>(dstbuf_in);
  const char* srcbuf = static_cast<const char*>(srcbuf_in);
  if (ndim < 1) {
    // sanity check fail, do nothing.
  }
  else if (ndim == 1 && srcstride[0] == 1 && dststride[0] == 1) {
    // degenerate case: both source and destination stride is 1,
    // so copy all in one go.
    memcpy(dstbuf, srcbuf, size[0]*itemsize);
  }
  else if (ndim == 1) {
    // degenerate case: one-dimensional data
    // copy each element. Inefficient but should not get here
    // in the only use case (8x8 subtiling) that we have.
    // If it becomes a problem, use a switch on itemsize
    // 1, 2, 4, 8 to help the compiler.
    for (std::int64_t i = 0; i < size[0]*itemsize; i += itemsize)
      memcpy(&dstbuf[i*dststride[0]], &srcbuf[i*srcstride[0]], itemsize);
  }
  else {
    // general case: recursive
    for (std::int64_t i = 0; i < size[0]; ++i)
      permute(ndim - 1, size + 1,
              dststride + 1, dstbuf + i*dststride[0],
              srcstride + 1, srcbuf + i*srcstride[0],
              itemsize);
  }
}

/**
 * Convert a single brick from subtiled to standard (if remove=true)
 * or back.
 */
void
subtiling(const std::array<std::int64_t,3>& bricksize, std::int64_t itemsize,
          void *dst, const void* src, bool remove)
{
  const std::int64_t size[4]{
    bricksize[0]/8,
    bricksize[1]/8,
    8,
    8*bricksize[2]*itemsize};
  const std::int64_t standard_stride[4]{
    8*bricksize[1]*bricksize[2]*itemsize,
    8*bricksize[2]*itemsize,
    bricksize[1]*bricksize[2]*itemsize,
    1};
  const std::int64_t subtiled_stride[4]{
    8*bricksize[1]*bricksize[2]*itemsize,
    8*8*bricksize[2]*itemsize,
    8*bricksize[2]*itemsize,
    1};
  // Note that itemsize has already been accounted for in size and stride.
  // This means that permute() must copy the buffer as if it were char*.
  permute(4, size,
          remove ? standard_stride : subtiled_stride, dst,
          remove ? subtiled_stride : standard_stride, src,
          1);
}

} // namespace
