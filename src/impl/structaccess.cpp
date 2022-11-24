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

#include "structaccess.h"

/**
 * Convert a 64-bit integer stored on disk as vcs_uint64 to the in-memory
 * representation or back. Class vcs_uint64 consisted of two uint32_type
 * members representing the high and low (in that order) 32-bits of a 64-bit
 * unsigned integer. The two halves were stored little-endian on disk.
 *
 * This method includes conversion between big and little endian if needed,
 * byte swapping each of the 32-bit halves separately. Caveat, that
 * functonality has not been properly tested.
 *
 * The format was used by:
 *    - VCS (compressed ZGY) version 1 - 13
 *    - VBS (uncompressed ZGY) version 1
 *    - VCO (compressed object file) version 1
 *
 * Currently only VBS 1 is a concern for this code. There are 64 bit ints
 * in the alpha- and brick lookup tables and in the offset headers.
 */
void
InternalZGY::byteswapV1Long(std::int64_t *ptr, size_t n)
{
  union {
    std::int64_t i64;
    std::int32_t hilo[2];
  } in;
  while (n-- != 0) {
    in.i64 = *ptr;
    byteswapT(&in.hilo[0], 2);
    *ptr++ = (static_cast<std::int64_t>(in.hilo[0])<<32) | in.hilo[1];
  }
}

void
InternalZGY::byteswapV1Long(std::uint64_t *ptr, size_t n)
{
  union {
    std::uint64_t i64;
    std::uint32_t hilo[2];
  } in;
  while (n-- != 0) {
    in.i64 = *ptr;
    byteswapT(&in.hilo[0], 2);
    *ptr++ = (static_cast<std::uint64_t>(in.hilo[0])<<32) | in.hilo[1];
  }
}
