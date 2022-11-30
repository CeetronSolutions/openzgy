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

#include "declspec.h"

#include <array>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <type_traits>

/**
 * \file structaccess.h
 * \brief Low level tools for data conversion.
 */

namespace InternalZGY { // namespace StructAccess {
#if 0
}
#endif

/**
 * Convert a raw C pointer, possibly misaligned, to a safer std::array.
 * Callers will need to specify both template parameters. In theory
 * the compiler should be able to deduce T, but I haven't found a way
 * to do this while still explicitly providing N.
 */
template<typename T, std::size_t N>
std::array<T, N>
ptr_to_array(const T* in)
{
  std::array<T, N> result;
  memcpy(result.data(), in, sizeof(result));
  return result;
}

/**
 * Return the input scalar in such a way that &in is allowed to be
 * misaligned with respect to T. Note that I am not sure this is
 * good enough. A simple assignment instead of the memcpy will
 * trigger undefined behavior. A decent compiler should optimize
 * away the memcpy but also make sure that it doesn't use e.g.
 * not try to use e.g. SSE instructions that require alignment.
 * Worry: By defining the input argument as T&, will the compiler
 * already at that point make up its mind that in must be aligned?
 * https://pzemtsov.github.io/2016/11/06/bug-story-alignment-on-x86.html
 */
template<typename T>
T align(const T& in)
{
  T result;
  memcpy(&result, &in, sizeof(T));
  return result;
}

/**
 * For debugging, works as std::to_string but takes care to also
 * output std::int8_t and std::uint8_t as numbers, even though those
 * types are probably typedef'd to signed/unsigned char.
 */
template<typename T, std::size_t N>
std::string
array_to_string(const std::array<T, N>& a)
{
  std::stringstream ss;
  ss << "(";
  for (std::size_t ii=0; ii<N; ++ii)
    ss << std::to_string(a[ii]) << (ii == N-1 ? ")" : ", ");
  return ss.str();
}

/**
 * For debugging, format a std::array of a numeric type into a hex string.
 */
template<typename T, std::size_t N>
std::string
array_to_hex(const std::array<T, N>& a)
{
  std::stringstream ss;
  ss << std::hex << std::setfill('0') << std::right;
  for (std::size_t ii=0; ii<N; ++ii)
    ss << std::setw(2*sizeof(T)) << int(a[ii]) << (ii==N-1 ? "" : ",");
  return ss.str();
}

/**
 * For debugging, format a raw pointer to a numeric type into a string.
 * Caller needs to explicitly provide both template parameters.
 */
template<typename T, std::size_t N>
std::string
ptr_to_string(const T* a)
{
  return array_to_hex<T, N>(ptr_to_array<T, N>(a));
}

/**
 * For debugging, format a raw pointer to a numeric type into a hex string.
 * Caller needs to explicitly provide both template parameters.
 */
template<typename T, std::size_t N>
std::string
ptr_to_hex(const T* a)
{
  return array_to_hex<T, N>(ptr_to_array<T, N>(a));
}

/**
 * Byte swap between little-endian (as stored on file) and big-endian
 * (as used in some machines). The same function can be used both ways.
 * Try to make this work safely even if the pointer is misaligned.
 *
 * TODO-Low: Support for big-endian machines is experimental and there
 * are currently no tests. In other words, by definition it doesn't work.
 *
 * Avoid <endian.h> and <byteswap.h> for maximum portability. Use memcpy
 * to (hopefully) solve the misalignment issues.
 */
template<typename T>
static void
byteswapAlways(T *ptr, size_t n = 1)
{
  char in[sizeof(T)], out[sizeof(T)];
  for (size_t offset=0; offset<n; ++offset) {
    memcpy(&in[0], ptr, sizeof(T));
    for (size_t ii=0; ii<sizeof(T); ++ii)
      out[sizeof(T)-ii-1] = in[ii];
    memcpy(ptr, &out[0], sizeof(T));
    ++ptr;
  }
}

/**
 * Byte swap between little-endian (as stored on file) and big-endian
 * If and only if the current architecture is big-endian.
 */
template<typename T>
static void
byteswapT(T *ptr, size_t n = 1)
{
#if BIG_ENDIAN_ARCH
  byteswapAlways(ptr, n);
#endif
}

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
extern OPENZGY_TEST_API void byteswapV1Long(std::int64_t *ptr, size_t n = 1);

extern OPENZGY_TEST_API void byteswapV1Long(std::uint64_t *ptr, size_t n = 1);

template<typename T, typename U, int N>
std::array<T,N>
array_cast(const std::array<U,N>& in)
{
  std::array<T,N> result;
  for (int ii=0; ii<N; ++ii)
    result[ii] = static_cast<T>(in[ii]);
  return result;
}

} // namespace
