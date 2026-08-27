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

#pragma once

#include <array>
#include <ostream>
#include <sstream>

/** \file: arrayops.h
 *
 * Extend std::array<> with numeric types to support arithmetic operations
 * between two instances (operations done per element) and between an
 * instance and a scalar (the scalar operates on each element).
 *
 * \li  a @  b      -- implemented for +.- * /
 * \li  a @= b      -- NOT implemented
 * \li  a @  scalar -- implemented for + - * /
 * \li  a @= scalar -- NOT implemented
 * \li  scalar @ a  -- NOT implemented
 *
 * TODO-Low: Should I write a new class and add behavior to that instead of
 * trying to "improve" std::array? In particular the question is relevant
 * because I might only ever need std::array<std::int64_t,3> so my new
 * class won't even need to be templated. The downside is that this
 * type may end up visible in the API and I don't want to introduce yet
 * another custom defined type if I can help it.
 *
 * To make use of these operators you will need to use:
 * \code using namespace InternalZGY::ArrayOps; \endcode
 */

namespace InternalZGY { namespace ArrayOps {
#if 0
}}
#endif

// Operations on arrays. TODO-Low ugly code smell.
// I probably need to make a dedicated "index3_t" type
// that will support these operators.

template<typename T, std::size_t N>
static std::array<T,N> operator+(const std::array<T,N>& a, const std::array<T,N>& b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] + b[dim];
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator-(const std::array<T,N>& a, const std::array<T,N>& b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] - b[dim];
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator*(const std::array<T,N>& a, const std::array<T,N>& b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] * b[dim];
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator/(const std::array<T,N>& a, const std::array<T,N>& b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] / b[dim];
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator+(const std::array<T,N>& a, T b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] + b;
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator-(const std::array<T,N>& a, T b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] - b;
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator*(const std::array<T,N>& a, T b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] * b;
  return result;
}

template<typename T, std::size_t N>
static std::array<T,N> operator/(const std::array<T,N>& a, T b)
{
  std::array<T,N> result;
  for (std::size_t dim=0; dim<N; ++dim)
    result[dim] = a[dim] / b;
  return result;
}

}} // namespace

namespace InternalZGY { namespace Formatters {
template<typename T, std::size_t N>
static std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a)
{
  os << "(";
  for (std::size_t dim=0; dim<N; ++dim)
    os << a[dim] << (dim==N-1 ? "" : ", ");
  os << ")";
  return os;
}

// Kludge for g++11, my logger trick no longer works.
// Will need to explicitly convert vector to string.
template<typename T, std::size_t N>
static const std::string fmt(const std::array<T,N>& a)
{
  std::stringstream ss;
  ss << a;
  return ss.str();
}

}} // namespace
