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
#include "declspec.h"
#include <cstdint>
#include <cstddef>

/**
 * \file transform.h
 * \brief General coordinate conversion based on 3 control points.
 * INTERNAL, but exported for use in vds2zgy.
 */
namespace InternalZGY {

  /**
   * \brief General coordinate conversion based on 3 control points.
   */
  OPENZGY_API extern bool generalTransform(
     double AX0, double AY0,
     double AX1, double AY1,
     double AX2, double AY2,
     double BX0, double BY0,
     double BX1, double BY1,
     double BX2, double BY2,
     double *X,  double* Y, std::size_t length);

} // namespace
