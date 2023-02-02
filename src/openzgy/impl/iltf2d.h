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

//Adapted from: Zgy/Common/iltf2d
#pragma once

/**
 * \file iltf2d.h
 * \brief Deprecated. See transform.h
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * This class represents a linear transform T: R^2 -> R^2 from 2-dimensional
 * Euclidean space A to 2-dimensional Euclidean space B, defined implicitly by
 * three given tie points.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are meant to be short lived, used only inside some
 * method with no possibility of accessing it outside.
 */
class ImplicitLinearTransform2d
{
public:

  typedef double value_type;                    /** Type of values processed by this class. */
  typedef ImplicitLinearTransform2d this_type;  /** Type of this class. */

  /**
   * Tie-point definition.
   */
  struct TiePoint
  {
    value_type a[2];  /** Coordinates in space A. */
    value_type b[2];  /** Coordinates in space B. */

    /**
     * Construct from arguments.
     * @param a0 Coordinate in first dimension of space A.
     * @param a1 Coordinate in second dimension of space A.
     * @param b0 Coordinate in first dimension of space B.
     * @param b1 Coordinate in second dimension of space B.
     */
    TiePoint(value_type a0, value_type a1, value_type b0, value_type b1);
  };

  /**
   * Construct from three tie-points. Fails by throwing an exception if
   * the tie-points don't span a 2-dimensional space (i.e. lie along a
   * straight line in either space).
   * @param pt0 First tie-point.
   * @param pt1 Second tie-point.
   * @param pt2 Third tie-point.
   */
  ImplicitLinearTransform2d(const TiePoint& pt0,
                            const TiePoint& pt1,
                            const TiePoint& pt2);

  /**
   * @param other Instance to compare with.
   * @return True if other is equal to self, otherwise false.
   */
  bool operator==(const this_type& other) const = delete;

  /**
   * @param other Instance to compare with.
   * @return False if other is equal to self, otherwise true.
   */
  bool operator!=(const this_type& other) const = delete;

  /**
   * Convert from space A to space B.
   * @param b Pointer to buffer holding (at least) two elements of value_type which will receive the transformed coordinate in B.
   * @param a Pointer to buffer holding (at least) two elements of value_type serving as the A coordinate to transform.
   */
  void Apply(value_type* b, const value_type* a) const;
  void operator()(value_type* b, const value_type* a) const;

private:

  value_type m_matrix[2][3]; /** Homogeneous transform matrix used in A-to-B conversions, third row is implicitly [0  0  1]. */
};

} // end namespace
