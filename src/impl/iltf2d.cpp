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

#include "iltf2d.h"
#include <stdexcept>
#include <cmath>
#include <memory.h>

namespace InternalZGY {
#if 0
}
#endif

/*

NOTE: The math used (Trond Byberg, 2006-12-07)

  We represent the mapping R^2 -> R^2 in the form of a 3x3 homogeneous transform
  matrix T such that

    | v |     | x |
    |   |     |   |
    | w | = T*| y |                                                          (1)
    |   |     |   |
    | 1 |     | 1 |


  where T is implicitly defined by the known points:

    (x0, y0)  ->  (v0, w0)
    (x1, y1)  ->  (v1, w1)
    (x2, y2)  ->  (v2, w2)


  We now define normalized coordinate (s, t) such that for
  a general point (x, y) and its corresponding (v, w) we have


    | x |   | (x1 - x0)   (x2 - x0)    x0 | | s |     | s |
  	|   |   |                             | |   |     |   |
  	| y | = | (y1 - y0)   (y2 - y0)    y0 |*| t | = A*| t |                  (2)
  	|   |   |                             | |   |     |   |
  	| 1 |   | 0           0            1  | | 1 |     | 1 |

  and

    | v |   | (v1 - v0)   (v2 - v0)    v0 | | s |     | s |
  	|   |   |                             | |   |     |   |
  	| w | = | (w1 - w0)   (w2 - w0)    w0 |*| t | = B*| t |                  (3)
  	|   |   |                             | |   |     |   |
  	| 1 |   | 0           0            1  | | 1 |     | 1 |


  where A and B are known 3x3 matrices given by the tie-points.
  Solving (2) for (s, t) gives


    | s |          | x |
  	|   |          |   |
  	| t | = inv(A)*| y |                                                     (4)
  	|   |          |   |
  	| 1 |          | 1 |


  Combining (3) and (4) we can now tie (x, y) and (v, w) as


    | v |     | s |            | x |
  	|   |     |   |            |   |
  	| w | = B*| t | = B*inv(A)*| y |                                         (5)
  	|   |     |   |            |   |
  	| 1 |     | 1 |            | 1 |


  And (5) gives the solution to (1) as


    T = B*inv(A)                                                             (6)


  Writing out (6) explicitly becomes


        | (x1 - x0)   (x2 - x0)    x0 |
        |                             |
    T = | (y1 - y0)   (y2 - y0)    y0 |
        |                             |
        | 0           0            1  |

        |  (w2 - w0)   -(v2 - v0)    -v0*(w2 - w0) + w0*(v2 - v0) |
        |                                                         |
      * | -(w1 - w0)    (v1 - v0)    -w0*(v1 - v0) + v0*(w1 - w0) |
        |                                                         |
        |  0            0             1                           |

      * 1/((v1 - v0)*(w2 - w0) - (v2 - v0)*(w1 - w0))                        (7)

*/


/**
 * This class is used in the implementation of ImplicitLinearTransform2d.
 * \details Thread safety: Safe because it contains static methods only,
 */
class ImplicitLinearTransform2dImp
{
  ImplicitLinearTransform2dImp() = delete;
  ImplicitLinearTransform2dImp(const ImplicitLinearTransform2dImp&) = delete;
  ImplicitLinearTransform2dImp& operator=(const ImplicitLinearTransform2dImp&) = delete;
public:
  typedef double value_type;

  /**
   * Multiply two 3x3 homogeneous matrices (i.e. 2x3 explicit representation)
   * @param result  Matrix to write result to.
   * @param left    Matrix to serve as left-hand of multiplication.
   * @param right   Matrix to serve as right-hand of multiplication.
   */
  static void Mult(value_type result[2][3], const value_type left[2][3], const value_type right[2][3]);
};

void ImplicitLinearTransform2dImp::Mult(value_type result[2][3], const value_type left[2][3], const value_type right[2][3])
{
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 3; ++j) {

      // accumulate over the two first (explicit) rows of right
      result[i][j] = 0;
      for (size_t k = 0; k < 2; ++k) {
        result[i][j] += left[i][k]*right[k][j];
      }

      // add last (implicit) row of right (1 on diagonal, 0 off diagonal)
      if (j == 2) {
        result[i][j] += left[i][2];
      }
    }
  }
}

ImplicitLinearTransform2d::TiePoint::TiePoint(value_type a0, value_type a1, value_type b0, value_type b1)
{
  a[0] = a0;
  a[1] = a1;
  b[0] = b0;
  b[1] = b1;
}

ImplicitLinearTransform2d::ImplicitLinearTransform2d(const TiePoint& pt0,
                                                        const TiePoint& pt1,
                                                        const TiePoint& pt2)
{
  /*
  NOTE: See mathematical deduction at top of file for explanation.
  */

  // find the normalized-to-A and normalized-to-B homogeneous matrices
  value_type
    A[2][3] =
    {
      { pt1.a[0] - pt0.a[0], pt2.a[0] - pt0.a[0], pt0.a[0] },
      { pt1.a[1] - pt0.a[1], pt2.a[1] - pt0.a[1], pt0.a[1] }
    },
    B[2][3] =
    {
      { pt1.b[0] - pt0.b[0], pt2.b[0] - pt0.b[0], pt0.b[0] },
      { pt1.b[1] - pt0.b[1], pt2.b[1] - pt0.b[1], pt0.b[1] }
    };

  // find the inverse of A
  value_type
    Ainv[2][3] =
    {
      {  A[1][1], -A[0][1], -A[0][2]*A[1][1] + A[1][2]*A[0][1] },
      { -A[1][0],  A[0][0],  A[0][2]*A[1][0] - A[1][2]*A[0][0] }
    };
  {
    // check 2x2 determinant
    value_type
      Adet2x2 = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (Adet2x2 == 0) {
      // TODO-Low: Throw exception here?
      memset(m_matrix, 0, sizeof(m_matrix));
      return;
    }

    // divide by 2x2 determinant
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        Ainv[i][j] /= Adet2x2;
      }
    }
  }

  // calculate final homogeneous transform matrices by multiplying B*Ainv
  ImplicitLinearTransform2dImp::Mult(m_matrix, B, Ainv);
}

void ImplicitLinearTransform2d::Apply(value_type* b, const value_type* a) const
{
  // apply homogeneous transform matrix
  b[0] = m_matrix[0][0]*a[0] + m_matrix[0][1]*a[1] + m_matrix[0][2];
  b[1] = m_matrix[1][0]*a[0] + m_matrix[1][1]*a[1] + m_matrix[1][2];
}

void ImplicitLinearTransform2d::operator()(value_type* b, const value_type* a) const
{
  Apply(b, a);
}

} // end namespace
