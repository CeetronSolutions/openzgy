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

#include "transform.h"
#include <cmath>

namespace InternalZGY {
#if 0
}
#endif

/*********************************************************************
 *
 * Purpose:
 *
 *      We have two coordinate systems (e.g. a grid coordinate
 *      system and a world (projection) coordinate system), but
 *      we don't know the general transformation between them.
 *
 *      What we do have is three arbitrary points, and their
 *      positions in both coordinate systems.
 *
 *      This routine uses these three points to transform
 *      points from one system to the other.
 *
 * Description:
 *
 *      Let p0, p1, p2 be the 3 points, and let a = p1 - p0,
 *      and b = p2 - p0. An arbitrary point q can then be
 *      written as: q = p0 + s*a + t*b; where (s,t) are
 *      the "p"-coordinates of point q.
 *
 *      (s,t) are independant of which coordinate system
 *      we use to express p, so if we solve for s,t in
 *      one system, we can then directly compute q in
 *      the other system.
 *
 *      To see why (s,t) don't depend on system, write the
 *      equation as: q - p0 = s*a + t*b.
 *      Let T(v) be the transformation (some combination of
 *      scaling and rotation) that takes a point from A to B.
 *      A linear transformation is defined by these rules:
 *       1. T(cv) = cT(v)     - invariant under multiplication with a constant
 *       2. T(v+w) = T(v) + T(w) - invariant under vector addition
 *
 *      If we apply this, we get:
 *        T(q-p0) = T(s*a + t*b)
 *                = T(s*a) + T(t*b)   - rule 2
 *                = s*T(a) + t*T(b)   - rule 1
 *      Ie; (s,t) are the "p"-coordinates of q-p0 in both
 *      coordinate systems. QED.
 *
 *      If we expand the above equation, we get:
 *        xq - xp0 = s*xa + t*xb
 *        yq - yp0 = s*ya + t*yb
 *      We then have two equations in two unknowns, which
 *      we solve using Cramer's rule, and then use
 *      (s,t) to compute q in the other system.
 *
 *********************************************************************/
OPENZGY_API bool generalTransform(
     double AX0, double AY0,
     double AX1, double AY1,
     double AX2, double AY2,
     double BX0, double BY0,
     double BX1, double BY1,
     double BX2, double BY2,
     double *X,  double* Y, std::size_t length)
{
  // Make everything relative to p0
  AX1 -= AX0; AY1 -= AY0;
  AX2 -= AX0; AY2 -= AY0;
  BX1 -= BX0; BY1 -= BY0;
  BX2 -= BX0; BY2 -= BY0;

  double det = AX1*AY2 - AX2*AY1; // The determinant

  if(std::abs(det) < 1.0e-6)
    return false; // colinear or coincident control points.

  for (std::size_t ii = 0; ii < length; ++ii) {
    double xq = X[ii] - AX0;
    double yq = Y[ii] - AY0;
    double s  = (xq*AY2 - AX2*yq)/det;
    double t  = (AX1*yq - xq*AY1)/det;
    X[ii] = BX0 + s*BX1 + t*BX2;
    Y[ii] = BY0 + s*BY1 + t*BY2;
  }
  return true;
}

} // namespace
