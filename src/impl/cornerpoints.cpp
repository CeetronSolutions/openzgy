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

//Adapted from: Zgy/FileAccess/OrderedCornerPoints

#include "cornerpoints.h"
#include "iltf2d.h"

#include <memory.h>

namespace InternalZGY {
#if 0
}
#endif

OrderedCornerPoints::Element::Element()
  : i(0), j(0), il(0), xl(0), x(0), y(0)
{
}

OrderedCornerPoints::Element::Element(index_type _i, index_type _j)
  : i(_i), j(_j), il(0), xl(0), x(0), y(0)
{
}

OrderedCornerPoints::Element::Element(index_type _i, index_type _j, annot_type _il, annot_type _xl, coord_type _x, coord_type _y)
  : i(_i), j(_j), il(_il), xl(_xl), x(_x), y(_y)
{
}

OrderedCornerPoints::OrderedCornerPoints
(
  annot_type il0,    annot_type ilinc,  size_t ilcnt,
  annot_type xl0,    annot_type xlinc,  size_t xlcnt,
  annot_type acp0il, annot_type acp0xl, coord_type acp0x, coord_type acp0y,
  annot_type acp1il, annot_type acp1xl, coord_type acp1x, coord_type acp1y,
  annot_type acp2il, annot_type acp2xl, coord_type acp2x, coord_type acp2y
)
{
  // create linear transform from annotation to coordinates
  ImplicitLinearTransform2d
    iltf2d(
      ImplicitLinearTransform2d::TiePoint(acp0il, acp0xl, acp0x, acp0y),
      ImplicitLinearTransform2d::TiePoint(acp1il, acp1xl, acp1x, acp1y),
      ImplicitLinearTransform2d::TiePoint(acp2il, acp2xl, acp2x, acp2y)
    );

  // fill in (i, j) indices for each corner
  m_element[Min0Min1] = Element(0,         0);
  m_element[Max0Min1] = Element(ilcnt - 1, 0);
  m_element[Min0Max1] = Element(0,         xlcnt - 1);
  m_element[Max0Max1] = Element(ilcnt - 1, xlcnt - 1);

  // calculate in (inline, crossline) and (x, y) for each corner
  for (size_t i = 0; i < 4; ++i) {

    // find (inline, crossline) annotation from (i, j) indices
    m_element[i].il = il0 + ilinc*m_element[i].i;
    m_element[i].xl = xl0 + xlinc*m_element[i].j;

    // find (x, y) coordinates from (inline, crossline) annotation using the implicit linear transform
    const ImplicitLinearTransform2d::value_type
      index[2] = { m_element[i].il, m_element[i].xl };  /** holds (inline, crossline) */
    ImplicitLinearTransform2d::value_type
      coord[2] = { 0, 0};                               /** receives (x, y) */
    iltf2d(coord, index);
    m_element[i].x = coord[0];
    m_element[i].y = coord[1];
  }
}

std::array<std::array<OrderedCornerPoints::coord_type,2>,4>
OrderedCornerPoints::index_coords() const
{
  return std::array<std::array<coord_type,2>,4> {{
    {(double)m_element[0].i, (double)m_element[0].j},
    {(double)m_element[1].i, (double)m_element[1].j},
    {(double)m_element[2].i, (double)m_element[2].j},
    {(double)m_element[3].i, (double)m_element[3].j}
  }};
}

std::array<std::array<OrderedCornerPoints::coord_type,2>,4>
OrderedCornerPoints::annot_coords() const
{
  return std::array<std::array<coord_type,2>,4> {{
    {m_element[0].il, m_element[0].xl},
    {m_element[1].il, m_element[1].xl},
    {m_element[2].il, m_element[2].xl},
    {m_element[3].il, m_element[3].xl}
  }};
}

std::array<std::array<OrderedCornerPoints::coord_type,2>,4>
OrderedCornerPoints::world_coords() const
{
  return std::array<std::array<coord_type,2>,4> {{
    {m_element[0].x, m_element[0].y},
    {m_element[1].x, m_element[1].y},
    {m_element[2].x, m_element[2].y},
    {m_element[3].x, m_element[3].y}
  }};
}

} // end namespace
