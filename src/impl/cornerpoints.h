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
#pragma once

#include "declspec.h"

#include <stddef.h>
#include <array>
#include <cstdint>

/**
 * \file cornerpoints.h
 * \brief Deprecated. See transform.h
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * This class can be used to calculate the map projection (x, y) coordinates
 * of the four corners of a cube from a set of arbitrary control points (ACP).
 * The result is ordered according to the Petrel Ordered Corner Points (OCP)
 * definition, which is as follows corresponding to bulk data access indices:
 *
 *   (          0,           0, ?)
 *   (size[0] - 1,           0, ?)
 *   (          0, size[1] - 1, ?)
 *   (size[0] - 1, size[1] - 1, ?)
 *
 * Ref: PetrelOrientationHandling
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are meant to be short lived, used only inside some
 * method with no possibility of accessing it outside.
 */
class OPENZGY_TEST_API OrderedCornerPoints
{
public:

  /** Ordering */
  enum
  {
    Min0Min1, /** Lowest index in dimension 0 (I), lowest index in dimension 1 (J). */
    Max0Min1, /** Highest index in dimension 0 (I), lowest index in dimension 1 (J). */
    Min0Max1, /** Lowest index in dimension 0 (I), highest index in dimension 1 (J). */
    Max0Max1  /** Highest index in dimension 0 (I), highest index in dimension 1 (J). */
  };

  typedef std::int64_t  index_type; /** Bulk-data index type. */
  typedef float   annot_type; /** Inline/crossline annotation index datatype. */
  typedef double  coord_type; /** Map projection coordinate datatype. */

  /** Represents one specific corner. */
  struct OPENZGY_TEST_API Element
  {
    index_type  i;    /** Bulk-data index along first dimension. */
    index_type  j;    /** Bulk-data index along second dimension. */
    annot_type  il;   /** Inline annotation index. */
    annot_type  xl;   /** Crossline annotation index. */
    coord_type  x;    /** Map projection X (Easting) coordinate. */
    coord_type  y;    /** Map projection Y (Northing) coordinate. */

    /** Default constructor, initializes all members to zero. */
    Element();

    /** Constructor. Initializes i and j from args, remaining members are set to zero. */
    Element(index_type _i, index_type _j);

    /** Constructor. Initializes all members from args. */
    Element(index_type _i, index_type _j, annot_type _il, annot_type _xl, coord_type _x, coord_type _y);
  };

  /**
   * Construct OCPs from cube extent in annotation space and three ACPs.
   * @param il0     Inline annotation corresponding to bulk data index 0 along dimension 0 (I).
   * @param ilinc   Inline annotation increment corresponding to bulk data index increment of 1 along dimension 0 (I).
   * @param ilcnt   Number of samples along dimension 0 (I).
   * @param xl0     Crossline annotation corresponding to bulk data index 0 along dimension 1 (J).
   * @param xlinc   Crossline annotation increment corresponding to bulk data index increment of 1 along dimension 1 (J).
   * @param xlcnt   Number of samples along dimension 1 (J).
   * @param acp0il  Inline annotation for first ACP.
   * @param acp0xl  Crossline annotation for first ACP.
   * @param acp0x   X (easting) map coordinate for first ACP.
   * @param acp0y   Y (northing) map coordinate for first ACP.
   * @param acp1il  Inline annotation for second ACP.
   * @param acp1xl  Crossline annotation for second ACP.
   * @param acp1x   X (easting) map coordinate for second ACP.
   * @param acp1y   Y (northing) map coordinate for second ACP.
   * @param acp2il  Inline annotation for third ACP.
   * @param acp2xl  Crossline annotation for third ACP.
   * @param acp2x   X (easting) map coordinate for third ACP.
   * @param acp2y   Y (northing) map coordinate for third ACP.
   */
  OrderedCornerPoints
    (
      annot_type il0,    annot_type ilinc,  size_t ilcnt,
      annot_type xl0,    annot_type xlinc,  size_t xlcnt,
      annot_type acp0il, annot_type acp0xl, coord_type acp0x, coord_type acp0y,
      annot_type acp1il, annot_type acp1xl, coord_type acp1x, coord_type acp1y,
      annot_type acp2il, annot_type acp2xl, coord_type acp2x, coord_type acp2y
    );

  /**
   * Get an OCP.
   * @param i OCP index.
   * @return Element corresponding to ith OCP.
   */
  const Element& operator[](size_t i) const { return m_element[i]; }

  /** \brief corner coordinates in index space. */
  std::array<std::array<coord_type,2>,4> index_coords() const;
  /** \brief corner coordinates in annotation space. */
  std::array<std::array<coord_type,2>,4> annot_coords() const;
  /** \brief corner coordinates in world space. */
  std::array<std::array<coord_type,2>,4> world_coords() const;

private:

  std::array<Element,4> m_element;
};

} // end namespace
