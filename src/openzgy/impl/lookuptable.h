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

#include "enum.h"

//#include <sstream>
//#include <algorithm>

#include <cstdint>
#include <vector>
#include <array>
#include <string>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file lookuptable.h
 * \brief Static methods to assist working with lookup tables.
 */

/**
 * \class LookupTable
 * \brief Static methods to assist working with lookup tables.
 *
 * \details I am uncertain whether I should extend this class.
 * Some alternatives are:
 *
 * \li (1) [current choice] the %LookupTable instance has no data members
 *     and only static methods. All context needs to be passed in.
 *     Which is tedious but makes it easier to see what the methods
 *     are actually doing.
 *
 * \li (2) The %LookupTable instance is meant to be short lived. It has
 *     references to the context it needs but no good memory management
 *     or ways to handle stale data.
 *
 * \li (3) The %LookupTable instance owns the actual lookup tables
 *     making this class more visible than just a helper for
 *     collecting lookup table related code. ZgyInternalMeta
 *     will hold a reference to a %LookupTable instead of the
 *     current pointer to an ILookupTableAccess
 *
 * \li (4) As (3) but use separate instances for alpha and brick.
 *     Make some changes to have the methods more generic.
 *     Implement the boilerplate code from ILookupTableAccess.
 *     The code then becomes a drop in replacement for LookupTableV0Access.
 *
 * Thread safety: With the current implementation the class only
 * contains static methods and no data. This means the class itself is
 * thread safe; the challenge is that the caller may need to protect
 * the data strucutures being passed in.
 */
class LookupTable
{
private:
#if 0
  // Actual lookup tables as read from file (alup, blup)
  // or computed (aend, bend for all current file versions)
  std::vector<std::uint64_t> _alup;
  std::vector<std::uint64_t> _aend;
  std::vector<std::uint64_t> _blup;
  std::vector<std::uint64_t> _bend;

  // Information derived from size, bricksize, and datatype.
  // This will never change so it is safe to cache.
  std::vector<std::array<std::int64_t,3>> _lodsizes;
  std::vector<std::int64_t> _alphaoffsets;
  std::vector<std::int64_t> _brickoffsets;
  std::int64_t _bytesperalpha;
  std::int64_t _bytesperbrick;

  // Scenario 4: If replacing the ILookupTableAccess
  // then some additional boilerplate code is needed.
  virtual void dump(std::ostream& out, const std::string& prefix = "");
  virtual void read(const std::shared_ptr<IFileADT>& file, std::int64_t offset, std::int64_t size);
  virtual void byteswap();
#endif

public:
  /** \brief Decoded contents of the lookup table for one brick or tile. */
  struct LutInfo
  {
    BrickStatus status;
    std::int64_t offset_in_file;
    std::int64_t size_in_file;
    std::uint32_t raw_constant; // Value stored in least significant bits.
    /** \brief Create an instance will all data fields filled in. */
    LutInfo(BrickStatus status_in, std::int64_t offset, std::int64_t size, std::uint32_t constant)
      : status(status_in)
      , offset_in_file(offset)
      , size_in_file(size)
      , raw_constant(constant)
    {
    }
  };

public:
  static std::vector<std::uint64_t> calcLookupSize(
      const std::vector<std::uint64_t>& lookup,
      std::int64_t eof, std::int64_t maxsize,
      bool *return_file_truncated,
      bool *return_bricks_overlap);

  static LutInfo getAlphaFilePosition(
      std::int64_t i, std::int64_t j, std::int64_t lod,
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& alphaoffsets,
      const std::vector<std::uint64_t>& alup,
      std::int64_t bytesperalpha);

  static LutInfo getAlphaFilePositionFromIndex(
      std::int64_t index,
      const std::vector<std::uint64_t>& alup,
      std::int64_t bytesperalpha);

  static LutInfo getBrickFilePosition(
      std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& brickoffsets,
      const std::vector<std::uint64_t>& blup,
      const std::vector<std::uint64_t>& bend,
      std::int64_t bytesperbrick);

  static LutInfo getBrickFilePositionFromIndex(
      std::int64_t index,
      const std::vector<std::uint64_t>& blup,
      const std::vector<std::uint64_t>& bend,
      std::int64_t bytesperbrick);

  static void setBrickFilePosition(
      std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
      const LutInfo& info,
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& brickoffsets,
      std::vector<std::uint64_t>* blup,
      std::vector<std::uint64_t>* bend);

  static std::int32_t usableBrickLOD(
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& brickoffsets,
      const std::vector<std::uint64_t>& blup,
      const std::vector<std::uint64_t>& bend);

  static bool hasBrickCompression(
      const std::vector<std::uint64_t>& blup,
      const std::vector<std::uint64_t>& bend);

  static void setNoBrickLOD(
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& brickoffsets,
      std::vector<std::uint64_t>* blup,
      std::vector<std::uint64_t>* bend);

  static std::int64_t getAlphaLookupIndex(
      std::int64_t i, std::int64_t j, std::int64_t lod,
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& alphaoffsets);

  static std::int64_t getBrickLookupIndex(
      std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
      const std::vector<std::array<std::int64_t,3>>& lodsizes,
      const std::vector<std::int64_t>& brickoffsets);

private:

  static std::string _formatPosition(
      std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod);

  static void _validatePosition(
      std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
      const std::vector<std::array<std::int64_t,3>>& lodsizes);

  static std::int64_t _getLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& offsets);

  static std::pair<std::int64_t, std::int64_t> _getBegAndSize(
      std::int64_t ix,
      const std::vector<std::uint64_t>& lup,
      const std::vector<std::uint64_t>& end,
      std::int64_t maxsize);
};
} // namespace
