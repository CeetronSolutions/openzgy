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

#include "lookuptable.h"
#include "exception.h"
#include "enum.h"

#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <array>
#include <cstdint>
#include <limits>

using OpenZGY::Errors::ZgyUserError;
using OpenZGY::Errors::ZgyInternalError;
using OpenZGY::Errors::ZgyFormatError;

namespace InternalZGY {
#if 0
}
#endif

namespace
{
struct TmpLookupEntry
{
  std::uint64_t offset;
  std::uint64_t endpos;
  std::uint64_t type;
  std::uint64_t ordinal;
  TmpLookupEntry() : offset(0), endpos(0), type(0), ordinal(0) {}
};
}

/**
 * Given an index => start_offset lookup table, produce an
 * index => end_offset table by assuming there are no holes
 * in the allocated data.
 *
 * The function understands constant-value and compressed blocks.
 *
 * If eof and maxsize are known, the code can also make the
 * following checks:
 *
 * Blocks that have a start offset > eof are unreadable and
 * should be ignored. Set them to start and end at eof.
 * The same applies to offsets of unknown type i.e. the most
 * significant bit is 1 but the most significant byte is
 * neither 0x80 (constant) nor 0xC0 (compressed).
 * Blocks ending past eof should be assumed to end at eof.
 *
 * Blocks that appear to be larger than an uncompressed block are
 * probably too large. This may be caused by holes in the allocated
 * data. Assume the block is the same size as an uncompressed block.
 * If a compressed block takes up more room than an uncompressed one
 * then the writer should simply refrain from compressing it.
 * But for extra robustness the code that makes use of this
 * information should be prepared to retry the access of the block
 * really turned out to be larger.
 *
 * This method might be called unconditionally on file open, or
 * called only if at least one compressed brick was found, or it
 * might be deferred until the first time we read a compressed brick.
 *
 * TODO-Low: If alpha tiles are present then both brick and alpha offsets
 * ought to be considered in the same pass. The way it will work now
 * is that for bricks, a few bricks will appear too large because
 * they are followed by some alpha tiles. This is harmless.
 * For aplha tiles the end offsets will be hopelessly wrong.
 * We will need to just assume 4 KB for those.
 *
 * TODO-Performance: The brick-end calculation can probably be skipped
 * if the file has no compressed bricks. In that case the test for truncated
 * file needs to be moved (fairly easy) and we would lose the test for
 * overlapped bricks.
 */
std::vector<std::uint64_t>
LookupTable::calcLookupSize(
    const std::vector<std::uint64_t>& lookup,
    std::int64_t eof, std::int64_t maxsize,
    bool *return_file_truncated,
    bool *return_bricks_overlap)
{
  // [out] arguments
  if (return_file_truncated != nullptr)
    *return_file_truncated = false;
  if (return_bricks_overlap != nullptr)
    *return_bricks_overlap = false;

  std::vector<TmpLookupEntry> entries;
  std::uint64_t ord = 0;
  const std::uint64_t codeshift(56);
  const std::uint64_t codemask(((std::uint64_t)0xFF) << codeshift);
  for (const auto& offset : lookup) {
    TmpLookupEntry tmp;
    tmp.offset = offset;
    tmp.ordinal = ord++;
    tmp.type = (offset >> codeshift) & 0xFF;
    switch(tmp.type) {
    case 0x00: break;
    case 0x80: tmp.offset = 0; break;
    case 0xC0: tmp.offset &= (~codemask); break;
    default: tmp.offset = 0; break; // ignore blocks of unknown type.
    }
    tmp.endpos = 0;
    entries.push_back(tmp);
  }

  std::sort(entries.begin(), entries.end(),
            [](const TmpLookupEntry& a, const TmpLookupEntry& b) {
              return a.offset < b.offset;
            });

  // Technically I could also check the end of the last block, but
  // that only works for uncompressed and only if maxsize is known.
  //minimum_file_eof = entries.back().offset;
  if (eof != 0 && !entries.empty() && entries.back().offset >= (std::uint64_t)eof)
    if (return_file_truncated != nullptr)
      *return_file_truncated = true;

  // The end of block i is the start of block i+1,
  // except the last block which ends at EOF, just use a huge number.
  // And except all blocks that are not real (i.e. offset 0)
  // which should all end at 0 as well. For consistency,
  for (auto it = entries.begin() + 1; it < entries.end(); ++it) {
    if ((it-1)->offset != 0) {
      if ((it-1)->offset == it->offset) {
        if (return_bricks_overlap != nullptr)
          *return_bricks_overlap = true;
        // TODO-Low: This is currently a fatal error. It might later
        // be allowed for testing purposes to have the same block
        // pointed to by multiple entries. In that case more work is
        // needed to compute endpos.
      }
      (it-1)->endpos = it->offset;
    }
  }

  entries.back().endpos = eof ? eof : std::numeric_limits<std::int64_t>::max();

  std::sort(entries.begin(), entries.end(),
            [](const TmpLookupEntry& a, const TmpLookupEntry& b) {
              return a.ordinal < b.ordinal;
            });

  if (eof != 0 && maxsize != 0) {
    for (auto it = entries.begin(); it < entries.end(); ++it) {
      // End can neither go past eof nor indicate a size > maxsize.
      // Note that maxsize cannot be max for its type; that will overflow.
      // TODO-Low decide on unsigned or signed, not a mix.
      it->endpos = std::min(it->endpos, std::min((std::uint64_t)eof, it->offset + (std::uint64_t)maxsize));
      // If size appears negative then make it zero.
      it->endpos = std::max(it->endpos, it->offset);
    }
  }

  std::vector<std::uint64_t> result;
  result.reserve(entries.size());
  for (const auto& it : entries)
    result.push_back(it.endpos);
  return result;
}

/**
 * Convert i,j,k + lod to a string for error reporting or logging.
 */
std::string
LookupTable::_formatPosition(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod)
{
  std::stringstream ss;
  ss << "(" << i << ", " << j << ", " << k << ") lod " << lod;
  return ss.str();
}

/**
 * Check that i,j,k,lod are all inside valid bounds. Throw if not.
 * If used to validate an alpha tile then k should be passed as 0.
 */
void
LookupTable::_validatePosition(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes)
{
  if (lod < 0 || i<0 || j<0 || k<0)
    throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                       " cannot be negative");
  const std::int64_t nlods = static_cast<std::int64_t>(lodsizes.size());
  if (lod >= nlods)
    throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                       " cannot have lod >= " + std::to_string(nlods));
  const std::array<std::int64_t,3>& size = lodsizes[lod];
  if (i >= size[0] || j >= size[1] || k >= size[2])
    throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                       " cannot be >= " +
                       _formatPosition(size[0], size[1], size[2], nlods-1));
}

std::int64_t
LookupTable::_getLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& offsets)
{
  _validatePosition(i, j, k, lod, lodsizes);
  const std::array<std::int64_t,3>& size = lodsizes[lod];
  const std::int64_t index = (offsets[lod] +
                              i +
                              (size[0] * j) +
                              (size[0] * size[1] * k));
  if (index < 0 || index >= offsets.back())
    throw ZgyInternalError("Internal error in _getLookupIndex: " +
                           _formatPosition(i, j, k, lod) +
                           " result " + std::to_string(index) +
                           " max " + std::to_string(offsets.back()));
  return index;
}

/**
 * Normally called only from inside LookupTable, but can be useful
 * for dealing with other per alpha information such as dirty flags.
 */
std::int64_t
LookupTable::getAlphaLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& alphaoffsets)
{
  return _getLookupIndex(i, j, 0, lod, lodsizes, alphaoffsets);
}

/**
 * Normally called only from inside LookupTable, but can be useful
 * for dealing with other per brick information such as dirty flags.
 */
std::int64_t
LookupTable::getBrickLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets)
{
  return _getLookupIndex(i, j, k, lod, lodsizes, brickoffsets);
}

/**
 * Use only for bricks known to be compressed. Unless testing.
 * Return both the file offset of the compressed block and its size.
 * The latter is a hint that should not be fully trusted.
 *
 * For testing pass a huge maxsize, this skips the test for too large
 * bricks etc. which means I can use the tool to find holes in the
 * allocated data.
 *
 * Repeats the clipping against uncompressed brick/tile size
 * (no compressed data is allowed to be larger than that)
 * because maxsize might not have been available at that point.
 * Could also have repeated the test for eof, but that should
 * always have been available.
 */
std::pair<std::int64_t, std::int64_t>
LookupTable::_getBegAndSize(
    std::int64_t ix,
    const std::vector<std::uint64_t>& lup,
    const std::vector<std::uint64_t>& end,
    std::int64_t maxsize)
{
  // TODO-Low consider deferring _calc_lookupsize() until first needed.
  // The simple version below id NOT THREADSAFE.
  // if self._metadata._blup._lookend is None:
  //   self._metadata._blup._lookend = (
  //        self._metadata._blup._calc_lookupsize(self._lookup, eof, maxsize))

  const std::uint64_t raw_beg = lup[ix];
  const std::uint64_t raw_end = end[ix];
  const std::uint8_t type = static_cast<uint8_t>((raw_beg >> 56) & 0xff);
  const std::uint64_t addrmask = 0x00ffffffffffffff;
  // TODO-Low static_cast<std::int64_t>(~(static_cast<std::uint64_t>(0xFF) << 56));

  if (type == 0x80 || raw_beg < 2)
    return std::make_pair(0, 0);
  else if (type == 0 || type == 0xC0) {
    const std::uint64_t beg = raw_beg & addrmask;
    const std::uint64_t end = std::max(beg, raw_end);
    return std::make_pair(beg, std::min(end - beg, (std::uint64_t)maxsize));
  }
  else
    throw ZgyFormatError("Unknown alpha- or brick type " + std::to_string(type));
}

/**
 * \brief Get file offset for the specified alpha tile. Input i/j/lod.
 *
 * \details
 * Convert directly from i/j/lod to the entry for this tile.
 * The entry gives the type, file offset, constant value, etc.
 * Using this method hides knowledge of the interal "index" space
 * as well as how the entry is encoded into a single 64-bit integer.
 */
LookupTable::LutInfo
LookupTable::getAlphaFilePosition(
    std::int64_t i, std::int64_t j, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& alphaoffsets,
    const std::vector<std::uint64_t>& alup,
    std::int64_t bytesperalpha)
{
  const std::int64_t index =
    getAlphaLookupIndex(i, j, lod, lodsizes, alphaoffsets);
  return getAlphaFilePositionFromIndex(index, alup, bytesperalpha);
}

/**
 * \brief Get file offset for the specified alpha tile. Input linear index,
 *
 * \details
 * Using this method implies knowledge of how (i,j,k,lod) maps to
 * an internal "index" position. It might be better to use the
 * corresponding method without "Index" in its name.
 */

LookupTable::LutInfo
LookupTable::getAlphaFilePositionFromIndex(
    std::int64_t index,
    const std::vector<std::uint64_t>& alup,
    std::int64_t bytesperalpha)
{
  const std::uint64_t pos = alup[index];
  const std::uint8_t type = static_cast<uint8_t>((pos >> 56) & 0xff);
  if (pos == 0)
    return LutInfo(BrickStatus::Missing, 0, 0, 0);
  else if (pos == 1)
    return LutInfo(BrickStatus::Constant, 0, 0, 0);
  else if (type == 0x80)
    return LutInfo(BrickStatus::Constant, 0, 0, pos & 0xff);
  else if (type == 0xC0)
    throw ZgyInternalError("Compressed alpha tiles not yet implemented");
  else if (type & 0x80)
    throw ZgyFormatError("Unknown alpha type " + std::to_string(type));
  else
    return LutInfo(BrickStatus::Normal, pos, bytesperalpha, 0);
}

/**
 * \brief Get file offset or constant-value for the specified brick.
 * Input i/j/k/lod.
 *
 * \details
 * Convert directly from i/j/k/lod to the entry for this brick.
 * The entry gives the type, file offset, constant value, etc.
 * Using this method hides knowledge of the interal "index" space
 * as well as how the entry is encoded into a single 64-bit integer.
 */
LookupTable::LutInfo
LookupTable::getBrickFilePosition(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    const std::vector<std::uint64_t>& blup,
    const std::vector<std::uint64_t>& bend,
    std::int64_t bytesperbrick)
{
  const std::int64_t ix =
    getBrickLookupIndex(i, j, k, lod, lodsizes, brickoffsets);
  return getBrickFilePositionFromIndex(ix, blup, bend, bytesperbrick);
}

/**
 * \brief Get file offset or constant-value for the specified brick.
 * Input linear index.
 *
 * \details
 * Using this method implies knowledge of how (i,j,k,lod) maps to
 * an internal "index" position. It might be better to use the
 * corresponding method without "Index" in its name.
 */
LookupTable::LutInfo
LookupTable::getBrickFilePositionFromIndex(
    std::int64_t ix,
    const std::vector<std::uint64_t>& blup,
    const std::vector<std::uint64_t>& bend,
    std::int64_t bytesperbrick)
{
  const std::int64_t pos = blup[ix];
  const std::uint8_t type = static_cast<uint8_t>((pos >> 56) & 0xff);
  if (pos == 0)
    return LutInfo(BrickStatus::Missing, 0, 0, 0);
  else if (pos == 1)
    return LutInfo(BrickStatus::Constant, 0, 0, 0);
  else if (type == 0x80)
    return LutInfo(BrickStatus::Constant, 0, 0,
                   static_cast<std::uint32_t>(pos & 0xffffffff));
  else if (type == 0xC0) {
    std::pair<std::int64_t, std::int64_t> pair =
      _getBegAndSize(ix, blup, bend, bytesperbrick);
    return LutInfo(BrickStatus::Compressed, pair.first, pair.second, 0);
  }
  else if (type & 0x80)
    throw ZgyFormatError("Unknown brick type " + std::to_string(type));
  else
    return LutInfo(BrickStatus::Normal, pos, bytesperbrick, 0);
}

/**
 * \brief Set file offset or constant-value for the specified brick.
 */
void
LookupTable::setBrickFilePosition(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const LutInfo& info,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    std::vector<std::uint64_t>* blup,
    std::vector<std::uint64_t>* bend)
{
  const std::int64_t addrmask = 0x00ffffffffffffff;
  const std::int64_t ix =
    getBrickLookupIndex(i, j, k, lod, lodsizes, brickoffsets);
  switch (info.status) {
  case BrickStatus::Missing:
    (*blup)[ix] = 0;
    (*bend)[ix] = 0;
    break;
  case BrickStatus::Constant:
    (*blup)[ix] = (static_cast<std::uint64_t>(1)<<63) | (info.raw_constant & 0xffffffff);
    (*bend)[ix] = 0;
    break;
  case BrickStatus::Normal:
    (*blup)[ix] = info.offset_in_file & 0x7fffffffffffffff;
    (*bend)[ix] = (*blup)[ix] + info.size_in_file;
    break;
  case BrickStatus::Compressed:
    (*blup)[ix] = ((static_cast<std::uint64_t>(0xC0) << 56) |
                   (info.offset_in_file & addrmask));
    (*bend)[ix] = info.offset_in_file + info.size_in_file;
    break;
  default:
    throw ZgyFormatError("Unknown brick enum " + std::to_string((int)info.status));
  }
}

/**
 * \brief Check whether finalize() has written low resolution data.
 * \returns Number of levels with data. Including the full resolution level.
 *
 * \details
 * If low resolution data has been generated then all its bricks should
 * exist. Brick status might be "Constant" but should bever be "Missing".
 *
 * Only test the single brick at the highest LOD level. Other LOD levels
 * are by convention treated as empty when the highest level is empty.
 * This is to allow re-use of disk space. See setNoBrickLOD.
 * Level 0 (full resolution) is by convention always assumed to exist
 * even if all its bricks are "Missing".
 *
 * After the next format update there might be a specific bit pattern
 * in the lookup table for "logically deleted but here is the disk area
 * that used to be allocated".
 */
std::int32_t
LookupTable::usableBrickLOD(
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    const std::vector<std::uint64_t>& blup,
    const std::vector<std::uint64_t>& bend)
{
  //return lodsizes.size() <= 1 || blup[0]==0 ? 1 : lodsizes.size(); // Kludge
  return (lodsizes.size() <= 1) ? 1 :
    LookupTable::getBrickFilePosition
    (0, 0, 0, lodsizes.size()-1,
     lodsizes, brickoffsets, blup, bend,
     0).status == BrickStatus::Missing ? 1 :
    (std::int32_t)lodsizes.size();
}

bool
LookupTable::hasBrickCompression(
      const std::vector<std::uint64_t>& blup,
      const std::vector<std::uint64_t>& bend)
{
  constexpr std::int64_t bytesperbrick{1}; // Bogus, we only care about type.
  for (std::size_t ix = 0; ix < blup.size(); ++ix) {
    LookupTable::LutInfo info =
      LookupTable::getBrickFilePositionFromIndex(ix, blup, bend, bytesperbrick);
    if (info.status == BrickStatus::Compressed)
      return true;
  }
  return false;
}

/**
 * \brief Logically erase all low resolution bricks.
 *
 * \details
 * It is sufficient to clear the highest lod level. Which is always
 * just one brick. And it just happens to always be the first entry in
 * the lookup table. It is possible but too ugly to use a huge
 * shortcut.
 *
 * Note that the disk space for the single brick being cleared will be
 * leaked. Other LOD levels are by convention treated as empty when the
 * highest level is empty. But the file pointers remain. So for an
 * uncompressed on-prem file the disk space for all but one entry
 * will be re-used the next time finalize() is run.
 *
 * If the entire file is just a singe brick there is low resolution.
 */
void
LookupTable::setNoBrickLOD(
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    std::vector<std::uint64_t>* blup,
    std::vector<std::uint64_t>* bend)
{
  if (lodsizes.size() > 1) {
    // (*blup)[0] = 0; // Too much of a kludge.
    LookupTable::setBrickFilePosition
      (0, 0, 0, lodsizes.size()-1,
       LookupTable::LutInfo(BrickStatus::Missing, 0, 0, 0),
       lodsizes, brickoffsets, blup, bend);
  }
}

} // namespace
