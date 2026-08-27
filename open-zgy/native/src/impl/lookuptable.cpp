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
#include "../exception.h"
#include "enum.h"

#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <array>
#include <cstdint>
#include <limits>
#include <iostream>

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

  if (!entries.empty())
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
 * Or silently return false if caller asked for that.
 * If used to validate an alpha tile then k should be passed as 0.
 */
bool
LookupTable::_validatePosition(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    bool throw_on_error)
{
  if (lod < 0 || i<0 || j<0 || k<0) {
    if (throw_on_error)
      throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                         " cannot be negative");
    else
      return false;
  }
  const std::int64_t nlods = static_cast<std::int64_t>(lodsizes.size());
  if (lod >= nlods) {
    if (throw_on_error)
      throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                         " cannot have lod >= " + std::to_string(nlods));
    else
      return false;
  }
  const std::array<std::int64_t,3>& size = lodsizes[lod];
  if (i >= size[0] || j >= size[1] || k >= size[2]) {
    if (throw_on_error)
      throw ZgyUserError("Requested brick " + _formatPosition(i, j, k, lod) +
                         " cannot be >= " +
                         _formatPosition(size[0], size[1], size[2], nlods-1));
    else
      return false;
  }
  return true;
}

/**
 * If throw_on_error = true, the function guarantees that the returned index
 * is between 0 and offsets.back(). If throw_on_error = false then the
 * return can additionally be negative, but still not larger than the
 * last valid brick.
 */
std::int64_t
LookupTable::_getLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& offsets,
    bool throw_on_error)
{
  if (!_validatePosition(i, j, k, lod, lodsizes, throw_on_error))
    return -1;
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
    const std::vector<std::int64_t>& alphaoffsets,
    bool throw_on_error)
{
  return _getLookupIndex(i, j, 0, lod, lodsizes, alphaoffsets, throw_on_error);
}

/**
 * Normally called only from inside LookupTable, but can be useful
 * for dealing with other per brick information such as dirty flags.
 */
std::int64_t
LookupTable::getBrickLookupIndex(
    std::int64_t i, std::int64_t j, std::int64_t k, std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    bool throw_on_error)
{
  return _getLookupIndex(i, j, k, lod, lodsizes, brickoffsets, throw_on_error);
}

/**
 * Get all indices for the parents (i.e. lod-1) of the specified bricks.
 * Input indices are block numbers relative to the given lod, NOT lod-1.
 *
 * If throw_on_error = true, all input bricks need to be within the
 * legal range. And the lod cannot be 0 (where there are no parents).
 * If throw_on_error = false, outside is simply ignored.
 *
 * In both cases it is not an error for the input brick(s) to have
 * less than 8 parents each. Due to being close to the survey edge.
 */
std::vector<std::int64_t>
LookupTable::getParentBrickIndices(
    std::int64_t i0, std::int64_t j0, std::int64_t k0,
    std::int64_t ni, std::int64_t nj, std::int64_t nk,
    std::int64_t lod,
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    bool throw_on_error)
{
  if (throw_on_error) {
    _validatePosition(i0, j0, k0, lod, lodsizes, /*throw=*/true);
    _validatePosition(i0*2, j0*2, k0*2, lod-1, lodsizes, /*throw=*/true);
    _validatePosition(i0+ni-1, j0+nj-1, k0+nk-1, lod, lodsizes, /*throw=*/true);
  }
  std::vector<std::int64_t> result;
  result.reserve(ni*nj*nk*8);
  std::int64_t idx;
  for (std::int64_t subi = i0 * 2; subi < (i0 + ni) * 2; ++subi)
    for (std::int64_t subj = j0 * 2; subj < (j0 + nj) * 2; ++subj)
      for (std::int64_t subk = k0 * 2; subk < (k0 + nk) * 2; ++subk)
        if ((idx = _getLookupIndex
             (subi, subj, subk, lod-1, lodsizes, brickoffsets, false)) >= 0)
          result.push_back(idx);
  return result;
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
    const std::uint64_t my_beg = raw_beg & addrmask;
    const std::uint64_t my_end = std::max(my_beg, raw_end);
    return std::make_pair(my_beg, std::min(my_end - my_beg, (std::uint64_t)maxsize));
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
    getAlphaLookupIndex(i, j, lod, lodsizes, alphaoffsets, /*throw=*/true);
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
    getBrickLookupIndex(i, j, k, lod, lodsizes, brickoffsets, /*throw=*/true);
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
    getBrickLookupIndex(i, j, k, lod, lodsizes, brickoffsets, /*throw=*/true);
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
 * \brief Forward to LookupTable::getFinalizedStatus()
 */
bool
LookupTable::hasBrickLOD(
    const std::vector<std::array<std::int64_t,3>>& lodsizes,
    const std::vector<std::int64_t>& brickoffsets,
    const std::vector<std::uint64_t>& blup,
    const std::vector<std::uint64_t>& bend)
{
  const LookupTable::FinalizedStatus status = LookupTable::getFinalizedStatus
    (lodsizes, brickoffsets, blup, bend);

  switch (status) {
  case LookupTable::FinalizedStatus::SingleBrick:
  case LookupTable::FinalizedStatus::ConstantFile:
  case LookupTable::FinalizedStatus::IsFinalized:
    return true;

  default:
  case LookupTable::FinalizedStatus::NotFinalized:
    return false;
  }
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
 * \brief Check whether the file is known to be all constant.
 *
 * \details
 *
 * Unlike doing a readconst() on the entire survey, this function also
 * checks the low resolution bricks. Unlike readconst() it will not
 * return the actual sample value.
 *
 * There are a few unlikely cases where the function returns false in
 * spite of the file being all constant. The opposite case will not
 * happen; a true return can be trusted.
 *
 * \internal
 *
 * See also filestats_nocache(), it does much of the same.
 * See also ZgyInternalBulk::readConstantValue().
 *
 * Implementation: Count the number of tiles or bricks between
 * first_lod and last_lod (inclusive) that have bulk data allocated on
 * the file. Also do a superficial check of whether all the lod(s)
 * have a single value.
 *
 *  - File has never been written to: Always return true.
 *
 *  - File has been written exactly once, setting the entire survey
 *    to a constant value: Always return true.
 *
 *  - Other all constant surveys *might* return true, see below.
 *
 *  - Part of the file is written with a constant value, part has
 *    never been written, and the constant value happens to match the
 *    default value. TODO-WIP-BrickedAPI: BUG: Returns a false negative.
 *    The problem: For float data and for unscaled data, 0, 1, and 0x80xx
 *    all mean the same thing: constant zero. But for scaled integer data
 *    0 (missing) might be different from 1 and 0x80xx (storage zero).
 *    And this code is at a too low level to get the correct float value.
 *    To avoid too surprising behavior, I won't fix this for float either.
 *
 *    This can happen e.g. when finalizing a file that was never
 *    written to, using LodMode::Rebuild. The file ends up treated as
 *    "not finalized". Which is wrong. And very counter intuitive
 *    since the application wanted to make really sure the file git
 *    finalized. Ironically, passing LodMode::Never does the right
 *    thing in this case. See unit test "reopen.setversion"
 *
 *    ZgyInternalBulk::readConstantValue() would have handled this
 *    because it actually decodes the raw constant.
 *
 *  - The file has at one point contained real data but was later
 *    cleared. Return a false negative, because all allocated bricks
 *    are assumed to still be non-constant.
 *
 *    TODO-WIP-BrickedAPI: This can cause a user visible problem. If
 *    Finalize wasn't run until lod0 had been reset to become all
 *    const, the file would look like it had data but no finalize. It
 *    would not be allowed to update this file with LodMode::Early. If
 *    we had done a thorough check by reading all real bricks then the
 *    problem would go away. At a ridiculous performance penalty.
 *
 *    Unit test "reopen.track_changes" spotted this issue.
 *    Issue also noted in ZgyWriter::checkValidLodMode()
 */
bool /*static*/
LookupTable::hasAllBricksSameValue(
     int first_lod, int last_lod,
     const std::vector<std::int64_t>& offsets,
     const std::vector<std::uint64_t>& lup,
     const std::vector<std::uint64_t>& lupend)
{
  constexpr std::int64_t bytesperbrick{1}; // Bogus, we only care about type.

  std::int64_t missing_count{0};
  std::int64_t missing_count_fullres{0};
  std::int64_t constant_count{0};
  std::int64_t normal_count{0};
  std::int64_t compressed_count{0};
  std::uint32_t raw_const_value{0};
  bool is_all_same{true};

  // Note that higher lods have lower offsets. Lod nlods-1 has offest 0.
  // The end of lod0 is a special case. It starts at offsets[0] and
  // ends at offsets.size()-1, because the total size of the lookup
  // table has been apended to the end of offsets as a convenience.
  const std::int64_t nlods = offsets.size() - 1;
  const std::int64_t lo_ix =
    (last_lod < 0 || last_lod >= nlods) ? 0 : offsets[last_lod];
  const std::int64_t hi_ix =
    (first_lod <= 0) ? offsets[nlods] : // end (not inclusive) lod 0
    (first_lod >= nlods) ? 0 : // result will be zero
    offsets[first_lod-1];
  const std::int64_t lo_lod0_ix = offsets[0];
  const std::int64_t hi_lod0_ix = offsets[nlods];

  //std::cout << "Search in lut " << lo_ix << " to " << hi_ix
  //          << " (max " << lup.size() << ")"
  //          << " for lod " << first_lod << " to " << last_lod << " inclusive"
  //          << std::endl;
  if (lo_ix < 0 || hi_ix < 0 || hi_ix > offsets[nlods])
    throw ZgyInternalError("hasAllBricksSameValue index out of range");

  for (std::int64_t ix = lo_ix; ix < hi_ix; ++ix) {
    const LookupTable::LutInfo info =
      LookupTable::getBrickFilePositionFromIndex(ix, lup, lupend,bytesperbrick);
    switch (info.status) {
    case BrickStatus::Missing:
      ++missing_count;
      if (ix >= lo_lod0_ix && ix < hi_lod0_ix)
         ++missing_count_fullres;
      if (constant_count != 0)
        is_all_same = false;
      break;
    case BrickStatus::Constant:
      if (constant_count == 0) {
        raw_const_value = info.raw_constant;
      }
      else {
        if (raw_const_value != info.raw_constant)
          is_all_same = false;
      }
      if (missing_count)
        is_all_same = false;
      ++constant_count;
      break;
    case BrickStatus::Normal:
      ++normal_count;
      is_all_same = false;
      break;
    case BrickStatus::Compressed:
      ++compressed_count;
      is_all_same = false;
      break;
    default:
      break;
    }
  }

  // Mitigation for the case where the file has never been written,
  // not even with a constant value, but lowres has been generated.
  // Since bricks cannot change back to Missing (except for the
  // last lowres brick, which is another kludge) we can assume the
  // lowres is not stale.
  if (!is_all_same) {
    if (first_lod == 0 && missing_count_fullres == hi_lod0_ix - lo_lod0_ix) {
      //std::cout << "@@ Fuzzy: This file is const because all "
      //          << missing_count_fullres
      //          << " lod0 bricks are missing."
      //          << std::endl;
      is_all_same = true;
    }
  }

  return is_all_same;
}

/**
 * Fuzzy logic to tell us whether the file contains lowres data or not,
 * i.e. whether the last write used LodMode::Never or not. This is
 * surprisingly different until the next ZGY file format update,
 * where there should be an explicit field for "usable lods".
 *
 * The heuristic is not exact. Case in point: hasAllBricksSameValue()
 * doesn't detect a brick that was first written with real data, then
 * overwritten with a constant. That leaves the brick inflated, and it
 * is way too expensive to read and check for that. I doubt that this
 * will happen during normal usage. Similarly, a mix of Constant and
 * Missing bricks will currently be seen as "not all the same".
 *
 * Another obscure case is when lod0 contains real data, but the
 * decimation algorithm made all the lowres, even lod1, a constant.
 * Technically we could allow both reading all lowres and allow
 * re-open with any LodMode. Just like the ConstantFile case.
 * The case is probably too obscure to worry about.
 *
 * One corner case needs to be handled, but this needs to be done in
 * _close_internal() and not here. If the top lodN is Constant, the
 * code assumes the file is finalized with some aggressive algorithm
 * that removed all samples. But the survey might also have been
 * cleared, setting the top lodN, and then written fullres without
 * finalizing. The code in _close_internal() must then switch the top
 * lod0 (or perhaps all lowres) to become Missing *if* the file was
 * actually openen with LodMode::Never. This is the last point where
 * we know for sure how the file was written.
 *
 * After the next format update there might be a specific bit pattern
 * in the lookup table for "logically deleted but here is the disk area
 * that used to be allocated". If so, updates may be needed here.
 *
 * LOD levels 1..nlods-2 are by convention treated as missing when the
 * highest level is empty, even when actual data is allocated. This is
 * to allow re-use of disk space. Currently, the application isn't
 * allowed to clear lowres by switching an already finalized file to
 * LodMode::Never. So the disk space issue is N/A. But clearing the
 * top level with setNoBrickLOD() might still be needed.
 *
 * Consider the following quite reasonable case:
 * - Application creates a new file with LodMode::Never.
 * - writeconst() the default sample value to the entire survey.
 * - write() real data without generating lowres.
 *
 * A subsequent open for read will believe that lowres exists. It is
 * valid for decimation of non-const data end up being const. To
 * prevent the wrong conclusion, the top lowres brick, or maybe all
 * lowres bricks, need to be changed from Constant to Missing.
 * That needs to be done when the file is closed. Because that is
 * the last place we know that user wrote with LodMode::Never.
 *
 * Returns:
 *   - NotFinalized   => Only lod 0 usable. Next update needs Never or Rebuild.
 *   - IsFinalized    => All lods usable. Next update cannot use Never.
 *   - SingleBrick    => Single brick survey. Only lod 0. LodMode irrelevant.
 *   - ConstantFile   => All lods usable. Next update can use any mode.
 *
 * Called by checkValidLodMode when updating a file, and ZgyMeta::nlods() via
 * LookupTable::hasBrickLOD() when deciding how many lod levels can be read.
 */
LookupTable::FinalizedStatus
LookupTable::getFinalizedStatus(
     const std::vector<std::array<std::int64_t,3>>& lodsizes,
     const std::vector<std::int64_t>& brickoffsets,
     const std::vector<std::uint64_t>& blup,
     const std::vector<std::uint64_t>& bend)
{
  const std::int32_t nlods = static_cast<std::int32_t>(lodsizes.size());
  if (nlods == 1)
    return LookupTable::FinalizedStatus::SingleBrick;

  // Get brick (0, 0, 0, last_lod)
  const BrickStatus top_status =
    LookupTable::getBrickFilePosition
    (0,0,0, lodsizes.size()-1, lodsizes, brickoffsets, blup, bend, 0).status;

  // The all-same test is expensive. And if the top brick has real data
  // then it can never succeed. The call is still a performance problem
  // when reading from a file that doesn't have lowres (argument check),
  // The only code that bypasses hasBrickLOD and calls getFinalizedStatus()
  // directly is checkValidLodMode() which is only called once per file
  // open. The result of calling hasBrickLOD() ought to be cached in some
  // cases. This must be done by the caller, sincd LookupTable is a static
  // class.
  if (top_status!=BrickStatus::Normal && top_status!=BrickStatus::Compressed) {
    const bool all_const = InternalZGY::LookupTable::hasAllBricksSameValue
      (0, -1, brickoffsets, blup, bend);
    if (all_const)
      return LookupTable::FinalizedStatus::ConstantFile;
  }

  switch (top_status) {
  default: // should not happen
  case BrickStatus::Missing:
    return LookupTable::FinalizedStatus::NotFinalized;
  case BrickStatus::Constant:
  case BrickStatus::Normal:
  case BrickStatus::Compressed:
    return LookupTable::FinalizedStatus::IsFinalized;
  }
}

/**
 * \brief Count the number of bricks that are physically allocated.
 *
 * \details
 *
 * Can be used to allow or disallow opening for update with compression.
 * Similar to hasAllBricksSameValue(), but here we aren't worried about
 * contents. Only about bricks (typically lowres) where we might end
 * up overwriting an already allocated brick.
 */
std::int64_t /*static*/
LookupTable::countAllocatedBricks(
     const std::vector<std::int64_t>& offsets,
     const std::vector<std::uint64_t>& lup,
     const std::vector<std::uint64_t>& lupend)
{
  std::int64_t count{0};
  for (std::int64_t ix = 0; ix < offsets.back(); ++ix) {
    const LookupTable::LutInfo info =
      LookupTable::getBrickFilePositionFromIndex(ix, lup, lupend, 1);
    switch (info.status) {
    case BrickStatus::Normal:
    case BrickStatus::Compressed:
      ++count;
      break;
    default:
      break;
    }
  }
  return count;
}

/**
 * \brief Logically erase all low resolution bricks.
 *
 * \details
 *
 * ONLY use this for the special case in _close_internal(), to prevent
 * readers from believing that lowres is present because of a clear
 * survey. In this case there OpenZGY no longer allows deleting the
 * lowres of an existing finalized file. That would mess up the
 * single-pass logic. The _close_internal code, in contrast, is only
 * triggered when the file doesn't have lowres, it just looks that
 * way.
 *
 * It is sufficient to clear the highest lod level. Which is always
 * just one brick. And it just happens to always be the first entry in
 * the lookup table. It is possible but too ugly to use a huge
 * shortcut.
 *
 * The following is N/A when used in the special case.
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
