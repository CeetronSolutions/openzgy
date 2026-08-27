// Copyright 2017-2024, Schlumberger
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

/**
 * \file sanity.cpp
 *
 * Methods to do sanity checks on brick- and survey size and data type
 * both when creating a new file (because the application might try to
 * do something silly) and when opening an existing file (because the
 * file might have been corrupted, either maliciously or by accident).
 * Throw on error.
 *
 * Corrupted settings that only cause incorrect results, not crashes
 * or confusing exceptions at a later time, need to be checked here.
 *
 * The checks should be run as early as possible. Before any bad
 * parameters can cause trouble.
 *
 * Other checks might be added later. It is tidier to have all the
 * checks in one place, here in class SanityCheck. But maybe not
 * practical.
 *
 * Rules:
 *
 *   - Data type must be one of the supported (int8, int16, float).
 *
 *   - Brick- and survey size cannot be zero or negative.
 *
 *   - Brick- and survey size in each dimension must not be so large
 *     that computing the total size in bytes will cause a numerical
 *     overflow.
 *
 *   - Brick size in each dimension must be a power of 2.
 *
 *   - Attempting to read past EOF is an error. Checking this is
 *     currently done in FileADT::_validate_read instead of up front.
 *     That might be sufficient
 *
 *     Not checked up front: Survey size, which indirectly controls
 *     the size of both lookup tables, must not be so large that the
 *     lookup tables become ridiculously big. The test for reading
 *     past EOF might be good enough. Earlier checking might be
 *     possible in OffsetHeaderV2Access::calculate(). That can be done
 *     earlier but still doesn't guarantee that it is early enough.
 *
 *   - Brick- and survey sizes should not be so large that they cause
 *     some function to run out of memory. Or cause the application to
 *     run so slowly that it appears to be hung. This is a fuzzy
 *     requirement. Currently limit the total brick size in bytes to 2
 *     GB. Which also limits the total brick size to fit in an int32.
 *     Do not limit survey size except as specified above. If it turns
 *     out that we do run out of memory, just hope that the code gets
 *     an orderly std::bad_alloc later.
 *
 *   - Silently fixed: vunitfactor and hunitfactor must be positive
 *     to avoid a divide by zero. Check both on create and open.
 *     Given that most applications ignore those settings amyway,
 *     it might be best to quietly replace bad values with 1.
 *     Handled in InfoHeaderV2Access::{h,v}unitfactor() for read
 *     and in ZgyInternalMeta::validateCreateNewFile()
 *     and ZgyInternalMeta::validateOpenForUpdate() for write.
 *
 *   - Could be more robust: The string header size, which only
 *     applies to opening an existing file version two or later, must
 *     not be negative or ridiculously large. Something that causes
 *     any of the headers to extend past the current end of file is
 *     clearly wrong. Check on open only. The code might handle this
 *     already, by the "is past EOF" check, also for "negative" size,
 *     because the size as stored on file is declared as uint32. If
 *     the size is widened to int64 or before being used (to prevent
 *     overflow) and if "is past EOF" is done early enough then we are
 *     good. An explicit, earlier check might be more robust.
 *
 *   - The string header, which only applies to opening an existing
 *     file version two or later, should contain at least 5 strings.
 *     Silently treat missing strings as empty. Silently add a
 *     terminating null if the buffer doesn't have one. Checked in
 *     InfoHeaderV2Access::calculate_read().
 *
 *   - No brick or alpha pointer can start at eof or later. Currently
 *     LookupTable::calcLookupSize() checks this. That is probably
 *     early enough.
 *
 *   - No bricks should overlap unless doing some strictly local
 *     debugging. Currently the LookupTable::calcLookupSize() function
 *     checks this. Temporarily remove this check, strictly for
 *     debugging, to allow manually crafting a huge survey but with a
 *     reasonable size on disk by storing many pointers to each
 *     physical brick. This would be to do certain stress tests that
 *     might be difficult otherwise.
 */

#include "sanity.h"
#include "enum.h"
#include "types.h"
#include "environment.h"
#include "../exception.h"

#include <string>
#include <sstream>
#include <limits>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Check that a single number is positive and a power of two.
 * Currently this is a requirement for the ZGY brick size.
 * Any other size is untested and will probably not work.
 * The code is not efficient, but it is portable and really
 * not used much.
 */
bool
SanityCheck::isPowerOfTwo(std::int64_t n)
{
  for (int shift=0; shift<63; ++shift)
    if (n == static_cast<std::int64_t>(1) << shift)
      return true;
  return false;
}

/**
 * Return true if these unsigned inputs cannot be safely multiplied
 * together. The maximum result size is am argument and should not
 * be more than can fit in an unsigned integer. It could be less.
 * E.g. if the caller suspects that the result might be cast to a
 * signed integer at some point.
 */
bool
SanityCheck::willUnsignedMultiplyOverflow(
     std::uint64_t a,
     std::uint64_t b,
     std::uint64_t lim)
{
  return (a != 0 && b != 0 && (a > lim / b || b > lim / a));
}

/**
 * Return true if these signed inputs cannot be safely multiplied
 * together. Also return true if not all inputs are greater than zero.
 */
bool
SanityCheck::isNotPositiveOrWillMultiplyOverflow(
     std::int64_t a,
     std::int64_t b,
     std::uint64_t lim)
{
  return (a <= 0 || b <= 0 ||
          willUnsignedMultiplyOverflow
          (static_cast<std::uint64_t>(a),
           static_cast<std::uint64_t>(b),
           lim));
}

/**
 * Return true if these signed inputs cannot be safely multiplied
 * together. Also return true if not all inputs are greater than zero.
 * Could have been templated on array size, but only dims=3 is needed.
 */
bool
SanityCheck::isNotPositiveOrWillMultiplyOverflow(
     const std::array<std::int64_t,3>& in,
     std::int64_t bytes_per_sample,
     std::uint64_t lim)
{
  return (isNotPositiveOrWillMultiplyOverflow(in[0], in[1], lim) ||
          isNotPositiveOrWillMultiplyOverflow(in[0] * in[1], in[2], lim) ||
          isNotPositiveOrWillMultiplyOverflow(in[0] * in[1] * in[2], bytes_per_sample, lim));
}

/**
 * Throw if this brick size is not acceptable.
 *
 * Note that even though the brick size is ussually passed as
 * int64[3], the code is unlikely to work if the product won't fit in
 * an int32. This is enforced.
 *
 * The actual limit might be lower due to memory constraints. This is
 * not enforced because it is too difficult to figure out what that
 * limit should be.
 *
 * Officially there is no support for 2d data and library probably
 * doesn't work so a size < 4 in any direction should only be used for
 * experiments. Setting OPENZGY_ALLOW_2D removes this particular
 * restriction.
 */
void
SanityCheck::checkValidBrickSize(
     const std::array<std::int64_t,3>& in,
     const std::int64_t bytes_per_sample)
{
  static std::int64_t minsize =
    Environment::getNumericEnv("OPENZGY_ALLOW_2D", 0) > 0 ? 1 : 4;
  constexpr std::uint64_t maxsize_total =
    static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max());
  std::stringstream err;
  for (std::size_t ii=0; ii<in.size(); ++ii) {
    if (in[ii] < minsize || !isPowerOfTwo(in[ii])) {
      err << "bricksize must be >= " << minsize << " and a power of 2.";
      throw OpenZGY::Errors::ZgyUserError(err.str());
    }
  }
  if (isNotPositiveOrWillMultiplyOverflow(in, bytes_per_sample, maxsize_total)) {
    err << "Total brick size in bytes must be <= 0x"
        << std::hex << maxsize_total << std::dec;
    throw OpenZGY::Errors::ZgyUserError(err.str());
  }
}

/**
 * Throw if this survey size is not acceptable.
 *
 * The size in each dimension needs to fit in an int32. This is a
 * limitation in the physical file format. If the format is changed,
 * the tests here may need to be relaxed. As well as carefully
 * auditing and testing that the rest of the library can handle that.
 *
 * The actual limit might be lower due to memory constraints. This is
 * not enforced because it is too difficult to figure out what that
 * limit should be.
 */
void
SanityCheck::checkValidSurveySize(
     const std::array<std::int64_t,3>& in,
     const std::int64_t bytes_per_sample)
{
  constexpr std::uint64_t maxsize_total =
    static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
  constexpr std::uint64_t maxsize_dim =
    static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max());
  std::stringstream err;
  for (std::size_t ii=0; ii<in.size(); ++ii) {
    if (in[ii] < 1 || static_cast<std::uint64_t>(in[ii]) > maxsize_dim) {
      err << "survey size[" << ii << "] is " << in[ii]
          << " but must be between 1 and 0x" << std::hex << maxsize_dim
          << std::dec;
      throw OpenZGY::Errors::ZgyUserError(err.str());
    }
  }
  if (isNotPositiveOrWillMultiplyOverflow(in, bytes_per_sample, maxsize_total)) {
    err << "Total survey size in bytes must be <= 0x"
        << std::hex << maxsize_total << std::dec;
    throw OpenZGY::Errors::ZgyUserError(err.str());
  }
}

/**
 * Do sanity checks on brick- and survey size and data type both when
 * creating a new file (because the application might try to do
 * something silly) and when opening an existing file (because the
 * file might have been corrupted, either maliciously or by accident).
 * Throw on error.
 *
 * See the sanity.cpp file documentation for more details.
 *
 * This function should be run as early as possible. Before any bad
 * parameters can cause trouble.
 *
 * For create, the sizes are found in bricksize_, size_, datatype_
 * in ZgyWriterArgsV3Impl. For open existing files, use bricksize(),
 * size(), datatype() in IInfoHeaderAccess
 */
void
SanityCheck::checkValidBrickAndSurveySize(
     const std::array<std::int64_t,3>& bricksize,
     const std::array<std::int64_t,3>& size,
     RawDataType dtype)
{
  switch (dtype) {
  case RawDataType::SignedInt8:
  case RawDataType::SignedInt16:
  case RawDataType::Float32:
    break;
  default:
    throw OpenZGY::Errors::ZgyUserError("Datatype must be int8, int16, or float");
  }
  const std::int64_t bytes_per_sample = RawDataTypeDetails(dtype).size;
  checkValidBrickSize(bricksize, bytes_per_sample);
  checkValidSurveySize(size, bytes_per_sample);
}

} // namespace
