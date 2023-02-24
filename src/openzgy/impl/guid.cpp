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

#include "guid.h"

#include <random>
#include <chrono>
#include <mutex>
#include <string.h>

namespace InternalZGY {
#if 0
}
#endif

GUID::GUID() : _data(generate())
{
}

GUID::GUID(const guid_bytes_t& in) : _data(in)
{
}

GUID::GUID(std::nullptr_t) : _data{0}
{
}

std::string
GUID::toString() const
{
  return format(this->_data);
}

std::string
GUID::makeGUID()
{
  return GUID(nullptr).toString();
}

/**
 * Copy the raw bytes of the guid to the specified pointer.
 * The length is redundant. It should always be 16 but is
 * included because it is then more obvious at the calling
 * site what is going on.
 */
void
GUID::copyTo(std::uint8_t *ptr, std::int64_t len)
{
  memcpy(ptr, _data.data(), std::min((std::int64_t)sizeof(_data), len));
}

/**
 * Generate a random GUID/UUID, according to the instructions at
 * https://en.wikipedia.org/wiki/Universally_unique_identifier.
 *
 * The result is returned as an array of bytes but *not* in
 * big-endian order as specified by RFC 4122. Instead use the
 * piecewise little-endian format that would allow casting the
 * buffer to a struct UUID on a little-endian nachine.
 *
 * See doc/uuid.md.
 *
 * Implementing general support for all GUID versions and variants is
 * ridiculously complicated. So, go the easy route and generate a
 * random one. Then the only trick is to tell any reader that the guid
 * is stored big-endian (variant 1) and is in fact a random number
 * (version 4).
 *
 * TODO-Worry: Is the entropy of the random seed good enough?
 */
GUID::guid_bytes_t
GUID::generate()
{
  static std::mutex mutex;
  std::lock_guard<std::mutex> lk(mutex);
  guid_bytes_t result;
  typedef std::default_random_engine engine_t;
  static engine_t generator
    {static_cast<engine_t::result_type>
     (std::chrono::system_clock::now().time_since_epoch().count())};
  static std::uniform_int_distribution<int> distribution(0,255);
  for (std::size_t ii = 0; ii < result.size(); ++ii){
    int one = distribution(generator);
    result[ii] = static_cast<std::uint8_t>(one);
  }

  // The variant is encoded into byte 8 i.e. the non-byteswapped part.
  // The version is encoded into the highest nibble in the int16 field
  // in bytes 6 and 7. This is where I need to take the odd
  // byteswapping into account. The highest nibble is now in
  // the 7th byte, not the 6th.
  result[8] = (result[8] & 0x3f) | 0x80; // variant 1 (DCE)
  result[7] = (result[7] & 0x0f) | 0x40; // version 4 (random)
  return result;
}

/**
 * \brief Convert a ZGY stored guid to a user readable string.
 *
 * \details
 * The expected input is a 16-byte array that matches what is stored
 * in the ZGY file i.e. a little-endian struct GUID. This is
 * not standard. RFC 4122 states that big-endian is expected.
 * This method should only be used in code that reads ZGY files.
 *
 * If you copy the code to use it in a different context you need
 * to remove the explicit byteswapping. I still cannot guarantee
 * that the logic will then match the standard apis in Linux and
 * Windows, or even whether those two are even compatible by each
 * other. But I think all three will match.
 * See doc/uuid.md for more details.
 *
 * Note that I have not defined a specific type for the guid;
 * it is just a std::array. So I cannot overload operator<< to
 * automatically invoke GUID::format().
 *
 * The contents of byte 8 high nibble i.e. (result[8]>>4):
 * -    0..7    -- NCS, a very old format.
 * -    8,9,a,b -- variant 1, RFC 4122/DCE 1.1.
 * -    c,d     -- variant 2, old Microsoft.
 * -    e,f     -- reserved.
 *
 * Most systems today use variant 1.
 * Variant 1 should be converted from a struct UUID to big-endian if a
 * binary format is desired. Variant 2 should according to the wiki page
 * be byteswapped differently but that information appears to be incorrect
 * or for a different context. Note that byte 8 is in the not byteswapped
 * part. Byte 6, indicating the version, is in the byteswapped part and
 * will be found in byte 7 if the guid struct was not converted to big-endian.
 * The remaining variants are according to the RFC not expected to be
 * interoperable between systems. So the binary to sting conversion
 * should not matter. Treating them as variant 1 should be good enough.
 */
std::string
GUID::format(const GUID::guid_bytes_t& guid_in)
{
  GUID::guid_bytes_t guid = guid_in;
  //bool variant_2 = (guid_in[8] & 0xE0) == 0xC0;
  // Mimic the behavior in rpcrt4.dll on little-endian machines,
  // always convert as if the input was cast from a struct UUID
  // to a byte array.
  if (true) {
    // First part byteswaps as an uint32_t
    guid[0] = guid_in[3];
    guid[1] = guid_in[2];
    guid[2] = guid_in[1];
    guid[3] = guid_in[0];
    // Second and third part byteswaps as two uint16_t
    guid[4] = guid_in[5];
    guid[5] = guid_in[4];
    guid[6] = guid_in[7];
    guid[7] = guid_in[6];
    // Remaining two parts are not byteswapped.
  }
  // RFC 4122 requires lower case output, case-insensitive input.
  const char hex[17]{"0123456789abcdef"};
  char result[37]{0};
  char *cp = &result[0];
  char *end = &result[sizeof(result)-1];
  for (std::size_t ii = 0; ii < guid.size() && cp < end; ++ii) {
    *cp++ = hex[(guid[ii] >> 4) & 0x0f];
    *cp++ = hex[(guid[ii] >> 0) & 0x0f];
    if (ii == 3 || ii == 5 || ii == 7 || ii == 9)
      *cp++ = '-';
  }
  *cp++ = '\0';
  return std::string(result);
}

std::ostream&
Formatters::operator<<(std::ostream& os, const GUID& guid)
{
  os << guid.toString();
  return os;
}

} // namespace
