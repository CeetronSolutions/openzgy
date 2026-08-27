// Copyright 2017-2023, Schlumberger
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

#include "iocontextimpl.h"
#include "environment.h"
#include <string>
#include <cstdint>
#include <iostream>

namespace InternalZGY {
#if 0
}
#endif

SeismicStoreIOContextV2Impl::SeismicStoreIOContextV2Impl(const OpenZGY::SeismicStoreIOContext* parent)
  : parent_(parent)
{
  using InternalZGY::Environment;

  // Allow passing headers in an environment variable. Use only if you
  // are desperate. Keys and values are separated by a colon. Use no
  // spaces before and after the colon. Multiple headers are separated
  // by backquotes. If a value is supposed to include a backquote,
  // this won't work. Just remember that the envirinment variable is a
  // kludge anyway. And don't complain.

  std::string hh = Environment::getStringEnv("OPENZGY_HTTP_HEADERS") + "`";
  std::size_t beg = 0;
  std::size_t end = hh.find_first_of('`');
  while (end != std::string::npos) {
    std::string one = hh.substr(beg, end - beg);
    std::size_t colon = one.find(':');
    if (colon != std::string::npos && colon != 0)
      http_headers_.insert
        (std::make_pair(one.substr(0, colon),
                        one.substr(colon+1, std::string::npos)));
    beg = end + 1;
    end = hh.find_first_of('`', beg);
  }
}

} // namespace
