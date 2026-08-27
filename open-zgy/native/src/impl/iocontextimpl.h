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

/**
 * \file iocontext.h
 * \brief Backend specific context.
 *
 * Class IOContext and derivatives are used to hold backend specific
 * information such as authorization tokens.
 */

#pragma once

#include <cstdint>
#include <string>
#include <map>
#include <functional>
#include <memory>

namespace OpenZGY {
  class SeismicStoreIOContextV2;
  class SeismicStoreIOContext;
}

namespace InternalZGY {
#if 0
}
#endif

class SeismicStoreIOContextV2Impl
{
  friend OpenZGY::SeismicStoreIOContextV2;

public:
  typedef std::map<std::string,std::string> http_headers_t;

private:
  const OpenZGY::SeismicStoreIOContext* parent_;
  http_headers_t  http_headers_;

public:
  explicit SeismicStoreIOContextV2Impl(const OpenZGY::SeismicStoreIOContext*);

public:
  const http_headers_t& ctxt_httpHeaders() const {return http_headers_; }

  // At some point, all data members should be moved to the V2 pimpl.
  // In the InternalZGY namespace, only the impl type will be known.
  // This will be a breaking change in the API. The accessors can be
  // iplemented already, but this will give a messy circular dependency
  // between iocontext.h and iocontextimpl.h. Or the need to move the
  // body of each of the methods out to a .cpp file.
  //const std::string& ctxt_sdurl()    const {return parent_->_sdurl; }
  //const std::string& ctxt_sdapikey() const {return parent_->_sdurl; }
  //const std::string& ctxt_sdtoken()  const {return parent_->_sdurl; }
};

} // namespace
