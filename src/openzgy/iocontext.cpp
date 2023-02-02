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

#include "iocontext.h"
#include "api.h" // Only for IZgyUtils
#include "impl/environment.h"

#include <sstream>

namespace OpenZGY {
#if 0
}
#endif

namespace
{
  std::string _redact(const std::string& s)
  {
    return s.size() < 20 ? s : s.substr(0,4) + "..." + s.substr(s.size()-4);
  }
}

/**
 * The class needs at least one virtual member to make dynamic_cast work.
 */
IOContext::~IOContext()
{
}

/** \cond SSTORE */

SeismicStoreIOContext::SeismicStoreIOContext()
{
  using InternalZGY::Environment;

  // Chicken and egg problem, segsplit() and writethreads() currently
  // have some special interaction where each of them looks at the
  // other value.
  _segsplit = 1;
  _writethreads = 1;

  sdurl(Environment::getStringEnv("OPENZGY_SDURL"));
  sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"));
  sdtoken(Environment::getStringEnv("OPENZGY_TOKEN"), "");
  maxsize(Environment::getNumericEnv("OPENZGY_MAXSIZE_MB", 64));
  maxhole(Environment::getNumericEnv("OPENZGY_MAXHOLE_MB", 2));
  aligned(Environment::getNumericEnv("OPENZGY_ALIGNED_MB", 0));
  segsize(Environment::getNumericEnv("OPENZGY_SEGSIZE_MB", 256));
  segsplit(Environment::getNumericEnv("OPENZGY_SEGSPLIT", 8));
  iothreads(Environment::getNumericEnv("OPENZGY_NUMTHREADS_IO", 1));
  cputhreads(Environment::getNumericEnv("OPENZGY_NUMTHREADS_CPU", 1));
  writethreads(Environment::getNumericEnv("OPENZGY_NUMTHREADS_WRITE", 1));
  legaltag(Environment::getStringEnv("OPENZGY_LEGALTAG"));
  writeid(Environment::getStringEnv("OPENZGY_WRITEID"));
  seismicmeta(Environment::getStringEnv("OPENZGY_SEISMICMETA"));
  setRoAfterWrite(Environment::getNumericEnv("OPENZGY_RO_AFTER_WRITE", 1) > 0);
  forceRoBeforeRead(Environment::getNumericEnv("OPENZGY_RO_BEFORE_READ", 0) > 0);
  forceRwBeforeWrite(Environment::getNumericEnv("OPENZGY_RW_BEFORE_WRITE", 0) > 0);
  retryCount(Environment::getNumericEnv("OPENZGY_SD_BACKOFF", -1));

  if (this->_writethreads != 1)
    segsplit(1);
}

SeismicStoreIOContext&
SeismicStoreIOContext::credentialsFrom(std::shared_ptr<IZgyUtils> utils)
{
  this->sdtokencb([utils](){return utils->idtoken();});
  return *this;
}

std::string
SeismicStoreIOContext::toString() const
{
  std::stringstream ss;
  ss << "SD context:\n"
     << "  sdurl:    \"" << _sdurl << "\"\n"
     << "  sdapikey: \"" << _redact(_sdapikey) << "\"\n"
     << "  sdtoken:  \"" << _redact(_sdtoken) << "\"\n"
     << "  maxsize:  " << _maxsize / (1024*1024) << " MB\n"
     << "  maxhole:  " << _maxhole / (1024*1024) << " MB\n"
     << "  aligned:  " << _aligned / (1024*1024) << " MB\n"
     << "  segsize:  " << _segsize / (1024*1024) << " MB\n"
     << "  segsplit: " << _segsplit << "\n"
     << "  threads:  " << _iothreads << " I/O, "
     << _cputhreads   << " CPU, "
     << _writethreads << " Write\n"
     << "  legaltag: \"" << _legaltag << "\"\n"
     << "  writeid:  \"" << _writeid << "\"\n"
     << "  seismeta: \"" << _seismicmeta << "\"\n"
     << "  romode:   "
     << (_set_ro_after_write? "ro_after_write":" ")
     << (_force_ro_before_read? "ro_before_read":" ")
     << (_force_rw_before_write? "rw_before_write":" ")
     << "\n";
  return ss.str();
}

std::shared_ptr<IOContext>
SeismicStoreIOContext::clone() const
{
  return std::make_shared<SeismicStoreIOContext>(*this);
}

/** \endcond */

} // namespace
