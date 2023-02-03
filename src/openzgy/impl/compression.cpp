// Copyright 2017-2020, Schlumberger
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

#include "compression.h"
#include "exception.h"

#include <algorithm>
#include <sstream>

/**
 * \file compression.cpp
 * \brief Compression factory.
 */

namespace InternalZGY {
#if 0
}
#endif

/**
 * \class CompressFactoryImpl
 * \brief Registry of known compress and decompress algorithms.
 *
 * The compression and decompression algorithms are completely
 * separate but we might as well handle both in the same class.
 *
 * Thread safety: Not thread safe by design.
 * Registering plug-ins is thread safe.
 */

/**
 * \brief Private method to get the compressor registry.
 *
 * The reason for this somewhat convoluted implementation, compared to
 * a "static [de]comp_reg_t" member declated in the class, is that the
 * instance gets initialized in a static constructor *and* members
 * gets added to in another static constructor. It is normally not
 * possible to control the order of those. Using this function forces
 * the instance to initialize before we use it.
 */
CompressFactoryImpl::comp_reg_t&
CompressFactoryImpl::_getCompressRegistry()
{
  static comp_reg_t instance;
  return instance;
}

/**
 * \copydoc _getCompressRegistry
 */
CompressFactoryImpl::decomp_reg_t&
CompressFactoryImpl::_getDecompressRegistry()
{
  static decomp_reg_t instance;
  return instance;
}

/**
 * \copydoc _getCompressRegistry
 */
std::mutex&
CompressFactoryImpl::_getMutex()
{
  static std::mutex mutex;
  return mutex;
}

/**
 * \brief Register a function that will be called to create a
 * compressor function.
 * Thread safety: The factory is synchronized using a lock.
 */
void
CompressFactoryImpl::registerCompressor(const std::string& name, const compfactory_t& fn)
{
  std::lock_guard<std::mutex> lk(_getMutex());
  if (fn)
    _getCompressRegistry()[name] = fn;
  else
    _getCompressRegistry().erase(name);
}

/**
 * \brief Register a function that can decompress one or more
 * types of compressed data.
 *
 * Pass an empty function in order to remove the registration.
 *
 * The registered functions will be called in reverse order of
 * registration until one of them indicates that it has decompressed
 * the data. The supplied name is only for information, and technically
 * doesn't even need to be unique. But it really should match the
 * name of the compressor.
 *
 * Thread safety: The factory is synchronized using a lock.
 */
void
CompressFactoryImpl::registerDecompressor(const std::string& name, const decompressor_t& fn)
{
  std::lock_guard<std::mutex> lk(_getMutex());
  if (fn) {
    std::pair<std::string, decompressor_t> e(name, fn);
    _getDecompressRegistry().insert(_getDecompressRegistry().begin(), e);
  }
  else {
    decomp_reg_t& reg = _getDecompressRegistry();
    decomp_reg_t old(reg);
    reg.clear();
    for (const auto& e : old)
      if (e.first != name)
        reg.push_back(e);
  }
}

/**
 * Return the names of all compressors known to the system.
 * This is primarily for logging, but might in principle be used
 * in a GUI to present a list of compressors to choose from.
 * The problem with that is how to handle the argument list.
 *
 * Thread safety: The factory is synchronized using a lock.
 */
std::vector<std::string>
CompressFactoryImpl::knownCompressors()
{
  std::lock_guard<std::mutex> lk(_getMutex());
  std::vector<std::string> result;
  for (const auto& it: _getCompressRegistry())
    result.push_back(it.first);
  std::sort(result.begin(), result.end());
  return result;
}

/**
 * Return the names of all decompressors known to the system.
 * This is primarily for logging.
 *
 * Thread safety: The factory is synchronized using a lock.
 */
std::vector<std::string>
CompressFactoryImpl::knownDecompressors()
{
  std::lock_guard<std::mutex> lk(_getMutex());
  std::vector<std::string> result;
  for (const auto& it: _getDecompressRegistry())
    result.push_back(it.first);
  std::sort(result.begin(), result.end());
  return result;
}

/**
 * Given a compressor name (passed in when creating a file)
 * return the actual compressor if it has been registered.
 *
 * Thread safety: The factory is synchronized using a lock.
 */
compressor_t
CompressFactoryImpl::getCompressor(const std::string& name, const std::vector<std::string>& args)
{
  {
    std::lock_guard<std::mutex> lk(_getMutex());
    if (_getCompressRegistry().count(name))
      return _getCompressRegistry()[name](args);
  } // unlock, because knownCompressors() re-locks.
  std::stringstream ss;
  ss << "Compression algorithm \"" << name << "\""
     << " not recognized. Must be one of (";
  int n{0};
  for (const std::string& name : knownCompressors())
    ss << (n++ ? "," : "") << name;
  ss << ")";
  throw OpenZGY::Errors::ZgyMissingFeature(ss.str());
}

/**
 * Loop over all registered decompressors and try to find one that
 * can handle this particular brick. Raises an error if none found.
 * See CompressPlugin.decompress() for parameter descriptions.
 *
 * Thread safety: The factory is synchronized using a lock.
 * The lock is dropped before the actual decompression is done.
 */
rawdata_t
CompressFactoryImpl::decompress(const rawdata_t& cdata, BrickStatus status, const index3_t& shape)
{
  rawdata_t result {nullptr, 0};
  if (status != BrickStatus::Compressed)
    throw OpenZGY::Errors::ZgyInternalError("Tried to decompress uncompressed data.");
  // TODO-Performance, need to copy the registry because the factory's lock
  // must not be held while decompressing. Maybe the decompress mechanism
  // should be changed to store the algorithm number in each brick so we
  // don't need to loop at this point. As a side benefit it would be simple
  // to list all compression algorithms used in a particular file.
  decomp_reg_t registry;
  {
    std::lock_guard<std::mutex> lk(_getMutex());
    registry = _getDecompressRegistry();
  }
  for (const auto& it : registry) {
    result = it.second(cdata, status, shape);
    if (result.first && result.second)
      return result;
  }
  throw OpenZGY::Errors::ZgyFormatError("Compression algorithm not recognized.");
}

} // namespace
