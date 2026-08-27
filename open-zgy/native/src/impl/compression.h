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

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <map>
#include <mutex>

#include "../declspec.h"
#include "enum.h"

/**
 * \file compression.h
 * \brief Factory for handling compression.
 */

namespace InternalZGY {
#if 0
}
#endif

// See the .cpp file for documentation.
class OPENZGY_TEST_API CompressFactoryImpl
{
public:
  //typedef std::function<rawdata_t(const rawdata_t&,const index3_t&)> compressor_t;
  typedef std::function<compressor_t(const std::vector<std::string>&)> compfactory_t;
  typedef std::function<rawdata_t(const rawdata_t&,BrickStatus,const index3_t&)> decompressor_t;

public:
  static void registerCompressor(const std::string& name,
                                 const compfactory_t& fn);
  static void registerDecompressor(const std::string& name,
                                   const decompressor_t& fn);
  static std::vector<std::string> knownCompressors();
  static std::vector<std::string> knownDecompressors();

  static compressor_t getCompressor(const std::string& name,
                                    const std::vector<std::string>& args);
  static rawdata_t decompress(const rawdata_t& cdata,
                              BrickStatus status,
                              const index3_t& shape);
private:
  typedef std::map<std::string, compfactory_t> comp_reg_t;
  typedef std::vector<std::pair<std::string, decompressor_t>> decomp_reg_t;
  static comp_reg_t& _getCompressRegistry();
  static decomp_reg_t& _getDecompressRegistry();
  static std::mutex& _getMutex();
};

} // namespace


