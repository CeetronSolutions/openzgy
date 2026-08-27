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

#include "../declspec.h"
#include <string>

/**
 * \file environment.h
 * \brief Reading environment variables.
 */
namespace InternalZGY {

/**
 * Tiny helper class to access the environment variables of the process.
 * Basically this is just a wrapper around getenv, but there are
 * some Linux/Windows differences and it is a bad idea to duplicate
 * that special handling too many places.
 *
 * Thread safety:
 * Modification by calls to ::putenv() may lead to a data race.
 * Application code is not expected to call that function at all,
 * and in that case concurrent calls are safe.
 *
 * This class used to provide putNumericEnv and putStringEnv
 * and there was also a helper class to add push/pop semantics.
 * Those are removed (yagni). Restore from git if you need them.
 */
class OPENZGY_TEST_API Environment
{
public:
  static int         getNumericEnv(const char *name, int dflt = 0);
  static std::string getStringEnv(const char *name, const char *dflt = 0);
private: // class has static members only, should not be instantiated.
  Environment() = delete;
  Environment& operator=(const Environment&) = delete;
};

} // namespace
