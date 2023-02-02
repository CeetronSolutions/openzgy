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

#include "environment.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

namespace InternalZGY {

/**
 * \brief Get value as a string.
 *
 * Get an environment variable as a std::string.
 * If the environment variable is empty or missing then the provided
 * default value is returned. Or the empty string if there is no default.
 * Attempting to distinguish between empty and missing variable is not
 * portable. So this function won't let you.
 * An environment variable with an empty name is always consdered unset.
 */
std::string Environment::getStringEnv(const char *name, const char *dflt)
{
  std::string result;
  if (!name || !*name)
    return dflt ? std::string(dflt) : std::string();
#ifdef _WIN32
   size_t requiredSize = 0;
   getenv_s(&requiredSize, NULL, 0, name);
   if (requiredSize > 0) {
     char *value = (char*)malloc(requiredSize * sizeof(char));
     if (value) {
       if (getenv_s(&requiredSize, value, requiredSize, name) == 0) {
         value[requiredSize-1] = '\0'; // ensure null terminated
         if (value[0] != '\0') {
           result = value;
         }
       }
       free(value);
     }
   }
#else
  const char *value = getenv(name);
  if (value != NULL && value[0] != '\0') {
    result = value;
  }
#endif
  return result.empty() ? dflt ? std::string(dflt) : std::string() : result;
}

/**
 * \brief Get value as a number.
 *
 * Get an environment variable as a number.
 * If the environment variable is missing or empty
 * then the supplied default value is returned instead.
 * If no default value was given then 0 is returned.
 * Trailing garbage is ignored, so if the variable
 * doesn't start with a digit the function silently
 * returns 0. I.e. NOT the provided default value.
 */
int Environment::getNumericEnv(const char *name, int dflt)
{
  // Note: If you are offended by the redundant creation of a std::string()
  // then it is possible to inline all of getStringEnv() here. But this should
  // not be a time critical function.
  std::string value = getStringEnv(name);
  if (value.empty())
    return dflt;
  else
    return atoi(value.c_str());
}

} // namespace
