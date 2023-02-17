// Copyright 2017-2022, SLB
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
#include <utility>

namespace InternalZGY { namespace WindowsTools {
#if 0
}}
#endif

/**
 * Debug function to return the number of threads created and
 * destroyed since the last call. Used to debug problems with OpenMP.
 * In particular, whether OpenMP is using a thread pool or not.
 * Caveat: Possibly too much overhead to allow leaving it in production code.
 */
extern std::pair<int, int> threadReportPair();

/**
 * As threadReportPair, but return a formatted string.
 */
extern std::string threadReportString();

/**
* Run a test that measures the time to create and destroy a thread.
* Apparently the result depends on where the test is run. e.g. Petrel
* takes 10 times longer than ZgyTool. This is why the test code needs
* to be inside the library. For ad-hoc testing a call to this function
* can be added e.g. in GenLodImpl::call(). 
*/
void threadTimeTest(int count);

}} // namespace