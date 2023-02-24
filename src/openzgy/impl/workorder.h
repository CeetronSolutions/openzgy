// Copyright 2017-2022, Schlumberger
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

#include <functional>
#include <cstdint>
#include <memory>
#include <string>

#include "../declspec.h"

namespace InternalZGY {
#if 0
}
#endif

class SummaryTimer;
class WorkOrder;

class WorkOrderRunner
{
public:
  typedef std::function<void(std::int64_t)> workstep_fn;
  enum class Provider { SERIAL=0, OMP, TP };
public:
  /**
   * Helper for executing a loop, possibly multi threaded.
   * Exceptions may be collected and the first exception,
   * if any, be re-thrown after the loop is done.
   */
  static OPENZGY_TEST_API void parallelFor(
       Provider provider,
       std::int64_t count,
       std::int64_t threadcount,
       const workstep_fn& workstep);

  static OPENZGY_TEST_API void parallelFor(
       Provider provider,
       std::int64_t count,
       const workstep_fn& workstep);

  static OPENZGY_TEST_API void parallelFor(
       std::int64_t count,
       std::int64_t threadcount,
       const workstep_fn& workstep);

  static OPENZGY_TEST_API void parallelFor(
       std::int64_t count,
       const workstep_fn& workstep);

  static OPENZGY_TEST_API Provider defaultProvider();
  static OPENZGY_TEST_API std::pair<int, int> debugPoolSize();
  static OPENZGY_TEST_API std::pair<int, int> debugPoolSize(Provider provider);
  static OPENZGY_TEST_API int debugThreadsInUse();
};

} // namespace
