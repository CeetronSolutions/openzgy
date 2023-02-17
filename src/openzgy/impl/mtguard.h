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

#pragma once

#include "../declspec.h"

#include <exception>
#include <atomic>
#include <functional>
#include <cstdint>
#include <string>
#include <vector>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Helper class to avoid exceptions escaping from inside an OpenMP loop.
 * If something goes wrong the loop will still run to completion
 * but each iteration after something went wrong will be very fast.
 *
 * If the OpenMP loop has an ordered section then you probably need
 * protect that section explicitly. Using a second call to run()
 * on the same instance.
 *
 * Example usage:
 *
 * \code
 *     void example() {
 *       MTGuard guard;
 *     #pragma omp parallel for
 *       for (std::int64_t ii = 0; ii < loopcount; ++ii) {
 *         guard.run([&](){
 *           // ACTUAL CODE GOES HERE
 *         });
 *       }
 *       guard.finished();
 *     }
 * \endcode
 */
class OPENZGY_TEST_API MTGuard
{
private:
  std::atomic<int> _errors;
  std::exception_ptr _first_ex;
  std::string        _debug_name;
  int                _debug_requested;
  bool               _debug_reported;
  std::vector<int>   _debug_threads_used;
  // cannot be inside class on Windows.
  //static thread_local int _debug_nesting;

public:
  explicit MTGuard();
  MTGuard(const std::string& name, int requested);
  ~MTGuard();
  MTGuard(const MTGuard&) = delete;
  MTGuard& operator=(const MTGuard&) = delete;
  bool failed() const;
  void run(const std::function<void()>& fn);
  void finished();

protected:
  void fail();
};

/**
 * Helper class that extends MTGuard to also handle progress reports
 * and cancellation from inside an OpenMP parallel region.
 *
 * Example usage:
 *
 * \code
 *     void examplewithprogress() {
 *       ProgressWithDots p1;
 *       MTGuardWithProgress guard(std::ref(p1), total);
 *     #pragma omp parallel for
 *       for (std::int64_t ii = 0; ii < loopcount; ++ii) {
 *         guard.run([&](){
 *           // ACTUAL CODE GOES HERE
 *           guard.progress();
 *         });
 *       }
 *       guard.finished();
 *     }
 *\endcode
 */
class OPENZGY_TEST_API MTGuardWithProgress: public MTGuard
{
public:
  typedef std::function<bool(std::int64_t, std::int64_t)> progress_fn;
private:
  progress_fn _progress;
  std::int64_t _total;
  std::atomic<std::int64_t> _done;
  int _last_done_by;

public:
  MTGuardWithProgress(const progress_fn progress, std::int64_t total);
  MTGuardWithProgress(const progress_fn progress, std::int64_t total, const std::string& name, int requested);
  MTGuardWithProgress(const MTGuardWithProgress&) = delete;
  MTGuardWithProgress& operator=(const MTGuardWithProgress&) = delete;
  void progress(std::int64_t steps = 1);
  void finished();
};

} // namespace
