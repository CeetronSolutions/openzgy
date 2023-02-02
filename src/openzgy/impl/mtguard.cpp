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

#include "mtguard.h"

#include <iostream>
#include <sstream>
#include <omp.h>

namespace InternalZGY {
#if 0
}
#endif

// Logically a member of class MTGuard, but Windows doesn't like that.
static thread_local std::atomic<int> MTGuard__debug_nesting(0);

MTGuard::MTGuard()
  : _errors(0)
  , _first_ex()
  , _debug_name("anonymous")
  , _debug_requested(0)
  , _debug_reported(true) // set false to enable reporting
{
}

/**
 * \brief Allow caller to provide more information for debugging and logging.
 * \details
 * It is ok to let the code permanently use this overloaded constructor,
 * although it does add more clutter to the code.
 */
MTGuard::MTGuard(const std::string& name, int requested)
  : _errors(0)
  , _first_ex()
  , _debug_name(name)
  , _debug_requested(requested)
  , _debug_reported(true) // set false to enable reporting
{
}

/**
 * When the Guard instance is destructed it is too late to throw any
 * pending exception to the application. Because *no* destructor
 * should *ever* throw an exception.
 */
MTGuard::~MTGuard()
{
  if (_first_ex) {
    try {
      std::rethrow_exception(_first_ex);
    }
    catch (const std::exception& ex) {
      std::cerr << "EXCEPTION inside OpenMP loop: " << ex.what() << std::endl;
    }
    catch (...) {
      std::cerr << "EXCEPTION inside OpenMP loop." << std::endl;
    }
  }
}

/**
 * The user of this class might do: if (guard.failed()) continue
 * at the start of the loop. If the loop only consists of calls
 * to run() then this is not necessary.
 */
bool
MTGuard::failed() const {
  return _errors.load() != 0;
}

/**
 * Do something, as specified in the input function, if everything is ok.
 * Catch any exceptions and if caught mark this guared instance as failed.
 * Subsequent calls to run() will then not do anything at all.
 * No exceptions will propagate out of run(). When calling finished()
 * the first exception, if any, is then thrown.
 */
void
MTGuard::run(const std::function<void()>& fn) {
  // Used to debug which OpenMT loops are run serially.
  // Only works if the loop uses MTGuard. And even in that
  // case your mileage may vary.
  if (!_debug_reported && omp_get_thread_num() == 0) {
    std::stringstream ss;
    ss << "OpenMP \"" << _debug_name << "\""
       << " level " << MTGuard__debug_nesting + 1
       << " wants " << _debug_requested
       << " got " << omp_get_num_threads()
       << " threads.\n";
    std::cerr << ss.str() << std::flush;
    _debug_reported = true;
  }
  if (_errors.load() == 0) {
    try {
      MTGuard__debug_nesting++;
      fn();
      MTGuard__debug_nesting--;
    }
    catch (...) {
      MTGuard__debug_nesting--;
      fail();
    }
  }
}

/**
 * Call ths outside the parallel region so that any pending
 * exception can be re-thrown. No-op if there was no exception.
 * Remember to call this before the instance goes out of scope.
 * Otherwise any exception just gets printed to std::cerr and
 * swallowed.
 */
void
MTGuard::finished() {
  if (_errors.load() != 0 && _first_ex) {
    std::exception_ptr throwme = _first_ex;
    _first_ex = std::exception_ptr();
    std::rethrow_exception(throwme);
  }
}

/**
 * To be called from inside a catch statement. Stores the exception
 * so it can be rethrown outside the parallel region.
 */
void
MTGuard::fail() {
  if (_errors.fetch_add(1) == 0) {
    _first_ex = std::current_exception();
    if (!_first_ex)
      _first_ex = std::make_exception_ptr(std::runtime_error("fail() with no current exception."));
  }
}

MTGuardWithProgress::MTGuardWithProgress(const progress_fn fn, std::int64_t total)
  : MTGuard()
  , _progress(fn)
  , _total(total)
  , _done(0)
  , _last_done_by(-1)
{
}

/**
 * \brief Allow caller to provide more information for debugging and logging.
 * \details
 * It is ok to let the code permanently use this overloaded constructor,
 * although it does add more clutter to the code.
 */
MTGuardWithProgress::MTGuardWithProgress(const progress_fn fn, std::int64_t total, const std::string& name, int requested)
  : MTGuard(name, requested)
  , _progress(fn)
  , _total(total)
  , _done(0)
  , _last_done_by(-1)
{
}

/**
 * Keep track of the progress i.e. how far we have come in the
 * parallel for. User callbacks are invoked but only from thread
 * zero. So the caller will see a monotonously increasing "done"
 * but will not be called for each step. If an error occurs then
 * no more callbacks are invoked, except possibly the last one
 * with done==total.
 *
 * The progress callback can return false if it wants to abort.
 * In this scenario (inside an OpemMP loop) we treat an exception
 * from inside the callback the same way.
 */
void
MTGuardWithProgress::progress(std::int64_t steps)
{
  if (!failed() && _progress && _total != 0 && steps > 0) {
    // All threads keep track of the total progress.
    const std::int64_t localdone = _done.fetch_add(steps)+steps;
    if (omp_get_thread_num() == 0) {
      try {
        if (!_progress(localdone, _total))
          throw std::runtime_error("aborted");
      }
      catch (const std::exception&) {
        fail();
      }
    }
    if (localdone == _total)
      _last_done_by = omp_get_thread_num();
  }
}

void
MTGuardWithProgress::finished() {
  // Test _last_done_by > 0 if you don't want to send the final
  // all-done if the process has been aborted or caught an exception.
  if (_progress && _total != 0 && _last_done_by != 0)
    if (!_progress(_total, _total))
      throw std::runtime_error("aborted");
  MTGuard::finished();
}

} // namespace
