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

#include "test_all.h"
#include "test_utils.h"
#include "../impl/workorder.h"
#include "../impl/timer.h"
#include "../impl/fancy_timers.h"

#include <thread>
#include <chrono>
#include <memory>
#include <atomic>
#include <map>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <omp.h>

#ifdef _WIN32
#define WINDOWS_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#endif

using namespace InternalZGY;

namespace InternalZGY {
  class WorkOrder;
}

namespace Test {
  class TestWorkOrder {
  public:
    static void test_vanilla(InternalZGY::WorkOrderRunner::Provider p);
    static void test_poolsize(InternalZGY::WorkOrderRunner::Provider p);
    static void test_chunksize(InternalZGY::WorkOrderRunner::Provider p);
    static void test_throw(InternalZGY::WorkOrderRunner::Provider p);
    static void test_single(InternalZGY::WorkOrderRunner::Provider p);
    static void test_last(InternalZGY::WorkOrderRunner::Provider p);
  };
}

namespace {
  // gcc 9 and above have gettid() in the system headers, but
  // devtoolset-10 doesn't. Use a different name to avoid conflicts.
#if defined(_WIN32)
  static int my_gettid() { return ::GetCurrentThreadId(); }
#elif defined(__GNUC__) && defined(SYS_gettid)
  static int my_gettid() { return syscall(SYS_gettid); }
#else
  // If the toolset has both SYS_gettid and gettid(), the former
  // will be used. No big deal, I think.
  static int my_gettid() { return gettid(); }
#endif
}

namespace {
  std::string format(InternalZGY::WorkOrderRunner::Provider p)
  {
    typedef InternalZGY::WorkOrderRunner::Provider Provider;
    switch (p)
    {
    case Provider::SERIAL: return "SERIAL";
    case Provider::OMP:    return "OMP";
    case Provider::TP:     return "TP";
    default:               return "Unknown";
    }
  }
}

/**
 * Check that all iterations in the loop are executed.
 * Check that multi-threading gives at least 3x speedup.
 */
void
Test::TestWorkOrder::test_vanilla(InternalZGY::WorkOrderRunner::Provider p)
{
  InternalZGY::Timer outertime(true, "test_vanilla");
  std::atomic<int> loops(0);
  std::atomic<int> innertime(0);
  const bool serial = (p == InternalZGY::WorkOrderRunner::Provider::SERIAL);
  const int iterations = (serial ? 42 : 1007);
  WorkOrderRunner::parallelFor(p, iterations, [&](std::int64_t)
    {
      int sleeptime = Test_Utils::random_delay(16);
      loops.fetch_add(1);
      innertime.fetch_add(sleeptime);
    });
  outertime.stop();

  const double itime = innertime; // milliseconds
  const double otime = outertime.getTotal() * 1000; // milliseconds
  const double speedup = (otime > 0 ? itime / otime : 0);
  if (verbose() > 0)
    std::cout << "\nProvider " << format(p)
              << " Times: "
              << "inner " << itime << " ms, "
              << "outer " << otime << " ms. "
              << "speedup " << speedup << " "
              << "loops " << loops.load() << ".\n"
              << "Pool:"
              << " size " << WorkOrderRunner::debugPoolSize(p).first
              << " used " << WorkOrderRunner::debugPoolSize(p).second
              << std::endl;

  // Must run the exact number of requested iterations.
  TEST_EQUAL(loops.load(), iterations);

  // Inner time is accumulated requested sleep time, which is
  // guaranteed to be within 0.5x and 2.0x the requested time.
  // Had it been measured instead then a slop value might be needed.
  TEST_CHECK(itime >= iterations * 8 && itime <= iterations * 32);

  // The Windows build servers we are currently using
  // have just 2 cores, so there is no way we can get
  // 3x speedup by multi-threading.
  if (p == InternalZGY::WorkOrderRunner::Provider::OMP && omp_get_max_threads() < 6)
    return;

  // If the MT loop doesn't speed up by at least 3x then this is a fail.
  if (serial)
    TEST_CHECK(speedup >= 0.5 && speedup < 1.5);
  else
    TEST_CHECK(speedup >= 3.0);

  // No current test machine has > 100 cores. Check for suspiciously
  // short outer time. If this test fails, it might be a false positive.
  TEST_CHECK(speedup <= 100.0);
}

/**
 * Check how many different threads were used for the task.
 */
void
Test::TestWorkOrder::test_poolsize(InternalZGY::WorkOrderRunner::Provider p)
{
  const bool serial = (p == InternalZGY::WorkOrderRunner::Provider::SERIAL);
  const int iterations = (serial ? 42 : 1007);
  std::map<int, int> threads;
  std::mutex mutex;
  WorkOrderRunner::parallelFor(p, iterations, [&](std::int64_t)
    {
      int tid = my_gettid();
      Test_Utils::random_delay(8);
      std::lock_guard<std::mutex> lk(mutex);
      threads.insert(std::pair<int,int>(tid, 0));
      ++threads[tid];
    });

  if (verbose() > 0) {
    const int me = my_gettid();
    std::cout << "\n   Threads used: " << threads.size()
              << ", team leader " << me << "\n"
              << "   count thread_id\n";
    for (const auto& it : threads)
      std::cout << "   "
                << std::setw(5) << it.second
                << "  " << it.first
                << (it.first == me ? "*" : " ")
                << "\n";
  }

  // Assumptions on how many threads were involved.
  // For very powerful test machines or very slow and heavily
  // loaded machines there might be spurious errors.

  // The Windows build servers we are currently using
  // have just 2 cores, so there is no way we can get
  // 3x speedup by multi-threading.
  if (p == InternalZGY::WorkOrderRunner::Provider::OMP && omp_get_max_threads() < 6)
    return;

  if (serial) {
    TEST_EQUAL(threads.size(), 1);
  }
  else {
    TEST_CHECK(threads.size() >= 3);
    TEST_CHECK(threads.size() <= 101);
  }
}

/**
 * Check how many different threads were used for the task.
 */
void
Test::TestWorkOrder::test_chunksize(InternalZGY::WorkOrderRunner::Provider p)
{
  const bool serial = (p == InternalZGY::WorkOrderRunner::Provider::SERIAL);
  const int iterations = 42;
  const int numthreads = 5; // i.e. chunksize 11
  std::map<int, int> threads;
  std::mutex mutex;
  std::atomic<int> loops(0);
  WorkOrderRunner::parallelFor(p, iterations, numthreads, [&](std::int64_t)
    {
      int tid = my_gettid();
      Test_Utils::random_delay(8);
      std::lock_guard<std::mutex> lk(mutex);
      threads.insert(std::pair<int,int>(tid, 0));
      ++threads[tid];
      loops.fetch_add(1);
    });

  if (verbose() > 0) {
    const int me = my_gettid();
    std::cout << "\n   Threads used: " << threads.size()
              << ", team leader " << me << "\n"
              << "   count thread_id\n";
    for (const auto& it : threads)
      std::cout << "   "
                << std::setw(5) << it.second
                << "  " << it.first
                << (it.first == me ? "*" : " ")
                << "\n";
  }

  TEST_EQUAL(loops.load(), iterations);

  // Assumptions on how many threads were involved.
  // For very powerful test machines or very slow and heavily
  // loaded machines there might be spurious errors.
  // Note that the team leader won't necessarily find any work
  // to do. But then it won't show up in "threads" either.
  if (serial) {
    TEST_EQUAL(threads.size(), 1);
  }
  else {
    TEST_CHECK(threads.size() > 1);
    TEST_CHECK(threads.size() <= 5);
  }
}

/**
 * Exceptions inside a step are caught, and the first of them will
 * be re-thrown as if it came from the thread that started the loop.
 */
void
Test::TestWorkOrder::test_throw(InternalZGY::WorkOrderRunner::Provider p)
{
  const std::string msg{"Life, the universe, and everything"};
  try {
    WorkOrderRunner::parallelFor(p, 1007, [&](std::int64_t num)
      {
        if (num == 42)
          throw std::runtime_error(msg);
        Test_Utils::random_delay(8);
      });
    TEST_CHECK(false && "Did not get expected exception.");
  }
  catch (const std::runtime_error& ex) {
    TEST_EQUAL(std::string(ex.what()), msg);
  }
}

/**
 * Asking for zero or one iterations might trigger corner cases or shortcuts.
 */
void
Test::TestWorkOrder::test_single(InternalZGY::WorkOrderRunner::Provider p)
{
  std::atomic<int> ping{0};
  WorkOrderRunner::parallelFor(p, 1, [&](std::int64_t num)
    {
      ping.fetch_add(1);
    });
  TEST_EQUAL(ping.load(), 1);

  WorkOrderRunner::parallelFor(p, 0, [&](std::int64_t num)
    {
      ping.fetch_add(1);
    });
  TEST_EQUAL(ping.load(), 1);
}

/**
 * Occasionally but rarely the team leader gets to finish the
 * last thread. Explicitly check this corner case.
 * Only really relevant for the TP case, bot there is no
 * harm in running it for all.
 */
void
Test::TestWorkOrder::test_last(InternalZGY::WorkOrderRunner::Provider p)
{
  std::atomic<int> ping{0};
  int team_leader = my_gettid();
  WorkOrderRunner::parallelFor(p, 4, [&](std::int64_t num)
    {
      if (my_gettid() == team_leader)
        Test_Utils::random_delay(24);
      else
        Test_Utils::random_delay(2);
      ping.fetch_add(1);
    });
  TEST_EQUAL(ping.load(), 4);
}

// Check the logging (or yagni?)
// Check progress and abort handling (or yagni?)

namespace {
  class Register
  {
  public:
    typedef InternalZGY::WorkOrderRunner::Provider Provider;
    Register()
    {
      register_test("workorder.vanilla_ser",   [](){Test::TestWorkOrder::test_vanilla(Provider::SERIAL);});
      register_test("workorder.vanilla_omp",   [](){Test::TestWorkOrder::test_vanilla(Provider::OMP);});
      register_test("workorder.vanilla_tp",    [](){Test::TestWorkOrder::test_vanilla(Provider::TP);});
      register_test("workorder.poolsize_ser",  [](){Test::TestWorkOrder::test_poolsize(Provider::SERIAL);});
      register_test("workorder.poolsize_omp",  [](){Test::TestWorkOrder::test_poolsize(Provider::OMP);});
      register_test("workorder.poolsize_tp" ,  [](){Test::TestWorkOrder::test_poolsize(Provider::TP);});
      register_test("workorder.chunksize_ser", [](){Test::TestWorkOrder::test_chunksize(Provider::SERIAL);});
      register_test("workorder.chunksize_omp", [](){Test::TestWorkOrder::test_chunksize(Provider::OMP);});
      register_test("workorder.chunksize_tp" , [](){Test::TestWorkOrder::test_chunksize(Provider::TP);});
      register_test("workorder.throw_ser",     [](){Test::TestWorkOrder::test_throw(Provider::SERIAL);});
      register_test("workorder.throw_omp",     [](){Test::TestWorkOrder::test_throw(Provider::OMP);});
      register_test("workorder.throw_tp",      [](){Test::TestWorkOrder::test_throw(Provider::TP);});
      register_test("workorder.single_ser",    [](){Test::TestWorkOrder::test_single(Provider::SERIAL);});
      register_test("workorder.single_omp",    [](){Test::TestWorkOrder::test_single(Provider::OMP);});
      register_test("workorder.single_tp",     [](){Test::TestWorkOrder::test_single(Provider::TP);});
      register_test("workorder.last_ser",      [](){Test::TestWorkOrder::test_last(Provider::SERIAL);});
      register_test("workorder.last_omp",      [](){Test::TestWorkOrder::test_last(Provider::OMP);});
      register_test("workorder.last_tp",       [](){Test::TestWorkOrder::test_last(Provider::TP);});
    }
  } dummy;
} // namespace for registration
