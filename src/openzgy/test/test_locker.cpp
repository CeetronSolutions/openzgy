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

#include "test_all.h"
#include "../impl/locker.h"

#include "test_utils.h"
#include "../api.h"
#include "../iocontext.h"
#include "../exception.h"
#include "../impl/environment.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <memory>
#include <cstdint>
#include <cmath>
#include <numeric>
#include <limits>
#include <thread>
#include <chrono>
#include <algorithm>
#include <functional>

using namespace OpenZGY;
using namespace OpenZGY::Formatters;
using Test_Utils::must_throw;

namespace Test_LOCKER {
#if 0
}
#endif

class TestLocker
{
public:
  enum class State {Dead=0, Starting, Running, Working, Stopping};
  typedef std::function<bool(int, const std::string&)> LoggerFn;
  enum{MAXTHREADS=15};

private:
  std::mutex mutex_;
  std::condition_variable cv_;
  State state_[MAXTHREADS];
  std::shared_ptr<InternalZGY::Locker> locker_;
  bool shutting_down_;
  bool exception_in_worker_;

public:
    explicit TestLocker(std::shared_ptr<InternalZGY::Locker> locker)
    : mutex_()
    , cv_()
    , state_{}
    , locker_(locker)
    , shutting_down_(false)
    , exception_in_worker_(false)
  {
  }

  ~TestLocker()
  {
    shutdown();
  }

  static std::string format(State state)
  {
    switch (state)
    {
    case State::Dead: return "Dead";
    case State::Starting: return "Starting";
    case State::Running: return "Running";
    case State::Working: return "Working";
    case State::Stopping: return "Stopping";
    default: return "State_" + std::to_string(int(state));
    }
  }

  int countActive()
  {
    std::unique_lock<std::mutex> lck(mutex_);
    int active{0};
    for (int id=0; id<MAXTHREADS; ++id)
      if (state_[id] != State::Dead)
        ++active;
    return active;
  }

  void shutdown()
  {
    int active{0};
    std::cout << "Shutdown start" << std::endl;
    shutting_down_ = true;
    locker_->discard();
    {
      std::unique_lock<std::mutex> lck(mutex_);
      for (int id=0; id<MAXTHREADS; ++id) {
        if (state_[id] != State::Dead) {
          std::cout << "Shutdown: Worker " << id << " is " << format(state_[id])<< std::endl;
          // Dead: Leave it. Stopping: Still need to wait until dead.
          ++active;
          state_[id] = State::Stopping; // in case it wasn't already.
        }
      }
      if (active) {
        std::cout << "Shutdown tells "
                  << active << " threads to stop" << std::endl;
        cv_.notify_all();
      }
    }
    if (active)
      for (int id=0; id<MAXTHREADS; ++id)
        waitForState(id, State::Dead, -1);
    std::cout << "Shutdown ended" << std::endl;
  }

  void setState(int id, State state)
  {
    std::cerr << "Worker " << id << ": "
              << "Change state from " << format(state_[id])
              << " to " << format(state) << std::endl;
    state_[id] = state;
    cv_.notify_all();
  }

  /**
   * The worker, to be run in a separate thread until somebody tells it to stop.
   */
  void worker(int id, int blocknum, bool write)
  {
    {
      std::unique_lock<std::mutex> lck(mutex_);
      setState(id, State::Running);
    }

    // The thread will say in the Running state while waiting for the
    // fine grained lock to be allocated. Then it switches to Working,
    // keeping the lock, until an external change to the state which
    // will cause fn to return and the lock released.
    auto fn =
      [this,id]()
      {
        std::unique_lock<std::mutex> lck(this->mutex_);
        this->setState(id, State::Working);
        while (this->state_[id] == State::Working)
          this->cv_.wait(lck);
      };
    if (write)
      locker_->runWithWriteLock(blocknum, 1024, false, fn);
    else
      locker_->runWithReadLock(blocknum, fn);

    {
      std::unique_lock<std::mutex> lck(mutex_);
      setState(id, State::Dead);
    }
  }

  /**
   * The worker, to be run in a separate thread until somebody tells
   * it to stop.
   */
  static void worker_s(TestLocker* self, int id, int blocknum, bool write)
  {
    try {
      self->worker(id, blocknum, write);
    }
    catch (const OpenZGY::Errors::ZgyAborted&) {
      if (self->shutting_down_) {
        // This is normal when a shutdown is forced.
        std::cout << "Ok: ** Thread " + std::to_string(id) + " aborted **\n" << std::flush;
      }
      else {
        std::cout << "ERROR: ** Thread " + std::to_string(id) + " aborted **\n" << std::flush;
        self->exception_in_worker_ = true;
      }
    }
    catch (...) {
      std::cout << "ERROR: ** Thread " + std::to_string(id) + " got exception **\n" << std::flush;
      self->exception_in_worker_ = true;
    }
  }

  void waitForState(int id, State state, int timeout)
  {
    std::string name = "Worker " + std::to_string(id) + ": ";
    std::cout << name << "waitForState(" << format(state_[id]) << "->" << format(state) << ")" << std::endl;
    std::unique_lock<std::mutex> lck(mutex_);
    if (state_[id] == state)
      return;
    std::cout << name << "Waiting for " << format(state_[id]) << "->" << format(state) << std::endl;
    while (state_[id] != state) {
      if (timeout > 0) {
        if (cv_.wait_for(lck, std::chrono::seconds(5)) == std::cv_status::timeout)
          throw std::runtime_error(name + "Did not respond");
      }
      else {
        cv_.wait(lck);
      }
      std::cout << name << "Now in " << format(state_[id]) << " want " << format(state) << std::endl;
    }
    std::cout << name << "waitForState(" << format(state) << ") returns" << std::endl;
  }

  void startNoWait(int id, int blocknum, bool write)
  {
    std::string name = "Worker " + std::to_string(id) + ": ";
    std::unique_lock<std::mutex> lck(mutex_);
    if (id < 0 || id >= MAXTHREADS || state_[id] != State::Dead)
      throw new std::runtime_error(name + "Cannot start now.");
    setState(id, State::Starting);
    std::thread(&TestLocker::worker_s, this, id, blocknum, write).detach();
  }

  void start(int id, int blocknum, bool write, State target_state)
  {
    startNoWait(id, blocknum, write);
    waitForState(id, target_state, 5);
  }

  void stopNoWait(int id)
  {
    std::string name = "Worker " + std::to_string(id) + ": ";
    std::cout << name << "Being told to stop ?" << std::endl;
    std::lock_guard<std::mutex> lck(mutex_);
    if (id < 0 || id >= MAXTHREADS)
      throw new std::runtime_error(name + "Cannot stop this thread.");
    if (state_[id] != State::Dead) {
      std::cout << name << "Being told to stop" << std::endl;
      setState(id, State::Stopping);
    }
  }

  void stop(int id)
  {
    stopNoWait(id);
    waitForState(id, State::Dead, 5);
  }

  void complete(int id)
  {
    // The thread is supposed to be unblocked now, so wait for it to
    // be marked as Working (when the locker wakes the thread) and
    // then tell it to shut down and make sue that happend. Directly
    // to shutdown gas a race condition because there might be a
    // Stopping -> Working transition.
    waitForState(id, State::Working, 5);
    stop(id);
  }

  bool exceptionInWorker() const
  {
    return exception_in_worker_;
  }
};

static bool
myLogger(int pri, const std::string& msg)
{
  return false;
}

/**
 * Simple test case with three threads all reading from brick 1.
 * Note that a test failure might well cause a deadlock.
 */
static void
test_simple()
{
  auto locker = std::make_shared<InternalZGY::Locker>(4, myLogger);
  TestLocker test(locker);
  test.start(1, 1, false, TestLocker::State::Working);
  test.start(2, 1, false, TestLocker::State::Working);
  test.start(3, 1, false, TestLocker::State::Working);
  std::cout << locker->toString() << std::endl;
  std::cout << locker->toString(1) << std::endl;
  TEST_EQUAL(locker->toString(), "Lock info: Read: 3 held, 0 waiting. Write: 0 held, 0 waiting.");
  std::cout << "==== Stopping all ====" << std::endl;
  test.stopNoWait(2);
  test.stopNoWait(1);
  test.stopNoWait(3);
  std::cout << "==== Waiting for all to die ====" << std::endl;
  test.waitForState(3, TestLocker::State::Dead, 5);
  test.waitForState(2, TestLocker::State::Dead, 5);
  test.waitForState(1, TestLocker::State::Dead, 5);
  TEST_EQUAL(test.countActive(), 0);
  if (test.countActive()) {
    std::cout << "==== Shutdown ====" << std::endl;
    test.shutdown();
  }
  std::cout << "==== Leaving scope ====" << std::endl;
}

/**
 * A semi-random mixture of read and write.
 * Note that a test failure might well cause a deadlock.
 */
static void
test_mixed()
{
  auto locker = std::make_shared<InternalZGY::Locker>(4, myLogger);
  TestLocker test(locker);
  typedef TestLocker::State State;
  test.start(1, 101, false, State::Working); // Thread 1 -- read block 101
  test.start(2, 102, false, State::Working); // Thread 2 -- read block 102
  test.start(3, 103, false, State::Working); // Thread 3 -- read block 103
  test.start(4, 103, true,  State::Running); // Thread 4 -- write block 103 // TODO set state "Waiting" when we wait, not if no wait needed. Otherwise it might move on from Running to Working.
                                           // write blocked by read lock.
  test.start(5, 104, true,  State::Working); // Thread 5 -- write block 104
  test.start(6, 105, true,  State::Working); // Thread 6 -- write block 105
  test.start(7, 101, false, State::Working); // Thread 7 -- another reader 101
  test.start(8, 103, false, State::Working); // Thread 8 -- read block 103
                                           // 103 now 2 readers, 1 waiting writer
  test.start(9, 104, true,  State::Running); // Thread 9 -- write block 104
                                           // Already write locked, so wait.
  test.start(10,105, false, State::Running); // Thread 10 -- read block 105
                                           // Blocked by writer.
  test.start(11,105, false, State::Running); // Thread 11 -- read block 105
                                           // Another read is blocked,
  TEST_EQUAL(locker->toString(101), "Read lock (2)"); // thread 1 and 7 have read locks
  TEST_EQUAL(locker->toString(102), "Read lock"); // thread 2 has read lock.
  TEST_EQUAL(locker->toString(103), "Read lock (2)"); // thread 3 and 8 have read locks, 4 is waiting for write.
  TEST_EQUAL(locker->toString(104), "Write lock"); // thread 5 writing, 9 is waiting to write.
  TEST_EQUAL(locker->toString(105), "Write lock"); // thread 6 is writing, 10 and 11 waiting to read.
  TEST_EQUAL(locker->toString(), "Lock info: Read: 5 held, 2 waiting. Write: 2 held, 2 waiting.");

  // It is an error for threads 4, 9, 10, or 11 to have slipped past
  // "Runnning" and entered "Working" at this point. That is tricky to
  // test explicitly, because the threads don't have a special case
  // "about to block on a fine grained lock. Wait a few seconds and
  // check. Don't need to check each thred; the summary should
  // suffice.
  std::this_thread::sleep_for(std::chrono::seconds(2));
  TEST_EQUAL(locker->toString(), "Lock info: Read: 5 held, 2 waiting. Write: 2 held, 2 waiting.");

  test.complete(1); // Now only 1 read lock left on 101
  test.start(12, 101, true, State::Running); // still 1 read lock on 101
  test.complete(7); // Mo more read locks on 101, so 12 wakes up.
  test.complete(12); // 101 is now released.
  TEST_EQUAL(locker->toString(101), "Idle");

  test.complete(2); // 102, which had no contention, is now released.
  TEST_EQUAL(locker->toString(102), "Idle");

  test.complete(8); // still one read lock on 103
  test.complete(3); // last read lock on 103. Thread 4 wakes up and writes 103.
  test.complete(4); // 103 written and ought to be free now.
  TEST_EQUAL(locker->toString(103), "Idle");

  test.complete(5); // Was writing 104, 9 can now start writing it.
  test.complete(9); // 104 should now be released
  TEST_EQUAL(locker->toString(104), "Idle");

  test.complete(6); // was writing 105. Readers 10 and 11 may wake.
  test.complete(10);
  test.complete(11); // Block 105 is now released.
  TEST_EQUAL(locker->toString(105), "Idle");
  TEST_EQUAL(locker->toString(), "Lock info: Idle");

  test.shutdown();
  // Delayed reporting of error.
  TEST_CHECK(!test.exceptionInWorker());
}

/**
 * Global read- and write locks.
 */
static void
test_global()
{
  auto locker = std::make_shared<InternalZGY::Locker>(4, myLogger);
  TestLocker test(locker);
  typedef TestLocker::State State;
  test.start(1, 101, false, State::Working); // Thread 1 -- read block 101
  test.start(2, -1,  false, State::Working); // Thread 2 -- global read lock
  test.start(3, 101, true,  State::Running); // Thread 3 -- wait to write 101
  TEST_EQUAL(locker->toString(), "Lock info: Read: 2 held, 0 waiting. Write: 0 held, 1 waiting.");
  test.complete(1);
  TEST_EQUAL(locker->toString(), "Lock info: Read: 1 held, 0 waiting. Write: 0 held, 1 waiting.");
  test.complete(2); // 3 wakes and starts writing
  TEST_EQUAL(locker->toString(), "Lock info: Read: 0 held, 0 waiting. Write: 1 held, 0 waiting.");
  test.start(4, -1,  true,  State::Running); // Thread 2 -- wait for global write lock
  TEST_EQUAL(locker->toString(), "Lock info: Read: 0 held, 0 waiting. Write: 1 held, 1 waiting.");
  test.complete(3);
  TEST_EQUAL(locker->toString(), "Lock info: Read: 0 held, 0 waiting. Write: 1 held, 0 waiting.");
  test.complete(4);
  TEST_EQUAL(locker->toString(), "Lock info: Idle");
}

class Register
{
public:
  Register()
  {
    // Disabled by default because on failure these tend to hang or
    // crash. And there might be failures due to race conditions in
    // the tests themselves.
    static bool enable_ = false; // Avoid "unused" warning.
    if (enable_) {
      register_test("locker.simple",           test_simple);
      register_test("locker.mixed",            test_mixed);
      register_test("locker.global",           test_global);
    }
  }
} dummy;

} // namespace
