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

#include "workorder.h"
#include "fancy_timers.h"
#include "mtguard.h"
#include "environment.h"

#include <functional>
#include <atomic>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cstdint>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>
#include <sstream>

namespace InternalZGY {
#if 0
}
#endif

class WorkOrder;
class WorkOrderOMP;
class WorkOrderTP;
class WorkOrderQueueTP;
class SummaryPrintingTimerEx;

/**
 * Helper for executing a loop, possibly multi threaded.
 * Exceptions may be collected and the first exception,
 * if any, be re-thrown after the loop is done.
 *
 * The trivial WorkOrder base class never does multi threading.
 */
class WorkOrder
{
public:
  typedef std::function<void(std::int64_t)> workstep_fn;

protected:
  std::int64_t total_;
  workstep_fn  workstep_;

public:
  WorkOrder(std::int64_t total, const workstep_fn& workstep);

  void run();
};

/**
 * Helper for executing a loop, possibly multi threaded.
 * Exceptions may be collected and the first exception,
 * if any, be re-thrown after the loop is done.
 *
 * OpenMP flavor. Leverages the existing class MtGuard for
 * handling exceptions.
 *
 * Using this helper is a bit different from #pragma omp
 * because the code to be executed is put inside a functor
 * instead of written inline. The loop iterator is passed
 * as the first argument to the functor. Calling omp_
 * functions directly is discouraged. Those will no longer
 * work if the code is switched to a different MT method.
 */
class WorkOrderOMP : public WorkOrder
{
public:
  WorkOrderOMP(std::int64_t total, const workstep_fn& workstep);
  void run();
};

/**
 * See WorkOrderTP::WorkOrderTP for a description:
 * \copydoc WorkOrderTP::WorkOrderTP
 */
class WorkOrderTP : public WorkOrder
{
private:
  std::int64_t              owner_id_; // Currently just for debugging.
  static std::atomic<std::int64_t> next_owner_id_;
  std::atomic<std::int64_t> started_; // fetch_add to allocate a new step
  std::atomic<std::int64_t> finished_; // only increment if allocated.
  std::atomic<int>          errors_; // if >0, just gobble up rest of steps.
  std::exception_ptr        first_ex_;
  std::mutex                order_done_mutex_;
  std::condition_variable   order_done_cv_;
  bool                      order_done_;

public:

  WorkOrderTP(std::int64_t total, const workstep_fn& workstep);

  ~WorkOrderTP();

  void fail();

  void rethrow(bool swallow);

  bool run1(bool owner, std::int64_t my_task = -1);

  bool more() const;

  static void runStatic(std::int64_t count, const workstep_fn& workstep);

  void run(
       const std::shared_ptr<WorkOrderTP>& self,
       std::shared_ptr<WorkOrderQueueTP> queue);
};

/**
 * See WorkOrderQueueTP::WorkOrderQueueTP for a description:
 * \copydoc WorkOrderQueueTP::WorkOrderQueueTP
 */
class WorkOrderQueueTP
{
private:
  std::mutex queue_mutex_;
  std::condition_variable queue_cv_;
  std::vector<std::shared_ptr<WorkOrderTP>> queue_;
  bool started_;
  std::atomic<int> shutdown_;
  std::atomic<int> debug_pool_size_;
  std::atomic<int> debug_busy_threads_;
private:
  std::shared_ptr<WorkOrderTP> getWorkOrder();
  std::shared_ptr<WorkOrderTP> waitForOrder();
public:
  WorkOrderQueueTP();
  ~WorkOrderQueueTP();
  static void startup(const std::shared_ptr<WorkOrderQueueTP>& queue);
  void push_back(const std::shared_ptr<WorkOrderTP>& item);
  void remove(const std::shared_ptr<WorkOrderTP>& item);
  void shutdown();
  std::pair<int, int> debugPoolSize()
  {
    // +1 to account for contribution from the thread leader.
    return std::pair<int, int>(debug_pool_size_ + 1, debug_busy_threads_ + 1);
  }
  int debugThreadsInUse()
  {
    // +1 to account for contribution from the thread leader.
    return debug_busy_threads_ + 1;
  }
  static void workerLoop(std::shared_ptr<WorkOrderQueueTP> queue);
  static std::shared_ptr<WorkOrderQueueTP> instance();
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Main entry point. Explicit provider, explicit thread count.
 *
 * \internal
 *
 * This is implemented as a wrapper around the version that will use
 * as many threads as available. It adds only chunking of steps.
 * Performance is *slightly* slower than making the entire implementation
 * aware of the chunking. Another minor issue is that when using the
 * OpenMP provider, the behavior with respect to chunking might be
 * slightly different.
 */
void
WorkOrderRunner::parallelFor(
     Provider provider,
     std::int64_t count,
     std::int64_t threadcount,
     const workstep_fn& workstep)
{
  threadcount = std::min(threadcount, count);
  const int chunksize = (int)((count + threadcount - 1) / threadcount);
  if (false && chunksize <= 1) {
     // Could disable this shortcut for better testing
    parallelFor(provider, count, workstep);
  }
  else {
    parallelFor(provider, threadcount, [&](std::int64_t num)
      {
        for (std::int64_t ii = num * chunksize; ii < (num+1) * chunksize && ii < count; ++ii)
          workstep(ii);
      });
  }
}

/**
 * Main entry point. Explicit provider, default thread count.
 */
void
WorkOrderRunner::parallelFor(
     Provider provider,
     std::int64_t count,
     const workstep_fn& workstep)
{
  // Short cut for when MT is not possible. Can be commented
  // out to get better test coverage of corner cases.
  //if (count <= 1)
  //  provider = Provider::SERIAL;

  switch (provider) {
  case Provider::SERIAL:
  {
    WorkOrder(count, workstep).run();
    break;
  }
  case Provider::OMP:
  {
    WorkOrderOMP(count, workstep).run();
    break;
  }
  case Provider::TP:
  {
    WorkOrderTP::runStatic(count, workstep);
    break;
  }
  default:
    throw std::runtime_error("Invalid choice of MT provider.");
  }
}

/**
 * Main entry point. Default provider, explicit thread count.
 */
void
WorkOrderRunner::parallelFor(
     std::int64_t count,
     std::int64_t threadcount,
     const workstep_fn& workstep)
{
  parallelFor(defaultProvider(), count, threadcount, workstep);
}

/**
 * Main entry point. Default provider, default thread count.
 */
void
WorkOrderRunner::parallelFor(
     std::int64_t count,
     const workstep_fn& workstep)
{
  parallelFor(count, count, workstep);
}

WorkOrderRunner::Provider
WorkOrderRunner::defaultProvider()
{
  static std::string str = InternalZGY::Environment::getStringEnv("OPENZGY_MT_PROVIDER", "TP");
  if (str == "SERIAL") return Provider::SERIAL;
  if (str == "OMP") return Provider::OMP;
  if (str == "TP") return Provider::TP;
  throw std::runtime_error("OPENZGY_MT_PROVIDER can be SERIAL, OMP, or TP.");
}

std::pair<int, int>
WorkOrderRunner::debugPoolSize()
{
  std::pair<int, int> result = debugPoolSize(defaultProvider());
  // Even when the thread pool is selected, some loops still
  // use OpenMP directly. Heuristic: If it looks like no TP loop
  // is active then look for an OpenMP one.
  if (result.second == 1 && defaultProvider() == Provider::TP)
    result = debugPoolSize(Provider::OMP);
  return result;
}

/**
 *
 * Estimate how many threads of any type are currently in use.
 *
 * Timers can use this to guess how long the task might use in wall clock time.
 * Under the assumption that the task scales perfectly with number of threads
 * and all current threads, both TP and OpenMP, are busy with this particular
 * task. Looking at both the estimate and the measured wall clock time outside
 * the MT loop might hint at how much overhead there is. The numbers might be
 * completely off if there are multiple unrelated tasks running multithreaded.
 */
int
WorkOrderRunner::debugThreadsInUse()
{
  return WorkOrderQueueTP::instance()->debugThreadsInUse() + omp_get_num_threads() - 1;
}

/**
 * Debug only! May be used for timers.
 *
 * Estimated how many threads are available and how many are in use.
 * Behavior inside nested loops or active loops in multiple threads
 * is undefined. For OpenMP, "in use" is the number of threads in
 * the current team, also counting idle threads waiting to join.
 * For TP, "in use" is the global number of busy threads.
 * This can be > 1 in the calling thread even of the MT code is
 * running somewhere else. Also the TP pool size will show as 1
 * until the first time it is used.
 */
std::pair<int, int>
WorkOrderRunner::debugPoolSize(Provider provider)
{
  switch (provider) {
  case Provider::SERIAL:
    return std::pair<int, int>(1, 1);
  case Provider::OMP:
    return std::pair<int, int>
      (omp_get_max_threads(), omp_get_num_threads());
  case Provider::TP:
      return WorkOrderQueueTP::instance()->debugPoolSize();
  default:
    throw std::runtime_error("Invalid choice of MT provider.");
  }
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

WorkOrder::WorkOrder(std::int64_t total, const workstep_fn& workstep)
    : total_(total)
    , workstep_(workstep)
{
}

/**
 * Execute the "parallel for" loop described by *this.
 * In the base class, no actual multi-threading will be used.
 */
void
WorkOrder::run()
{
  for (std::int64_t ii = 0; ii < total_; ++ii) {
    workstep_(ii);
  }
}

/////////////////////////////////////////////////////////////////////////////
///   WorkOrderOMP   ////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

WorkOrderOMP::WorkOrderOMP(
     std::int64_t total,
     const workstep_fn& workstep)
  : WorkOrder(total, workstep)
{
}

void
WorkOrderOMP::run()
{
  InternalZGY::MTGuard guard("WorkOrderOMP",
    total_ > std::numeric_limits<int>::max() ?
    std::numeric_limits<int>::max() :
    static_cast<int>(total_));
#pragma omp parallel for
  for (std::int64_t ii = 0; ii < total_; ++ii) {
    guard.run([&](){
	workstep_(ii);
      });
  }
  guard.finished();
};

/////////////////////////////////////////////////////////////////////////////
///   WorkOrderTP   /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

std::atomic<std::int64_t> WorkOrderTP::next_owner_id_(1);

/**
 * Helper for executing a loop, possibly multi threaded.
 * Exceptions may be collected and the first exception,
 * if any, may be re-thrown after the loop is done.
 *
 * \internal
 *
 * This version uses a home-grown thread pool implementation.
 * Also use some code from inside class MtGuard.
 * It isn't practical to use that class as-is.
 *
 * The team leader (the thread that started the parallel for):
 *
 *   - Describe the loop using a WorkOrderTP instance. Only one instance
 *     is created. This instance represents all iterations of this loop.
 *
 *   - Add the WorkOrderTP instance to the global WorkOrderQueueTP list
 *     to make it visible to workers in the thread pool.
 *
 *   - Fetch one step (one iteration) at a time from the WorkOrderTP,
 *     and call the workstep function to process this step.
 *     This will eventually complete the work even if no thread pool
 *     worker is available to help.
 *     The team leader only executes steps belonging to its own loop.
 *     It doesn't look at the global global WorkOrderQueueTP list
 *
 *   - When the last step has been started so no more can be fetched,
 *     wait on a condition variable belonging to the WorkOrderTP.
 *     The loop is finished when all iterations have not just started
 *     but also finished.
 *
 * Each background worker:
 *
 *   - Find a WorkOrderTP instance in the global WorkOrderQueueTP list
 *     where not all steps has been started.
 *     Wait on a condition variable in WorkOrderQueueTP if there are none.
 *
 *   - Fetch one step (one iteration) from this instance at a time,
 *     and call the workstep function to process this step.
 *     The worker can execute steps from any work order.
 *
 *   - When the last step has been started in this work order, go to
 *     the top and try to find another work order. Wait if there is
 *     currently no more work available. At that point, the only thing
 *     that can provide more work is that some code has added a new
 *     WorkOrderTP instance to the global WorkOrderQueue.
 *
 * No fair scheduling:
 *
 *     If there are more than one active loop, either started from different
 *     team leaders or due to a nested call, then the steps in the first loop
 *     have precedence. It would also be possible to do some round robin
 *     scheduling. I doubt that this is a good idea though.
 *
 * Scheduling of steps
 *
 *   - If the thread pool is disabled, the steps run sequentially.
 *
 *   - When the thread pool has N idle threads, the first N iterations
 *     will have static scheduling with each thread running exactly
 *     one step. If there are more iterations than that, the rest
 *     are scheduled dynamically. They get picked up by worker threads
 *     as the workers become available.
 *
 *   - The above isn't entirely accurate because of race conditions.
 *     E.g. if the steps complete very quickly then some worker threads
 *     might be able to grab more than one step because some other
 *     worker was slow to wake up. Just like myself in the morning.
 *     This should not be noticeable.
 *
 *   - The application cannot directly control the number of threads
 *     used to run a loop. The actual number of threads will be
 *     min(total,available_threads). If the application believes this
 *     will result in too small granularity then it needs to code a
 *     double loop with the count of the outer loop being the desired
 *     number of threads. Take care if the steps can take very
 *     different times to complete. Static scheduling might make
 *     the loop less than optimal.
 *
 *   - With static scheduling the progress indicator and the cancel
 *     feature are not useful. Progress will at most be reported twice:
 *     when the team leader finishes its task and when all threads
 *     are done. Even if the first of those returns "cancel", all the
 *     other threads will run to completion.
 *
 * The condition variable in the work order:
 *
 *     The cv is order_done_cv_ and the variable is the bool order_done_.
 *     order_done_mutex_ is only used for those two, so it will be held
 *     for a very short period. The actual state change from busy to all
 *     done involves atomic variables, so there is no explicit locking
 *     there. order_done can only be set true by the process (either the
 *     team leader or a worker thread) that got the honor of executing
 *     the last step of this work order.
 *
 * The condition variable in the queue:
 *
 *     The condition being tested for is that no work order instances
 *     have any steps that are not done or at least started. Once a work
 *     order has reached this state, it cannot go back. This means that
 *     idle workers can sleep until a new work order is added to the
 *     queue. It is sufficient to protect all access to the queue itself
 *     with a mutex. Work orders retrieved from the queue are safe to use
 *     even if they later get removed from the queue. Because they are
 *     smart pointers.
 *
 * There are a number of race conditions, each of which I believe are harmless.
 *
 *     When a worker finds a work order with more() returning true, it
 *     might see run1() returning false because some other process
 *     snapped up the step between the calls to more() and run1(). The
 *     worker's attempt to execute a step becomes a no-op. The next call
 *     to more() will return false. Which means we are back on track. One
 *     minor note: started_ is used to allocate a step to be processed.
 *     The race condition described here will cause started_ to become
 *     greater than the total number of steps. Both started_ == total_
 *     and started_ > total mean the same thing: No more iterations need
 *     to be started. If you want the actual number of steps that have at
 *     least been started, just use std::min(started_, total_).
 *
 *     When a worker catches an exception, it initially believes that
 *     this was the first exception seen. Because otherwise the step
 *     would not have been executed at all. In the catch block, however,
 *     it could find out that some other thread has already reported a
 *     "first exception". No problem; just swallow this exception.
 *
 * Deadlock analysis:
 *
 *   - The workstep code is responsible for not causing deadlocks itself,
 *     and for properly synchronizing shared data.
 *
 *   - No locks are held when executing the worksteps. This means that
 *     the code executing a step is free to place locks that can't
 *     cause deadlocks by themselves.
 *
 *   - For the locks owned by this module, only one lock will be held
 *     at a time, and only trivial code will be executed while the build
 *     is held.
 *
 * Issues with object lifetime:
 *
 *   - Lifetime of WorkOrderQueueTP: This is a singleton in a
 *     static variable, so it will get destructed during application
 *     shutdown. This might be problematic because threads in the pool
 *     might still be active. Mitigated by using a smart pointer that
 *     each worker thread holds a reference to.
 *
 *   - Lifetime of WorkOrderTP: Instances might not be destructed
 *     immediately when the "parallel for" loop is done. But I don't
 *     think this will be a problem. At the end of runStatic(), all
 *     the work described in the work order is done. But smart pointer
 *     references to the work order might still exist in worker
 *     threads executing something in this work order. And (but maybe
 *     not) in the global queue. The tasks should complete in
 *     microseconds. And this happens completely in the background.
 *
 *     Exiting the application when it is not idle (i.e. executing
 *     a loop in the background) might cause a crash. I suspect this
 *     is true for a lot of other code as well.
 *
 * Notes:
 *
 *   - Code that accesses started_++ and gets back a value between 0
 *     and total_-1 is always responsible for doing a corresponding
 *     finished_++. Even if there was an exception. If the increment
 *     is not done then the loop will hang.
 *
 *   - Code that accesses finished_++ and gets back exactly total_-1 is
 *     responsible for marking the work order as complete. This means
 *     setting order_done_ and waking up the team leader, if needed.
 *     It is possible that the team leader itself is the thread that
 *     got the honor of being the last one finishing a step.
 *
 *   - The thread that starts the last (i.e. total_-1) step won't
 *     necessarily be the same thread that marks the work order as
 *     complete.
 *
 * Risks:
 *
 *     Unfortunately there is a real risk of data races and/or deadlocks
 *     that I have not spotted. And issues with object lifetime. This kind
 *     of bug can be really difficult to spot in testing.
 *
 * ------------------------------------------------------------------------
 *
 * Checklist for thread safe access to data members. This is all
 * documented above, but the text below makes a complete list of data
 * members and then explains why they are supposed to be safe.
 *
 * Most data members in the WorkOrder base class are read only, i.e. only
 * set in the constructor. When it comes to invoking methods on these:
 *
 *   int64_t      total_       -- Scalar, no methods.
 *   workstep_fn  workstep_    -- Caller is required to make this thread safe.
 *
 * Member variables in WorkOrderTP:
 *
 *   int64_t owner_id_         -- Read only.
 *   int64_t next_owner_id_    -- Atomic
 *   int64_t started_          -- Atomic
 *   int64_t finished_         -- Atomic
 *   int     errors_           -- Atomic
 *   exception_ptr first_ex    -- Owned by he who changes errors_ from 0 to 1.
 *   mutex              order_done_mutex_ -- N/A
 *   condition_variable order_done_cv_    -- Protected by order_done_mutex_
 *   bool               order_done_       -- Protected by order_done_mutex_
 *
 * Member variables in WorkOrderQueueTP
 *
 *   All data members (queue_cv_ and queue_) are protected by queue_mutex_.
 *   waitForOrder(), push_back(), remove all() all lock queue_mutex_
 *   getWorkOrder() requires queue_mutex_ to be held on entry.
 *   instance() has a static initializer.
 *   workerLoop() is trickier:
 *   workerLoop() calls waitForOrder() which does set a lock.
 *   workerLoop() calls WorkOrderTP::run1() which should protect its own data.
 *
 * Checklist for data lifetime. This also duplicates explanations above.
 *
 * The thread pool is permanent. And even if it were to be shut down,
 * there are no outside references to it. So this should not have any
 * issues.
 *
 * The list inside WorkOrderQueueTP contains smart pointers. So, even
 * when the list gets emptied, any WorkOrder instances extracted from
 * the list will remain valid. The instance is a singleton with a
 * static constructor. It will be destructed on application exit.
 * I might consider deliberately leaking it instead.
 *
 * References to the WorkOrder created inside runStatic() might linger
 * in background processes, so the instance might not be destructed
 * immediately when the function ends. Those instances will all be
 * flagged as done. So they will not try to execute any steps. And
 * they will disappear quickly. Lingering references in the
 * WorkOrderQueueTP shouldn't be possible
 */
WorkOrderTP::WorkOrderTP(
     std::int64_t total,
     const workstep_fn& workstep)
  : WorkOrder(total, workstep)
  , owner_id_(next_owner_id_.fetch_add(1))
  , started_(0)
  , finished_(0)
  , errors_(0)
  , first_ex_()
  , order_done_mutex_()
  , order_done_cv_()
  , order_done_(false)
{
}

WorkOrderTP::~WorkOrderTP()
{
  // Be very careful putting any code here. The destructor might be
  // called later than you think.

  // The next line will swallow any unreported exception, after
  // logging it with logger_(). So, logger_ must remain valid to
  // call even when the application is shutting down.
  rethrow(true);
}

/**
 * To be called from inside a catch statement. Stores the exception
 * so it can be rethrown outside the parallel region.
 *
 * Normally the throw/catch will not execute at all if _errors != 0,
 * but a harmless race condition might give us _errors > 0 and this
 * code will swallow the exception
 */
void
WorkOrderTP::fail() {
  if (errors_.fetch_add(1) == 0) {
    first_ex_ = std::current_exception();
    if (!first_ex_)
      first_ex_ = std::make_exception_ptr
        (std::runtime_error("fail() with no current exception."));
  }
}

void
WorkOrderTP::rethrow(bool swallow)
{
  if (errors_.load() != 0 && first_ex_) {
    std::exception_ptr throwme = first_ex_;
    first_ex_ = std::exception_ptr();
    if (!swallow) {
      std::rethrow_exception(throwme);
    }
    else {
      try {
        std::rethrow_exception(throwme);
      }
      catch (const std::exception& ex) {
        std::cerr << "EXCEPTION inside parallel loop: " << ex.what() << std::endl;
      }
      catch (...) {
        std::cerr << "EXCEPTION inside parallel loop: " << std::endl;
      }
    }
  }
}

/**
 * \brief Run a single step of a single loop.
 *
 * \details Return false if there are no more steps that need to be
 * started in this loop. A false return value can be trusted. A true
 * return does not guarantee that there is more work. Because of
 * some harmless race conditions.
 *
 * Note, the function could also have returned true in the sense of
 * "a step was in fact executed" instead of "there are more steps".
 * That might be more orderly.
 */
bool
WorkOrderTP::run1(bool owner, std::int64_t my_task)
{
  if (my_task < 0)
    my_task = started_.fetch_add(1);
  if (my_task >= total_) {
    // There was a harmless race condition,
    // and somebody else gets the honor of
    // signalling order_filled. started_ is now too high, but nobody cares.
    return false;
  }
  if (errors_.load() == 0) {
    try {
      workstep_(my_task);
    }
    catch (...) {
      fail();
    }
  }
  if (finished_.fetch_add(1) == total_ - 1) {
    std::unique_lock<std::mutex> lck(order_done_mutex_);
    order_done_ = true;
    order_done_cv_.notify_all();
    return false; // Return no-more-left (false) or did-something (true) ?
  }
  return true;
}

/**
 * \brief Return true if there might be more steps waiting to start
 * in this WorkOrder.
 *
 * \details Emphasis on "might". The caller must be prepared for
 * "no more work" because of race conditions. This is by design.
 */
bool
WorkOrderTP::more() const
{
  // Note that started_ is actually "attempted started", it might be
  // greater than requested_
  return started_.load() < total_;
}

void
WorkOrderTP::runStatic(
     std::int64_t       count,
     const workstep_fn& workstep)
{
  auto wo = std::make_shared<WorkOrderTP>(count, workstep);
  std::shared_ptr<WorkOrderQueueTP> queue = WorkOrderQueueTP::instance();
  WorkOrderQueueTP::startup(queue);
  wo->run(wo, queue);
  // References to "wo" might linger in background processes
  // so the instance we created might not be destructed immediately.
  // Those instances will all be flagged as done. So they will
  // disappear quickly.
}

/**
 * \brief Main loop to be executed by the team leader.
 * \details For testing, it is allowed to pass nullptr as the queue.
 * There will then be no background processing.
 */
void
WorkOrderTP::run(const std::shared_ptr<WorkOrderTP>& self, std::shared_ptr<WorkOrderQueueTP> queue)
{
  // Check for degenerate case, needed because if count is 0
  // there will be no thread finishing the last step,
  // and the thread leader would wait indefinitely.
  if (total_ <= 0)
    return;

  if (total_ > 1 && queue != nullptr) {
    queue->push_back(self);
  }

  // Steps might already be executing in the thread pool.

  // See below. Iteration 0 could be  pre-allocated to the thread
  // calling run().
  //run1(true, 0);

  while (run1(true))
  {
  }

  // run1() returned false meaning all is done or in flight.
  // order_done_ might already been set by us, if our thread
  // got the honor of incrementing finished_ to total_.
  // Or it might have been set by somebody else, or it might
  // still be pending. Ideally a unit test should be able to
  // exercise all these cases. And don't forget to test the
  // total_ == 1 special case.
  {
    std::unique_lock<std::mutex> lck(order_done_mutex_);
    while (!order_done_)
      order_done_cv_.wait(lck);
  }

  if (total_ > 1 && queue != nullptr) {
    // Throws if not found.
    queue->remove(self);
  }

  // Steps might still be attempted executed from the thread pool.
  // But the workers won't find any work to be done. Proper use
  // of smart pointers is crucial to avoid lifetime issues.

  // Delayed throw of first exception seen inside loop, if any.
  rethrow(false);

  // Suggestion for a nice to have feature: Iteration 0 will always
  // be executed by the thread that called run(). This makes it
  // slightly more deterministic if the code has special handling
  // that only needs to run once and has been coded with if (it==0).
  //
  // Initialize  started_ to 1, assuming the team leader will eventually
  // get around to doing #0. And call run1(0) from the starting thread
  // before looping to pick up more jobs.

#if 0
  // DEBUG
  // Expect false positives for this consistency check. There is a race
  // between the thread having performed its duties and the thread
  // registers itself as available. There will also be warnings if
  // there are nested loops or loops started from unrelated threads.
  auto state = WorkOrderRunner::debugPoolSize();
  if (state.second != 1) {
    std::stringstream ss;
    ss << "Warning: MT mechanism is not idle: " << state.second << "/" << state.first << "\n";
    std::cerr << ss.str() << std::flush;
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////
///   WorkOrderQueueTP   //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Create a list of currently executing "parallel for" loops.
 * There is a corresponding thread pool to help speed up the loops.
 *
 * The list and the thread pool are currently singletons.
 * It would also be possible to have one set per open file.
 *
 * The class is only used together with WorkOrderTP, not
 * the vanilla or the OpenMP flavor.
 */
WorkOrderQueueTP::WorkOrderQueueTP()
  : queue_mutex_()
  , queue_cv_()
  , queue_()
  , started_(false)
  , shutdown_(0)
  , debug_pool_size_(0)
  , debug_busy_threads_(0)
{
}

WorkOrderQueueTP::~WorkOrderQueueTP()
{
  // If there is one pool per open file then the threads need
  // to be shut down now. If used as a singleton don't bother.
}

void
WorkOrderQueueTP::startup(const std::shared_ptr<WorkOrderQueueTP>& queue)
{
  // Function is not a normal member function because it needs the shared_ptr.
  // Default pool size is the available hardware threads, including
  // hyperthreading. Clipped to [4, 96] because extreme values are
  // likely to give worse performance. YMMV. Any count can be set using
  // OPENZGY_NUMTHREADS_POOL in the environment.
  const int poolsize = Environment::getNumericEnv
    ("OPENZGY_NUMTHREADS_POOL",
     std::max(4, std::min(96, (int)std::thread::hardware_concurrency())));

  std::unique_lock<std::mutex> lck(queue->queue_mutex_);
  if (!queue->started_) {
    queue->debug_pool_size_.fetch_add(poolsize);
    for (int ii=0; ii<poolsize; ++ii)
      std::thread(workerLoop, queue).detach();
    queue->started_ = true;
    // We probably don't need to keep references to the threads.
  }
}

/**
 * Add information about a new "parallel for" being executed.
 * There can be more than one, due to nested for loops and
 * possibly also unrelated functions already running in the
 * background e.g. via a Salmon worker.
 *
 * The threadpool is encouraged to help the thread that started
 * the "parallel for" but it is not required to do so. Currently
 * there is a constant number of worker threads available.
 * If the demand is higher then there is no attempt at fair
 * distribution of resources.
 */
void
WorkOrderQueueTP::push_back(const std::shared_ptr<WorkOrderTP>& item)
{
  std::unique_lock<std::mutex> lck(queue_mutex_);
  queue_.push_back(item);
  queue_cv_.notify_all();
}

void
WorkOrderQueueTP::remove(const std::shared_ptr<WorkOrderTP>& item)
{
  std::unique_lock<std::mutex> lck(queue_mutex_);
  auto it = std::find(queue_.begin(), queue_.end(), item);
  if (it != queue_.end())
    queue_.erase(it);
  else
    throw std::runtime_error("Object to be removed is missing");
}

/**
 * Ask the worker threads to shut down as soon as there is no more
 * work available.
 */
void
WorkOrderQueueTP::shutdown()
{
  shutdown_.fetch_add(1);
}

/**
 * Get a WorkOrder entry that might have more work left.
 * Due to race conditions it is possible that what seemed
 * like more work ended up being grabbed by some other thread.
 *
 * If the return is empty then the caller can safely assume
 * there is no more work to be done until a new WorkOrder
 * has been registered.
 *
 * The mutex must be HELD ON ENTRY.
 */
std::shared_ptr<WorkOrderTP>
WorkOrderQueueTP::getWorkOrder()
{
  for (std::shared_ptr<WorkOrderTP> it : queue_)
    if (it->more())
      return it;
  return nullptr;
}

/**
 * Get a WorkOrder entry that might have more work left.
 * Due to race conditions it is possible that what seemed
 * like more work ended up being grabbed by some other thread.
 *
 * Will not return until more work is available (non-null)
 * or the thread queue is being shut down (nullptr).
 */
std::shared_ptr<WorkOrderTP>
WorkOrderQueueTP::waitForOrder()
{
  std::shared_ptr<WorkOrderTP> wo;
  std::unique_lock<std::mutex> lck(queue_mutex_);
  while (!shutdown_.load() && (wo = getWorkOrder()) == nullptr)
    queue_cv_.wait(lck);
  return wo;
}

/**
 * \brief Main loop to be executed by threads in the pool.
 * \details
 * Logically this should have been an instance function, not a static.
 * But it needs to hold a smart pointer reference to the queue.
 * Not just the raw "this".
 */
void
WorkOrderQueueTP::workerLoop(std::shared_ptr<WorkOrderQueueTP> queue)
{
  //std::cerr << "Worker is starting up\n" << std::flush;
  std::shared_ptr<WorkOrderTP> wo;
  while ((wo = queue->waitForOrder()) != nullptr) {
    bool more = true;
    while (more) {
      queue->debug_busy_threads_.fetch_add(1);
      more = wo->run1(false, -1);
      queue->debug_busy_threads_.fetch_sub(1);
    };
  }
  //std::cerr << "Worker is shutting down\n" << std::flush;
}

/**
 * The WorkOrderQueueTP is used as a singleton, but has to be a smart
 * pointer. Because it can get destructed at application exit. While
 * there might still be worker threads holding a reference to it.
 *
 * If using a short lived instance instead of a singleton, for testing,
 * then it is more obvious that it needs to be a smart pointer.
 */
std::shared_ptr<WorkOrderQueueTP>
WorkOrderQueueTP::instance()
{
  static std::shared_ptr<WorkOrderQueueTP> instance_ =
    std::make_shared<WorkOrderQueueTP>();
  return instance_;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

} // namespace.
