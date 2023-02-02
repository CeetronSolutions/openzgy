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

#include "locker.h"


#include "sdinterface.h"
#include "environment.h"
#include "exception.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <memory.h>

namespace InternalZGY {
#if 0
}
#endif

Locker::Locker(int highwater, const ISDGenericDataset::LoggerFn& logger)
  : mutex_()
  , cv_()
  , global_readlocks_(0)
  , global_writelocks_(0)
  , locks_()
  , waiting_readers_(0)
  , waiting_writers_(0)
  , discarding_(false)
  , highwater_(highwater)
  , logger_(logger)
  , saved_exception_mutex_()
  , saved_exception_()
{
}

/**
 * PRIVATE:
 * Return the total number of locks and total number of write locks
 * pertaining to this block number. Or to any block in the system
 * if blocknum < 0. Both explicit and global locks are counted. The
 * mutex must be held by the caller.
 *
 * If result.first == 0, a write lock request does not need to wait.
 * If result.second == 0, a read lock request does not need to wait.
 */
std::pair<int,int>
Locker::countLocks(std::int64_t blocknum) const
{
  int rlocks = global_readlocks_;
  int wlocks = global_writelocks_;
  if (blocknum < 0) {
    for (const auto& it : locks_) {
      rlocks += it.second.readlocks;
      wlocks += it.second.writelocks;
    }
  }
  else {
    const auto it = locks_.find(blocknum);
    if (it != locks_.end()) {
      rlocks += it->second.readlocks;
      wlocks += it->second.writelocks;
    }
  }
  return std::make_pair(rlocks + wlocks, wlocks);
}

/**
 * PRIVATE:
 * Count the total number of block level write locks in the system.
 * This typically represents the number of busy worker threads,
 * which may need to be capped.
 */
int
Locker::countAllWriteLocks() const
{
  int result{0};
  for (const auto& it : locks_)
    result += it.second.writelocks;
  return result;
}

/**
 * PRIVATE:
 * Set, reset, or change a lock.
 * The mutex must be held by the caller.
 * The caller is responsible for checking that the requested change
 * is legal. Any consistency checks done here are to be considered
 * as asserts.
 */
void
Locker::updateLock(
     std::int64_t blocknum,
     int rlocks,
     int wlocks,
     std::int64_t writesize)
{
  if (blocknum < 0) {
    global_readlocks_ += rlocks;
    global_writelocks_ += wlocks;
  }
  else {
    locks_t::iterator it =
      locks_.insert(std::make_pair(blocknum, Entry{0,0,0})).first;
    int rl = it->second.readlocks + rlocks;
    int wl = it->second.writelocks + wlocks;
    if (rl < 0 || wl < 0 || wl > 1 || (rl > 0 && wl > 0))
      throw OpenZGY::Errors::ZgyInternalError("Internal error setting a lock.");
    it->second.readlocks += rlocks;
    it->second.writelocks += wlocks;
    if (wlocks > 0)
      it->second.writesize = writesize;
    if (it->second.readlocks <= 0 && it->second.writelocks <= 0)
      locks_.erase(it);
  }
}

/**
 * Expose the primary mutex. This is needed for code using the lower
 * level set* and release* methods that assume the lock is already held.
 */
std::mutex&
Locker::mutex()
{
  return mutex_;
}

/**
 * Set the system in "discarding" mode where writes become no-ops and
 * reads will fail.
 */
void
Locker::discard()
{
  std::unique_lock<std::mutex> lck(mutex_);
  discarding_ = true;
  cv_.notify_all();
}

/**
 * Set a read lock. Caller must hold the mutex, and caller must
 * guarantee that a lock that was acquired must eventually be
 * released even if exceptions were encountered. Prefer using
 * runWithReadLock instead. If the function returns false the file
 * is in an unusable state, the lock was not acquired and thus
 * should not be freed either.
 */
bool
Locker::setReadLock(std::unique_lock<std::mutex>& lck, std::int64_t blocknum)
{
  ++waiting_readers_;
  while (!discarding_ && countLocks(blocknum).second > 0)
    cv_.wait(lck);
  --waiting_readers_;
  if (!discarding_)
    updateLock(blocknum, 1, 0, 0);
  return !discarding_;
}

/**
 * As setReadLock() but for writes. Use runWithWriteLock() instead
 * where possible. I.e. everywhere that background threads are not involved.
 */
bool
Locker::setWriteLock(
     std::unique_lock<std::mutex>& lck,
     std::int64_t blocknum,
     std::int64_t blocksize,
     bool also_if_discarding)
{
  ++waiting_writers_;
  bool should_run{true}, wait_thread{false}, wait_lock{false};
  for (;;) {
    should_run = !discarding_ || also_if_discarding;
    wait_thread = countAllWriteLocks() >= highwater_;
    wait_lock = countLocks(blocknum).first > 0;
    if (!should_run || (!wait_thread && !wait_lock))
      break;
    cv_.wait(lck);
  }
  --waiting_writers_;
  if (should_run)
    updateLock(blocknum, 0, 1, blocksize);
  return should_run;
}

/**
 * This must always be called after a matching setReadLock() returned true.
 */
void
Locker::releaseReadLock(std::int64_t blocknum)
{
  updateLock(blocknum, -1, 0, 0);
  cv_.notify_all();
}

/**
 * This must always be called after a matching setWriteLock() returned true.
 */
void
Locker::releaseWriteLock(std::int64_t blocknum)
{
  updateLock(blocknum, 0, -1, 0);
  cv_.notify_all();
}

/**
 * Run the required inner function protected by a read lock.
 *
 * Throws an exception if the file is or becomes marked as "discard"
 * before the inner function starts executing. This is a mechanism
 * to clean up pending reader threads after an application has
 * decided that the file is corrupted. The application will normally
 * catch and ignore these errors at some level.
 */
void
Locker::runWithReadLock(std::int64_t blocknum, const std::function<void()>& fn)
{
  std::unique_lock<std::mutex> lck(mutex_);
  if (!setReadLock(lck, blocknum)) {
    throw OpenZGY::Errors::ZgyAborted("Operation was canceled, file is unusable.");
  }
  try {
    lck.unlock();
    fn();
    lck.lock();
    releaseReadLock(blocknum);
  }
  catch (...) {
    lck.lock();
    releaseReadLock(blocknum);
    throw;
  }
}

/**
 * As runWithReadLock but for exclusive a.k.a. write access.
 * Unlike the read case a "discard" state will just cause the
 * function to return false without calling the inner function,
 */
bool
Locker::runWithWriteLock(
     std::int64_t blocknum,
     std::int64_t blocksize,
     bool also_if_discarding,
     const std::function<void()>& fn)
{
  std::unique_lock<std::mutex> lck(mutex_);
  if (!setWriteLock(lck, blocknum, blocksize, also_if_discarding))
    return false;
  try {
    lck.unlock();
    fn();
    lck.lock();
    releaseWriteLock(blocknum);
  }
  catch (...) {
    lck.lock();
    releaseWriteLock(blocknum);
    throw;
  }
  return true;
}

/**
 * Return a summary of the lock state as a human readable string.
 */
std::string
Locker::toString() const
{
  std::pair<int,int> counts = countLocks(-1); // returns (read+write, write)
  if (!counts.first && !counts.second &&
      !waiting_readers_ && !waiting_writers_)
  {
    return "Lock info: Idle";
  }
  else
  {
    std::stringstream ss;
    ss << "Lock info: "
       << "Read: " << counts.first - counts.second << " held, "
       << waiting_readers_<< " waiting. "
       << "Write: " <<counts.second << " held, "
       << waiting_writers_ << " waiting.";
    return ss.str();
  }
}

/**
 * Return the lock state for one block a human readable string.
 * Does not report waiting threads, because that information
 * is not maintained per brick.
 */
std::string
Locker::toString(std::int64_t blocknum) const
{
  std::pair<int,int> counts = countLocks(blocknum);
  if (!counts.first && !counts.second)
    return "Idle";
  else if (counts.first == 1 && counts.second == 1)
    return "Write lock";
  else if (counts.first == 1 && counts.second == 0)
    return "Read lock";
  else if (counts.first > 1 && counts.second == 0)
    return "Read lock (" + std::to_string(counts.first) + ")";
  else {
    std::stringstream ss;
    ss << "Corrupted lock: "
       << "Read: " << counts.first - counts.second << ", "
       << "Write: " << counts.second << ".";
    return ss.str();
  }
}

/**
 * Save an exception thrown in some background thread so that it can
 * be reported lated. Technically this fumctionality should have been
 * in a class of its own. But I am lazy.
 */
void
Locker::setException(std::exception_ptr ex)
{
  std::lock_guard<std::mutex> lck(saved_exception_mutex_);
  if (!saved_exception_)
    saved_exception_ = ex;
}

/**
 * Rethrow an exception that happened in a background thread.
 * Some exceptions may be lost; the logic only guarantees
 * that one of them will be reported.
 */
void
Locker::throwPendingException()
{
  std::exception_ptr ex;
  {
    std::lock_guard<std::mutex> lck(saved_exception_mutex_);
    std::swap(ex, saved_exception_);
  }
  if (ex)
  {
    std::rethrow_exception(ex);
  }
}

} // namespace
