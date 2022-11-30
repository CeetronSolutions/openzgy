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

#include <condition_variable>
#include <functional>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include "declspec.h"

namespace InternalZGY {
#if 0
}
#endif

class ISDGenericDataset;

/**
 * Maintain block level read- and write locks for an open file.
 */
class OPENZGY_TEST_API Locker
{
private:
  struct Entry
  {
    int readlocks;
    int writelocks;
    std::int64_t writesize;
  };
  typedef std::map<std::int64_t, Entry> locks_t;
  std::mutex mutex_;
  std::condition_variable cv_;
  int global_readlocks_;
  int global_writelocks_;
  locks_t locks_;
  int waiting_readers_; // Only for debugging
  int waiting_writers_; // Only for debugging
  bool discarding_;
  std::int64_t highwater_;
  std::function<bool(int, const std::string&)> logger_;
  std::mutex saved_exception_mutex_;
  std::exception_ptr saved_exception_;

private:
  std::pair<int,int> countLocks(std::int64_t blocknum) const;
  int countAllWriteLocks() const;
  void updateLock(std::int64_t blocknum, int rlocks, int wlocks, std::int64_t writesize);

public:
  explicit Locker(int highwater, const std::function<bool(int, const std::string&)>& logger);
  Locker(const Locker&) = delete;
  Locker(const Locker&&) = delete;
  Locker& operator=(const Locker&) = delete;
  Locker& operator=(const Locker&&) = delete;
  std::mutex& mutex();
  void discard();
  bool setReadLock(std::unique_lock<std::mutex>& lck, std::int64_t blocknum);
  bool setWriteLock(
       std::unique_lock<std::mutex>& lck,
       std::int64_t blocknum,
       std::int64_t blocksize,
       bool also_if_discarding);
  void releaseReadLock(std::int64_t blocknum);
  void releaseWriteLock(std::int64_t blocknum);
  void runWithReadLock(std::int64_t blocknum, const std::function<void()>& fn);
  template<typename T> T runWithReadLockT(
       std::int64_t blocknum,
       const std::function<T()>& fn);
  bool runWithWriteLock(
       std::int64_t blocknum,
       std::int64_t blocksize,
       bool also_if_discarding,
       const std::function<void()>& fn);
  std::string toString() const;
  std::string toString(std::int64_t blocknum) const;
  void setException(std::exception_ptr ex);
  void throwPendingException();
};

/**
 * As runWithReadLock but runWithReadLockT<> is used for an inner
 * function that returns a value. Which then gets returned to our
 * caller.
 */
template<typename T>
T Locker::runWithReadLockT(std::int64_t blocknum, const std::function<T()>& fn)
{
  T result;
  runWithReadLock(blocknum, [&result,fn]() {result = fn();});
  return result;
}

} // namespace
