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

/**
 * \file impl/file_asyncio.cpp
 *
 * This is a file access layer that adds support for fully
 * asynchronous writes. Reads are still synchronous. It is meant
 * to be inserted above SDGenericDataset, below SeismicStoreFile.
 *
 * In short, allow writes to continue in a background thread without
 * the application realizing that this happened. Fancy locking is
 * needed in case the application tries to read information from the
 * file that is in reality still being written out.
 *
 * Note that while a SDGenericDataset is supposed to be thread safe,
 * it is unspecified what happens if multiple writes or both read and
 * write is done concurrently from separate threads. Or if one thread
 * issues a close while other threads continue working on the file.

 * The code here doesn't try to change those rules. It only protects
 * against race conditions that arise because of the background
 * writes. So for a given block there can be at most one background
 * thread writing it and either one write or multiple reads, but not
 * both, waiting for that block. This implies that the per-block
 * locking mechanism doesn't need to worry about fairness.
 *
 * Additional features:
 *  - A global read- or write lock (for all existing and future blocks)
 *    can be set when needed. E.g. to support close().
 *
 *  - Write locks are also limited by a per-file cap on active threads.
 *    TODO: Fair queuing when blocking on available threads.
 *    Issue is mostly relevant when mixing reads and writes in the
 *    same area of the file.
 *
 *  - A global "discarding" flag may be set if the file is known to be
 *    corrupted anyway. This will abort all waiting locks except for
 *    close() and either notify waiting threads or make them throw.
 *
 * Future:
 *  - Buffer pools (one for each open file) for writing.
 *  - Thread pools (one for each open file) for writing.
 *  - Less aggressive locking of e.g. getBlockSize(),
 *    use information saved along with the locks instead.
 *  - Named instead of numbered blocks makes the code more general.
 *  - Reads might be satisfied from a buffer currently owned by
 *    a worker thread and in the process of being out. (*)
 *  - Reads might be satisfied from a buffer in the buffer pool
 *    that has been released but not re-used yet and not made
 *    stale by some other write. (*)
 *
 * Non-requirements:
 *
 *  - A global thread pool for all open files would give better load
 *    balancing for applications that do concurrent access of many
 *    files. But it may be too complex to implement.
 *
 *  - Allow multiple write- and delete requests from concurrent
 *    threads. The last one submitted wins. The others get canceled.
 *    Deletes would end up with with some really subtle problems.
 *    write, delete, write might due to a cancellation end up as a
 *    write, write which means that now the second call to write needs
 *    to inform sdapi that the block might already exist. And... since
 *    the order of scheduling of threads is unknown, the result is
 *    unspecified anyway.
 *
 * Reasons for using (or not) this new code in OpenZGY:
 *
 * The delayed write feature is designed to be as transparent as
 * possible to the application that uses it, but there are some
 * caveats.
 *
 *  - Errors during write may be thrown from a different write or from
 *    close. For OpenZGY this is not a big issue because the code
 *    cannot recover from a failed write anyway. So the file is toast
 *    at that point.
 *
 *  - If the same file is opened both for read and write at the same
 *    time, or multiple times for write, this will very likely lead
 *    to corrupted data. Note that even without the async write this
 *    is pretty bad practice. For OpenZGY this is not a big issue
 *    because this scenario would fail in any case. It just fails
 *    harder now.
 *
 *  - In the current implementation, misaligned writes to the cloud
 *    might get worse performance than the old "segsplit"
 *    multithreading did. (*)
 *
 * Footnote (*)
 *
 * OpenZGY might see worse performance for misaligned writes when
 * using async write instead of the older segsplit mechanism. Some
 * read/modify/write cycles might need to read back data from a region
 * that is neither very fresh and part of the current segment nor so
 * old that it has been fully flushed to disk. Not only can the data
 * not be read from the in-memory output buffer, the library will also
 * need to wait for this buffer to reach the data store before it can
 * be read back.

 * The app should only be writing from one thread at a time, so the
 * issue here is to be allowed to access buffers in the detached
 * threads. For this to be safe the reading thread might be allowed to
 * set a read lock on the block it wants, in spite of it being write
 * locked, and get a shared pointer to the buffer, before grabbing the
 * region it wants. The idea is that this special read lock does not
 * allow reading from the data store, only from the buffer. A stretch
 * goal would be to do something similar with the buffer pool that can
 * hold unused buffers that happen to contain the correct data.
 */

#include "sdinterface.h"
#include "environment.h"
#include "locker.h"
#include "../exception.h"

#ifdef HAVE_SD // Most of the file

#ifndef _WIN32
#include <SDManager.h>
#include <SDGenericDataset.h>
#else
#include <SDAPI/SDManager.h>
#include <SDAPI/SDGenericDataset.h>
#endif

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

class LogEntryAndExit
{
  std::shared_ptr<Locker> tracker_;
  ISDGenericDataset::LoggerFn logger_;
  std::string name_;
  bool enabled_;
public:
  LogEntryAndExit(
       const std::shared_ptr<Locker>& tracker,
       const ISDGenericDataset::LoggerFn& logger,
       const std::string& name)
    : tracker_(tracker)
    , logger_(logger)
    , name_(name)
    , enabled_(logger(3, ""))
  {
    std::string msg = "asio +" + name_ + " " + tracker_->toString() + "\n";
    logger_(3, msg);
  }
  ~LogEntryAndExit()
  {
    std::string flag(std::uncaught_exception() ? "*" : "-");
    std::string msg = "asio " + flag + name_ + " " + tracker_->toString() + "\n";
    logger_(3, msg);
  }
};

/**
 * This is a file access layer that adds support for fully
 * asynchronous writes. Reads are still asynchronous. It is meant
 * to be inserted above ISDGenericDataset, below SeismicStoreFile.
 */
class SDGenericDatasetAsyncWrite : public ISDGenericDataset
{
  static std::atomic<std::int64_t> nextseq_;
  std::shared_ptr<Locker> tracker_;
  std::shared_ptr<ISDGenericDataset> relay_;
  LoggerFn logger_;
  SDGenericDatasetAsyncWrite(const SDGenericDatasetAsyncWrite&) = delete;
  SDGenericDatasetAsyncWrite& operator=(const SDGenericDatasetAsyncWrite&) = delete;

private:
  ISDGenericDataset& relay() { return *relay_; }
  const ISDGenericDataset& relay() const { return *relay_; }

public:
  SDGenericDatasetAsyncWrite(
       std::shared_ptr<ISDGenericDataset> sgds,
       std::int64_t highwater,
       const LoggerFn& logger)
    : tracker_(std::make_shared<Locker>((int)highwater, logger))
    , relay_(sgds)
    , logger_(logger)
  {
  }

  // Methods from ISDGenericDataset that are always passed thru

  virtual void open(
       seismicdrive::SDDatasetDisposition disposition,
       const std::unordered_map<std::string, std::string> &args) override
  {
    relay().open(disposition, args);
  }

  virtual std::string getConsistencyID() const override
  {
    return relay().getConsistencyID();
  }

  virtual std::string getSerializedContext() const override
  {
    return relay().getSerializedContext();
  }

  virtual bool getReadonlyMode() const override
  {
    return relay().getReadonlyMode();
  }

  virtual void setReadonlyMode(bool readonly) override
  {
    LogEntryAndExit lee(tracker_, logger_, "setRoMode:  ");
    tracker_->runWithReadLock(-1, [this,readonly]()
    {
      relay().setReadonlyMode(readonly);
    });
  }

  virtual void setExponentialRetryBackoffPolicy(
       const seismicdrive::ExponentialRetryBackoffPolicy *policy,
       const seismicdrive::HttpConnectionLink link) override
  {
    return relay().setExponentialRetryBackoffPolicy(policy, link);
  }

  virtual long long getSize() const override
  {
    LogEntryAndExit lee(tracker_, logger_, "getSize:    ");
    return tracker_->runWithReadLockT<long long>(-1, [this]()
    {
      return relay().getSize();
    });
  }

  void writeBlock(int block, const char *data, std::size_t nbytes, bool check_overwrite = false) override
  {
    LogEntryAndExit lee(tracker_, logger_, "writeBlock: ");
    tracker_->throwPendingException();
    // Careful. Need to replicate the insides of runWithWriteLock()
    // because the lock will end up being released by a different thread,
    std::unique_lock<std::mutex> lck(tracker_->mutex());
    if (!tracker_->setWriteLock(lck, block, nbytes, false))
      return; // Shutting down & discarding
    try {
      lck.unlock();
      startWriteWorker(tracker_, relay_, block, data, nbytes, check_overwrite);
      lck.lock();
      // The worker will do tracker->releaseWriteLock(block);
    }
    catch (...) {
      // Unlikely, but maybe the thread couldn't start.
      // In most cases the worker thread releases the lock, not us.
      lck.lock();
      tracker_->releaseWriteLock(block);
      throw;
    }
  }

  void readBlock(int block, char *data, size_t offset, size_t nbytes) const override
  {
    LogEntryAndExit lee(tracker_, logger_, "readBlock:  ");
    tracker_->throwPendingException();
    tracker_->runWithReadLock(block, [=]()
    {
      relay().readBlock(block, data, offset, nbytes);
    });
  }

  virtual void deleteBlock(const std::string &blockName) override
  {
    LogEntryAndExit lee(tracker_, logger_, "deleteBlock:");
    tracker_->throwPendingException();
    // TODO, turn blockName into an integer.
    // Or add deleteBlock(int) to the API.
    // Then I only need to wait on that particular block.
    tracker_->runWithWriteLock(-1, 0, false, [=]()
    {
      relay().deleteBlock(blockName);
    });
  }

  long long getBlockSize(int block) const override
  {
    LogEntryAndExit lee(tracker_, logger_, "getBlockSz: ");
    // TODO minor: Could actually get it from the lock enties.
    return tracker_->runWithReadLockT<long long>(block, [this,block]()
    {
      return relay().getBlockSize(block);
    });
  }

  std::vector<long long> getBlockSizes(const std::vector<std::string> &blockNames) const override
  {
    LogEntryAndExit lee(tracker_, logger_, "getBlkSizes:");
    return tracker_->runWithReadLockT<std::vector<long long>>(-1, [this,&blockNames]()
    {
      return relay().getBlockSizes(blockNames);
    });
  }

  uint64_t getBlockNum() const override
  {
    LogEntryAndExit lee(tracker_, logger_, "getBlockNum:");
    // TODO minor: Could actually get it from the lock entries.
    // But only if blocks are numbered sequentially.
    return tracker_->runWithReadLockT<uint64_t>(-1, [this]()
    {
      return relay().getBlockNum();
    });
  }

  /**
   * TODO add a discard boolean to xx_close(). Big ripple effect.
   */
  void close() override
  {
    LogEntryAndExit lee(tracker_, logger_, "close:      ");
    tracker_->throwPendingException();
    bool ok = tracker_->runWithWriteLock(-1, 0, true, [=]()
    {
      relay().close();
    });
    if (!ok) {
      // Should not happen because I asked to wait also if discarding.
      throw OpenZGY::Errors::ZgyInternalError("Internal error, cannot lock before close");
    }
  }

private:

  static void startWriteWorker(
       std::shared_ptr<Locker> tracker,
       std::shared_ptr<ISDGenericDataset> relay,
       const std::int64_t blocknum,
       const char *data,
       std::int64_t nbytes,
       bool check_and_overwrite)
  {
    // The thread that started the write will soon get back control.
    // At that point "data" will no longer be valid.
    // TODO: A buffer pool, because most segments will be the same size.
    // TODO: Even better, pass smart pointers in the api and don't copy.
    std::shared_ptr<char> buffer(new char[nbytes]);
    memcpy(buffer.get(), data, nbytes);

    // TODO: Using a thread pool might be more efficient. One per file
    // or even one global thread pool for better throttling. But the
    // latter opens up a huge can of worms where I/O in one file might
    // need to unblock I/O on another thread.
    std::thread(&SDGenericDatasetAsyncWrite::writeWorker,
                tracker, relay, blocknum,
                buffer, nbytes, check_and_overwrite).detach();
  }

  static void writeWorker(
       std::shared_ptr<Locker> tracker,
       std::shared_ptr<ISDGenericDataset> relay,
       const std::int64_t blocknum,
       std::shared_ptr<char> data,
       std::int64_t nbytes,
       bool check_and_overwrite)
  {
    try {
      relay->writeBlock((int)blocknum, data.get(), nbytes, check_and_overwrite);
    }
    catch(...) {
      tracker->setException(std::current_exception());
    }
    std::unique_lock<std::mutex> lck(tracker->mutex());
    tracker->releaseWriteLock(blocknum);
  }
};

std::shared_ptr<ISDGenericDataset>
ISDGenericDataset::injectAsyncInstance(
     const std::shared_ptr<ISDGenericDataset>& relay,
     std::int64_t highwater,
     const LoggerFn& logger)
{
  return std::make_shared<SDGenericDatasetAsyncWrite>(relay, highwater, logger);
}

} // namespace

#endif // HAVE_SD
