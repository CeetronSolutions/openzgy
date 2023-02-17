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

#include "../declspec.h"

#include <functional>
#include <string>

#ifdef HAVE_SD // Most of the file

#include "logger.h"
#include "file_consolidate.h"
#include "file_performance.h"
#include "file_parallelizer.h"
#include "fancy_timers.h"
#include "mtguard.h"
#include "../exception.h"
#include "../iocontext.h"
#include "environment.h"
#include "sdinterface.h"

#include <vector>
#include <memory>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <limits>
#include <string.h>
#include <algorithm>
#include <numeric>
#include <mutex>
#include <omp.h>

// It would have been nice to have similar include paths in Linux and Windows
// but that is a cosmetic issue only and it is only a problem in this file.
#ifndef _WIN32
#include <SDManager.h>
#include <SDDatasetDisposition.h>
#include <SDUtils.h>
#include <Constants.h>
#include <HttpContext.h>
#include <SDException.h>
#else
#include <SDAPI/SDManager.h>
#include <SDAPI/SDDatasetDisposition.h>
#include <SDAPI/SDUtils.h>
#include <SDAPI/Constants.h>
#include <SDAPI/HttpContext.h>
#include <SDAPI/SDException.h>
#endif

// SDAPI footprint, might not be accurate.
// Use e.g. to look for places where a try/catch might be needed.
//
// SDGenericDataset constructors, destructor
// \b(open|close|readBlock|writeBlock|deleteBlock|getBlockNum|getBlocksSize|getSize|getConsistencyID|getSerializedContext)\(
// Not used yet in C++, the first two supported in Python+SdGlue.
// \b((get|set)(MetaData|SeismicMeta|Tags|ReadOnlyMode))\(
//
// SDManager constructors, destructor
// \b(setAuthProviderFrom(String|File|ImpToken)|setLogStatus)\(

/** \cond SSTORE */

/**
 * \file: file_sd.cpp
 * \brief Low level I/O, Seismic Store.
 */

// IByteIO::Disposition -> OpenMode
//OpenMode:    { Closed = 0, ReadOnly, ReadWrite, Truncate, Delete, };
//Disposition: { Create, Read, Update, Delete, Access, List };

// Plans for utilizing code from byteio_sd.cpp in ZGY-Cloud.
//
// cherry pick code from here for read and write.
//     ByteIoSeismicDrive::closeWrite()
//
// bells and whistles for configuration.
//     slurp()
//     expandArgument()
//
// Not now. Don't keep a long lived cache.
//     class SDGenericDatasetCache
//
// Worry: CTag / ETag handling? Not now
//     ByteIoSeismicDrive::readETag()
//
// Might need soon.
//     checkAccess()
//
// Utilities, might need later.
// Need a way for API users to access these directly
//     getTag()
//     setTag()
//     getMetaData()
//     setMetaData()
//     getSeisMeta()
//     setSeisMeta()
//
// maybe later, for tokentype and for error messages.
// TODO-Low: if imp token is refreshed then keep those in a cache.
// Might need support from SDAPI for that.
//     parseJsonString()
//     decodeJwt()
//     getTokenType()

using OpenZGY::IOContext;
using OpenZGY::SeismicStoreIOContext;
namespace InternalZGY {
#if 0
}
#endif

class SDGenericDatasetWrapper;
class DatasetInformation;

/**
 * Access data in seismic store as a linear file even when the dataset
 * has multiple segments. There are some limitations on write.
 *
 *   - Writes starting at EOF are allowed, and will cause a new segment
 *     to be written.
 *
 *   - Writes starting past EOF, signifying a hole in the data, are not
 *     allowed.
 *
 *   - Writes starting before EOF are only allowed if offset,size exactly
 *     matches a previous write. This will cause that segment to be rewritten.
 *
 *   - Possible future extension: For the  last segment only offset
 *     needs to match. This means the last segment may be resized.
 *
 * For read the class provides a readv() method to do scatter/gather reads.
 * The code will then consolidate adjacent bricks to get larger brick size
 * sent to SDAPI. Optionally parallelize requests that cannot be consolidated.
 *
 * Thread safety when used for reading:
 * Designed to be thread safe as no internal data structures should change
 * after the file has been opened. Any lazy-evaluated information needs
 * to do appropriate locking.
 *
 * Thread safety when used for writing:
 * Not thread safe. See SeismicStoreFile::xx_write.
 *
 * Thread safety when closing a file:
 * Not thread safe.
 *
 * The class is noncopyable with copy and assign method deleted.
 * Not that the users could easily copy it anyway, as the class should
 * always be used via the IFileADT interface.
 */
class SeismicStoreFile : public FileADT
{
  SeismicStoreFile(const SeismicStoreFile&) = delete;
  SeismicStoreFile& operator=(const SeismicStoreFile&) = delete;
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;
public:
  SeismicStoreFile(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  virtual ~SeismicStoreFile();
  static std::shared_ptr<IFileADT> xx_make_instance(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  // Functions from IFileADT
  virtual void xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_readv(const ReadList& requests, bool parallel_ok=false, bool immutable_ok=false, bool transient_ok=false, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_close();
  virtual std::int64_t xx_eof() const;
  virtual std::vector<std::int64_t> xx_segments(bool complete) const override;
  virtual bool xx_iscloud() const override;
  // Functions from IFileBase
  virtual void deleteFile(const std::string& filename, bool missing_ok) const;
  virtual std::string altUrl(const std::string& filename) const;
  virtual std::string idToken() const;
private:
  void do_write_one(const void*  const data, const std::int64_t blocknum, const std::int64_t size, const bool overwrite);
  void do_write_many(const void*  const data, const std::int64_t blocknum, const std::int64_t size, const std::int64_t blobsize, const bool overwrite);
public:
  // Needed by SeismicStoreFileDelayedWrite.
  OpenMode _mode() const;
  // The raw SDGenericDataset is needed by SeismicStoreFileDelayedWrite
  // when opening a file for update.
  std::shared_ptr<SDGenericDatasetWrapper> datasetwrapper() const {return _dataset;}
  bool _mylogger(int priority, const std::string& message) const;
  bool _sslogger(int priority, const std::ios& ss) const;
  void _set_backoff(ISDGenericDataset* sdgd);
  std::shared_ptr<ISDGenericDataset>_open_dataset_ro(
       const std::shared_ptr<seismicdrive::SDManager>& manager,
       const std::string& filename,
       const std::unordered_map<std::string, std::string>& extra,
       bool sd_ds_log,
       const SeismicStoreIOContext *context);
  std::shared_ptr<ISDGenericDataset>_open_dataset_rw(
       const std::shared_ptr<seismicdrive::SDManager>& manager,
       const std::string& filename,
       bool truncate,
       const std::unordered_map<std::string, std::string>& extra,
       bool sd_ds_log,
       const SeismicStoreIOContext *context);
private:
  /**
   * This class is used by _split_by_segment to describe a request as seen by
   * seismic store.
   *
   * Thread safety:
   * Modification may lead to a data race. This should not be an issue,
   * because instances are only meant to be modified when created or
   * copied or assigned prior to being made available to others.
   */
  class RawRequest
  {
  public:
    std::int64_t blocknum;	// Seismic store segment a.k.a. block
    std::int64_t local_offset;	// File offset inside this blocknum
    std::int64_t local_size;	// How much to read from this block
    std::int64_t outpos;
    RawRequest(std::int64_t a_blocknum,
               std::int64_t a_offset,
               std::int64_t a_size,
               std::int64_t a_outpos)
      : blocknum(a_blocknum)
      , local_offset(a_offset)
      , local_size(a_size)
      , outpos(a_outpos)
    {
    }
  };
  typedef std::vector<RawRequest> RawList;
  RawList _split_by_segment(const ReadList& requests);
  void _cached_read(/*TODO-Low: seg, offset, view*/);
private:
  // TODO-Low: To improve isolation, the user visible context should
  // be copied into an equivalent InternalZGY::SeismicStoreConfig.
  // The downside is that it gets more tedious to maintain.
  std::shared_ptr<OpenZGY::SeismicStoreIOContext> _config;
  std::shared_ptr<SDGenericDatasetWrapper> _dataset;
  LoggerFn _logger;
  // As long as we don't inherit FileCommon we need our own timers.
  std::shared_ptr<SummaryPrintingTimerEx> _rtimer; // Access is thread safe
  std::shared_ptr<SummaryPrintingTimerEx> _wtimer; // Access is thread safe
};

/**
 * Improve on SeismicStoreFile, have it buffer large chunks of data before
 * writing it out to a new segment.
 *
 *   - Writes starting at EOF are allowed, and will buffer data in the
 *     "open segment" until explicitly flushed.
 *
 *   - Writes starting past EOF, signifying a hole in the data, are not
 *     allowed.
 *
 *   - Writes fully inside the open segment are allowed.
 *
 *   - Writes starting before the open segment are only allowed if
 *     offset,size exactly matches a previous write. This will cause that
 *     segment to be rewritten. As a corollary, writes canot span the
 *     closed segment / open segment boundary.
 *
 *   - Possible future extension: For the  last segment only offset
 *     needs to match. This means the last segment may be resized.
 *     Why we might want this: On opening a file with existing
 *     data bricks we might choose to read the last segment and
 *     turn it into an open segment. Then delete (in memory only)
 *     the last segment. When it is time to flush the data it gets
 *     rewritten. This allows adding bricks to a file, while still
 *     ensuring that all segments except first and last need to be
 *     the same size. Note that there are other tasks such as
 *     incrementally updating statistics and histogram that might
 *     turn out to be a lot of work.
 *
 *   - When used to create ZGY files, caller must honor the convention
 *     that all segments except the first and last must have the same size.
 *
 *   - Caveat: The fact that random writes are sometimes allowed, sometimes
 *     not depending on the segment number violates the principle of
 *     least surprise. And makes for more elaborate testing. For ZGY
 *     it is quite useful though. ZGY can recover from a ZgySegmentIsClosed
 *     exception by abandoning (leaking) the current block and write it
 *     to a new location. With a typical access pattern this will happen
 *     only occasionally.
 *
 * Thread safety:
 * Not thread safe by design, as it is only used for files opened for write.
 * If a file is opened for read/write (currently not supported)
 * or when running finalize it is safe to read data as long as no writes
 * can be pending.
 */
class SeismicStoreFileDelayedWrite : public FileADT
{
  std::shared_ptr<OpenZGY::SeismicStoreIOContext> _config;
  std::shared_ptr<SeismicStoreFile> _relay;
  std::vector<char> _open_segment;
  UsageHint _usage_hint;
  std::shared_ptr<SummaryPrintingTimerEx> _ctimer; // Access is thread safe

  SeismicStoreFileDelayedWrite(const SeismicStoreFileDelayedWrite&) = delete;
  SeismicStoreFileDelayedWrite& operator=(const SeismicStoreFileDelayedWrite&) = delete;

public:
  SeismicStoreFileDelayedWrite(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  virtual ~SeismicStoreFileDelayedWrite();
  static std::shared_ptr<IFileADT> xx_make_instance(const std::string& filename, OpenMode mode, const OpenZGY::IOContext *iocontext);
  // Methods from IFileBase
  virtual void deleteFile(const std::string& name, bool missing_ok) const override
  {
    _relay->deleteFile(name, missing_ok);
  }
  virtual std::string  altUrl(const std::string& name) const override
  {
    return _relay->altUrl(name);
  }
  virtual std::string  idToken() const override
  {
    return _relay->idToken();
  }
  // Methods from IFileADT:
  virtual void xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_readv(const ReadList& requests, bool parallel_ok=false, bool immutable_ok=false, bool transient_ok=false, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint=UsageHint::Unknown) override;
  virtual void xx_close() override;
  virtual std::int64_t xx_eof() const override;
  virtual std::vector<std::int64_t> xx_segments(bool complete) const override;
  virtual bool xx_iscloud() const override;

private:
  void _reopen_last_segment();
  void _flush_part(std::int64_t this_segsize);
  void _flush(bool final_call);
};

/////////////////////////////////////////////////////////////////////////////
//    class DatasetInformation   ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * \brief Cached information for a SDGenericDataset.
 *
 * Copied nearly verbatim from byteio_sd.cpp in the old accessor.
 *
 * All fields will be filled in when any of them are first needed.
 * An exception is thrown if some of the information is inaccessible.
 * For new datasets, the information may either be set explicitly by
 * the writer or it may be set on the first write.
 * For existing datasets the sdapi is queried.
 * Even if the dataset is opened read/write the information is not
 * expected to change. Some other mechanism needs to deal with
 * stale data in that case.
 *
 * Thread safety:
 * All access SDGenericDatasetWrapper::info_ will be protected by
 * SDGenericDatasetWrapper::mutex_. Methods that expect to change data
 * (currently only updateOnWrite()) will need some special handling.
 * TODO-High: there is still a snag here.
 */
class DatasetInformation
{
private:
  // Number of blocks on file.
  // On write, this will include the holes if data is written out of order.
  std::int64_t block_count_;

  // Size of first block. On write this might need to be explicitly set if the
  // writer decides to output blocks non sequentially. It is not possible to
  // just assume block 0 is the same size as block 1 because block 1 might
  // be the last block on file and this might have a different size.
  std::int64_t block0_size_;

  // Size of blocks that are neither the first nor last on the file.
  // Zero if there is only a single block on the file.
  std::int64_t block1_size_;

  // Size of the last block, which is allowed to be smaller.
  // Will be zero if there is just one block.
  std::int64_t last_block_size_;

  // For debugging
  SeismicStoreFile::LoggerFn logger_;

public:
  explicit DatasetInformation(const SeismicStoreFile::LoggerFn& logger);
  explicit DatasetInformation(ISDGenericDataset* sdgd, const SeismicStoreFile::LoggerFn& logger);
  std::string toString() const;

public:
  std::int64_t blockCount() const { return block_count_; }
  std::int64_t block0Size() const { return block0_size_; }
  std::int64_t block1Size() const { return block1_size_; }
  std::int64_t lastBlockSize() const { return last_block_size_; }

public:
  std::int64_t totalSize() const;
  std::vector<std::int64_t> allSizes(bool complete) const;
  void getLocalOffset(std::int64_t offset, std::int64_t size, std::int64_t *blocknum, std::int64_t *local_offset, std::int64_t *local_size) const;
  void checkOnWrite(std::int64_t blocknum, std::int64_t blocksize) const;
  void updateOnWrite(std::int64_t blocknum, std::int64_t blocksize);
};

DatasetInformation::DatasetInformation(const SeismicStoreFile::LoggerFn& logger)
  : block_count_(0)
  , block0_size_(0)
  , block1_size_(0)
  , last_block_size_(0)
  , logger_(logger ? logger : LoggerBase::emptyCallback())
{
}

/**
 * Create and poplulate an instance.
 * Catching and translating exceptions from SDAPI is done by caller.
 */
DatasetInformation::DatasetInformation(ISDGenericDataset* sdgd, const SeismicStoreFile::LoggerFn& logger)
  : block_count_(0)
  , block0_size_(0)
  , block1_size_(0)
  , last_block_size_(0)
  , logger_(logger ? logger : LoggerBase::emptyCallback())
{
  // Note that sdapi is a bit confusing with respect to signed/unsigned.
  // getBlockNum() returns an unsigned (uint64_t).
  // getSize() and getBlockSize() return long long, may return -1 as error.
  // I will treat all of them as signed, since there is no way that
  // there can be more than 2^63 blocks.
  long long nblocks = sdgd->getBlockNum();
  long long nbytes = sdgd->getSize();
  std::vector<std::string> names;
  if (nblocks >= 1)
    names.push_back("0");
  if (nblocks >= 2)
    names.push_back("1");
  if (nblocks >= 3)
    names.push_back(std::to_string(nblocks - 1));
  std::vector<long long> sizearray;
  if (nblocks <= 0) {
  }
#if 0 // Chicken...
  else if (nblocks == 1) {
    // Not using getSize() because I do not trust it.
    // SDAPI requires a hint check_and_overwrite=true when writing
    // a block that might exist already. If the hint is missing
    // then getSize() will return the wrong result. OpenZGY has
    // no control over who wrote the file. Defensive programming
    // says to just don't use the getSize() until the current
    // behavior (which I consider a bug) is fixed. Suggestion:
    // would it be possible to scan the file on exit and fix up
    // any incorrect size?
    sizearray.push_back(nbytes);
  }
#endif
  else {
    sizearray = sdgd->getBlockSizes(names);
  }
  if (nblocks < 0)
    throw OpenZGY::Errors::ZgyInternalError("Unable to get block count for SDGenericDataset");
  for (const auto size : sizearray)
    if (size <= 0)
      throw OpenZGY::Errors::ZgyInternalError("Unable to get segment size for SDGenericDataset");

  this->block_count_ = nblocks;
  this->block0_size_ = nblocks >= 1 ? sizearray.at(0) : 0;
  this->block1_size_ = nblocks >= 2 ? sizearray.at(1) : 0;
  this->last_block_size_ = nblocks >= 3 ? sizearray.at(2) : this->block1_size_;

  // Not throwing exceptions here, because sdgd->getSize() is not reliable.
  // But the problem *could* be that the block sizes of all the middle blocks
  // is not the same. Which would be a real problem.
  if (nbytes >= 0 && nbytes != (long long)totalSize())
    this->logger_(0, "Dataset has inconsistent size");

  if (this->logger_(1, ""))
    this->logger_(1, toString());
}

std::string
DatasetInformation::toString() const
{
  std::stringstream ss;
  ss << "Total " << totalSize()
     << " in " << blockCount() <<" segments"
     << " of size "<< block0Size()
     << ", " << block1Size()
     << ", " << lastBlockSize();
  return ss.str();
}

/**
 * Return the total file size, including any holes caused by bricks
 * not written yet. The result currently computes the answer based on
 * the individual block sizes assuming all blocks except the first and
 * last will have the same size. This is more reliable than asking
 * sdapi for the total size. And if the block size assumption is false
 * then we are completely hosed anyway. It is possible to verify the
 * assumption but it is probably too expensive to do so.
 *
 * Thread safety: Safe because instances are immutable once constructed
 * and made available.
 */
std::int64_t
DatasetInformation::totalSize() const
{
  switch (block_count_) {
  case 0: return 0;
  case 1: return block0_size_;
  default: return (block0_size_ +
                   (block_count_ - 2) * block1_size_ +
                   last_block_size_);
  }
}

/**
 * Return the total file size broken down into segments, not including
 * the "open" segment which DatasetInformation doesn't know about.
 */
std::vector<std::int64_t>
DatasetInformation::allSizes(bool complete) const
{
  switch (block_count_) {
  case 0: return std::vector<std::int64_t>{};
  case 1: return std::vector<std::int64_t>{block0_size_};
  case 2: return std::vector<std::int64_t>{block0_size_, last_block_size_};
  default: {
    std::vector<std::int64_t> result;
    result.push_back(block0_size_);
    result.push_back(block1_size_);
    if (complete)
      for (int ii = 0; ii < block_count_ - 3; ++ii)
        result.push_back(block1_size_);
    result.push_back(last_block_size_);
    return result;
  }
  }
}

/**
 * Do consistency checks before data is written.
 * Some of these might be redundant due to checks in the caller.
 * E.g. caller raises SegmentIsClosed if attempting to write "backwards".
 * Throws exceptions on error.
 */
void
DatasetInformation::checkOnWrite(std::int64_t blocknum, std::int64_t blocksize) const
{
  if (blocknum < 0 || blocknum > std::numeric_limits<int>::max()) {
    throw OpenZGY::Errors::ZgyInternalError("Cannot write block " + std::to_string(blocknum) + " because it is out of range.");
  }
  if (blocknum != 0 && block0_size_ == 0) {
    throw OpenZGY::Errors::ZgyInternalError("Cannot write block " + std::to_string(blocknum) + " before size of block 0 is known.");
  }
  if (blocksize < 1) {
    throw OpenZGY::Errors::ZgyInternalError("Cannot write less that 1 byte.");
  }
  if (blocknum == 0) {
    // Write or overwrite block 0.
    if (block0_size_ != 0 && block0_size_ != blocksize)
      throw OpenZGY::Errors::ZgyInternalError("Cannot change the size of block zero");
  }
  else if (blocknum + 1 < block_count_) {
    // Overwrite a block in the middle.
    if (block1_size_ != blocksize)
      throw OpenZGY::Errors::ZgyInternalError("Blocks must have the same size except the first and last");
  }
  else if (blocknum + 1 == block_count_) {
    // Overwrite the last block, which is not block 0.
    // TODO-Low: Technically I might have allowed this.
    // If update is to be supported then I probably need to.
    if (blocksize != last_block_size_)
      throw OpenZGY::Errors::ZgyInternalError("Cannot change the size when re-writing the last block");
  }
  else if (blocknum == block_count_) {
    // Append a block, which is not block 0.
    // This is the new "last block" which means what we previous thought
    // was the last block is in fact a middle one. Which means the size
    // of the former last block must be the same as block 1.
    // (Note that the former last block might actually *be* block 1).
    if (block1_size_ != 0 && blocksize > block1_size_)
      throw OpenZGY::Errors::ZgyInternalError("block " + std::to_string(blocknum) + " is too large");
    if (block_count_ != 0 && last_block_size_ != block1_size_)
      throw OpenZGY::Errors::ZgyInternalError("block " + std::to_string(block_count_ -1) + " had wrong size");
  }
  else {
    // Trying to write sparse data.
    throw OpenZGY::Errors::ZgyInternalError("block " + std::to_string(blocknum) + " written out of sequence");
  }
}

/**
 * Update cached size information after data is successfully written.
 * checkOnWrite() must have been called already.
 *
 * Thread safety: NOT thread safe.
 * Do not invoke SDGenericDatasetWrapper::info()->updateOnWrite() directly.
 * Call the thread safe SDGenericDatasetWrapper::updateOnWrite() instead.
 * That one wll make sure the smart pointer being updated is unique.
 */
void
DatasetInformation::updateOnWrite(std::int64_t blocknum, std::int64_t blocksize)
{
  if (blocknum == 0) {
    // Write or overwrite block 0.
    block0_size_ = blocksize;
    block_count_ = std::max(block_count_, blocknum + 1);
  }
  else if (blocknum + 1 < block_count_) {
    // Overwrite a block in the middle.
    // Size cannot change, so do nothing.
  }
  else if (blocknum + 1 == block_count_) {
    // Overwrite the last block, which is not block 0.
    // Redundant if checkOnWrite() forbids changing size.
    last_block_size_ = blocksize;
  }
  else if (blocknum == block_count_) {
    // Append a block which is not block 0.
    if (block1_size_ == 0)
      block1_size_ = blocksize;
    last_block_size_ = blocksize;
    block_count_ = std::max(block_count_, blocknum + 1);
  }
}

/**
 * Given a linear offset and count, convert this to a block local address.
 * The returned count may be smaller than requested in case the request
 * crossed a segment boundary. In that case the caller will need to read
 * in a loop.
 *
 * "blocks" refer to the block number in Seismic Store, not the potentially
 * larger logical blocks used by SeismicStoreFileDelayedWrite.
 *
 * Postcondition: If blocknum is returned as the last block, local_size
 * will be returned as requested size. If this were not so, the calling
 * function would be likely to loop forever.
 */
void
DatasetInformation::getLocalOffset(std::int64_t offset, std::int64_t size, std::int64_t *blocknum, std::int64_t *local_offset, std::int64_t *local_size) const
{
  if (offset < 0 || size < 0)
    throw OpenZGY::Errors::ZgyInternalError("Offset and size cannot be negative.");
  else if (size > std::numeric_limits<std::int64_t>::max() - offset)
    throw OpenZGY::Errors::ZgyInternalError("Overflow in offset + size.");
  else if (offset + size > totalSize()) {
    if (this->logger_(1, "")) {
      std::stringstream ss;
      ss << "Reading past EOF: read("
         << "off=" << offset
         << ", size=" << size
         << ", end=" << offset+size
         << ") dataset: " << toString()
         << std::endl;
      this->logger_(1, ss.str());
    }
    throw OpenZGY::Errors::ZgyInternalError("Reading past EOF");
  }

  if (block_count_ <= 1) {
    // In first block which is also last block.
    *blocknum = 0;
    *local_offset = offset;
    *local_size = size;
  }
  else if (block0_size_ == 0 || block1_size_ == 0) {
    // Can only happen when writing. Files open for read will have the sizes
    // set up when the DatasetInformation is first accessed.
    throw OpenZGY::Errors::ZgyInternalError("getLocalOffset() called before size is known.");
  }
  else if (offset < block0_size_) {
    // In first block, and there are more blocks following.
    *blocknum = 0;
    *local_offset = offset;
    *local_size = std::min(size, block0_size_ - offset);
  }
  else {
    const std::int64_t bnum = std::min(block_count_ - 1, ((offset - block0_size_) / block1_size_) + 1);
    const std::int64_t segment_start = block0_size_ + (bnum - 1) * block1_size_;
    *blocknum = bnum;
    *local_offset = offset - segment_start;
    if (bnum + 1 < block_count_)
      *local_size = std::min(size, block1_size_ - (offset - segment_start));
    else
      *local_size = size; // In last block.
  }
}

/////////////////////////////////////////////////////////////////////////////
//    class SDGenericDatasetWrapper   ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Wrapper around seismicdrive::SDGenericDataset.
 * Instances are typically managed by smart pointers.
 * Copied nearly verbatim from byteio_sd.cpp in the old accessor.
 *
 * What this adds compared to a raw SDGenericDataset:
 *
 * - We keep a smart pointer reference to the SDManager that was used
 *   to create the SDGenericDataset. The sdapi documentation does not
 *   state this explicitly, but I assume it is unhealthy to destruct
 *   the manager before the dataset is closed. And using an application
 *   global manager is not an option due to privilege separation.
 *
 * - The instance remembers whether it was created for create, read, etc.
 *
 * - We may later add wrappers for SDGenericDataset members for caching
 *   purposes (block size, file size, etc.) in case the sdapi doesn't.
 *
 * - We may later add wrappers for the other members to tweak the
 *   interface, map error messages, add logging, etc. Or to add
 *   additional checking such as not reading past eof / end of block.
 *
 * Thread safety: The class itself is thread safe. The data being wrapped
 * might not be.
 *
 * All mutable data is protected by this->mutex_. Access methods for
 * those return a smart pointer. If info() is called on one thread
 * while another thread is doing a write operation (which is actually
 * not allowed) then it is unspecified whether the returned value is
 * before or after the write. It is also unspecified whether
 * updateOnWrite() will have any effect in that case. The information
 * is still consistent though, and may be used for logging etc.
 *
 * CAVEAT: Make sure the returned pointer remains in scope long
 * enough. E.g. this is NOT safe, and might crash every 7 years or so.
 * Except on customer sites where it may crash every 7 minutes.
 *
 *    seismicdrive::SDGenericDataset& dataset = *wrapper->dataset(); // NO!!!
 *    foo(dataset); // the returned smart pointer is already deleted.
 *
 * TODO-Low: Yagni: virgin_ is not used. It is related to the CTag mechanism.
 * It is still being discussed whether that needs to be ported from the
 * old accessor. It might not be needed if we go for immutable ZGY files.
 */
class SDGenericDatasetWrapper
{
  std::shared_ptr<seismicdrive::SDManager> manager_;
  std::shared_ptr<ISDGenericDataset> dataset_;
  std::shared_ptr<const DatasetInformation> info_;
  OpenMode disposition_;
  bool virgin_; // If true, the cached CTag should be ok.
  SeismicStoreFile::LoggerFn logger_;
  mutable std::mutex mutex_; // Protect all members.
  OpenZGY::SeismicStoreIOContext::tokencb_t tokenrefresh2_;
  std::string tokenmessage_;

public:
  typedef std::shared_ptr<SDGenericDatasetWrapper> Ptr;
  SDGenericDatasetWrapper(std::shared_ptr<seismicdrive::SDManager> manager,
                          std::shared_ptr<ISDGenericDataset> dataset,
                          OpenMode disp,
                          const SeismicStoreFile::LoggerFn& logger)
    : manager_(manager), dataset_(dataset), disposition_(disp), virgin_(true)
    , logger_(logger ? logger : LoggerBase::emptyCallback())
    , mutex_()
    , tokenrefresh2_()
    , tokenmessage_("Token not initialzed")
  {
  }
  ~SDGenericDatasetWrapper();
  std::shared_ptr<ISDGenericDataset> dataset() {
    std::lock_guard<std::mutex> lk(mutex_);
    return dataset_;
  }
  std::shared_ptr<seismicdrive::SDManager> manager() {
    std::lock_guard<std::mutex> lk(mutex_);
    return manager_;
  }
  OpenMode disposition() const {
    // This is immutable, except in close() which is not threadsafe
    // anyway, so no lock is needed.
    return disposition_;
  }
  std::shared_ptr<const DatasetInformation> info() {
    std::lock_guard<std::mutex> lk(mutex_);
    if (!info_) {
      try {
      switch (disposition_) {
      case OpenMode::Truncate:
        info_.reset(new DatasetInformation(logger_));
        break;
      case OpenMode::ReadOnly:
      case OpenMode::ReadWrite:
        info_.reset(new DatasetInformation(dataset_.get(), logger_));
        break;
      case OpenMode::Closed:
      default:
        throw OpenZGY::Errors::ZgyInternalError("DatasetInformation: Dataset not open.");
      }
      }
      catch (const std::exception& ex) {
        throwCloudException(ex, "Initialize");
      }
    }
    return info_;
  }

  /**
   * Close the underlying SDGenericDataset.
   * If an exception occurs then the wrapper will still see
   * the dataset as closed. And dataset() will return empty.
   *
   * Caveat: There was at one point a bug where SDGenericDataset::close()
   * threw an appropriate exception (e.g. credentials timed out) then
   * destructing the dataset could also throw. Exceptions from a
   * destructor usually causes the program to terminate. Especially in
   * C++11 and later. There are workarounds with multiple try/catch that
   * *might* help but it should be a lot easier to fix the underlying bug.
   */
  void wrapper_close(bool set_readonly)
  {
    auto victim = dataset_;
    OpenMode disp = disposition_;
    dataset_.reset();
    disposition_ = OpenMode::Closed;
    if (victim) {
      if (set_readonly &&
          (disp == OpenMode::Truncate || disp == OpenMode::ReadWrite))
        victim->setReadonlyMode(true);
      victim->close();
      victim.reset(); // Any throw from SDGenericDataset dtor happens here.
    }
  }

  void updateDataset(std::shared_ptr<ISDGenericDataset> dataset) {
    std::lock_guard<std::mutex> lk(mutex_);
    dataset_ = dataset;
    info_.reset();
    // This would give a slight performance boost, saving a single call
    // to checkCTag() after a file has changed. This happens so rarely
    // that it isn't worth the extra testing.
    //virgin_ = true;
  }

  /**
   * Adjust the information after data has been written to the cloud.
   *
   * Thread safety:
   * Multiple concurrent writers will have a race condition, but that
   * scenario is expressly forbidden anyway. Smart pointers from
   * info() held by callig code are race free because they are
   * immutable. updateOnWrite() makes a new pointer.
   */
  void updateOnWrite(std::int64_t blocknum, std::int64_t blocksize)
  {
    std::shared_ptr<DatasetInformation> updated;
    updated.reset(new DatasetInformation(*this->info()));
    updated->updateOnWrite(blocknum, blocksize);
    info_ = updated;
  }

  bool isTouched() const {
    std::lock_guard<std::mutex> lk(mutex_);
    return !virgin_;
  }

  bool touch() {
    // Called to inform us that getCTag() might return stale data.
    //
    // if this->virgin_, we have just created the wrapper which means
    // we have just opened it. Only "List" entries don't get opened
    // immediately. But those aren't cached so they are not relevant
    // here. Since we just created it, we can call getCTag() without
    // calling checkCTag() first. Saving one database hit.
    //
    // The flag must be reset at least when an existing dataset accessor
    // is being re-used to open a file again. Currently it will be reset
    // much more often (every time the accessor is used to read data)
    // and it won't be set true again by updateDataset(). This will
    // simplify testing and has little impact on performance.
    std::lock_guard<std::mutex> lk(mutex_);
    bool old = virgin_;
    virgin_ = false;
    return old != virgin_;
  }

  /**
   * Called from the catch block after calling some method in SDAPI.
   * If no token was provided, this is an "I told you so" error.
   * It could have been reported earlier but it is possible that
   * the cloud library somehow managed to handle authenticcation
   * itself. The actual exception reported by SDAPI is in that case
   * less interesting.
   *
   * A missing token is probably a user error (with "user" in this
   * case being the application) and is reported as such.
   *
   * If the exception has alredy been wrapped by a ZgyException then
   * it will be re-thrown without modification.
   */
  void throwCloudException(const std::exception& ex, const char *message) {
    if (this->logger_(1, "")) {
      std::stringstream ss;
      ss << "Oops (" << message << ")" << " ("
         << tokenmessage_ << "): " << ex.what();
      this->logger_(1, ss.str());
    }
    if (dynamic_cast<const OpenZGY::Errors::ZgyError*>(&ex))
      throw;
    else if (!tokenmessage_.empty())
      throw OpenZGY::Errors::ZgyUserError(tokenmessage_);
    else
      throw OpenZGY::Errors::ZgyInternalError
        (std::string(message) + ": Seismic Store: " + std::string(ex.what()));
    // TODO-Low: Should I just re-throw the sdapi error instead?
    // TODO-Low: Either way, be more consistent about the answer.
  }

  /**
   * Set credentials from a single string, supporting any
   * authorization mode except "Callback". Calling this method gives
   * less compile time checks than invoking one of the specific
   * setAuthProviderXXX() methods. On the positive side the caller
   * need not worry about which authorization schemes were built into
   * SDAPI. Also the application need only pass around a single
   * authorization string.
   *
   * The input string is "<type>:arg1:arg2:..." with the number of
   * arguments depending on <type>. If no delimiter is found a type
   * of "String" is assumed. The function assumes that none of the
   * arguments themselves contain the delimiter.
   *
   * The list of allowed modes depends on how SDAPI was built.
   * It would have been really nice if SDManager could expose
   * this function itself. And even better, use a plug-in mechanism
   * so the check can be done at runtime.
   *
   * A "Client Credentials Grant" means there is no user involved.
   * Any process that knows the secrets can create access tokens
   * indefinitely (or until the secret is revoked).
   *
   * Types that require HAVE_INTERNAL_SD use SAuth credentials
   * and will interact with $SAUTH_TKSVC_URL, not the sdurl. E.g.
   * https://p4d.csi.cloud.slb-ds.com/v2/token.
   *
   * Caveat: For completeness this needs to be duplicated in
   * native/sdglue/sdgluemodule.cpp and needs to be called from the
   * low level test_sdutils. But the former case is for the pure
   * python implementation that is no longer supported and the
   * latter can be achieved by a fairly small kludge.
   */
  static void setAuthProvider(seismicdrive::SDManager *mgr,
                              const std::string& argstring,
                              const SeismicStoreFile::LoggerFn& logger)
  {
    std::istringstream is(argstring);
    std::string part;
    std::vector<std::string> args;
    while (std::getline(is, part, ':'))
      args.push_back(part);
    if (args.size() <= 1)
      args = std::vector<std::string>{"String", argstring};

    const std::string cred = "Credentials \"" + args[0] + "\" ";
    if (logger && logger(3, "")) {
      std::stringstream redacted;
      for (const std::string& s : args) {
        if (s.size() < 16)
          redacted << ":" << s;
        else
          redacted << ":" << s.substr(0, 5) << "..." << s.substr(s.size()-5);
      }
      logger(3, "Credentials: " + redacted.str().substr(1));
    }

    if (args[0] == "String") {
      if (args.size() != 2)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 1 argument");
      mgr->setAuthProviderFromString(args[1]);
    }
    else if (args[0] == "ImpToken" || args[0] == "imptoken") {
      if (args.size() != 2)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 1 argument");
      // The only difference between "ImpToken" and "String"
      // is that SDAPI will try to refresh the former if it
      // has timed out. TODO: Methinks SDAPI ought to have
      // figured out itself whether the token allows refresh.
      mgr->setAuthProviderFromImpToken(args[1]);
    }
 #ifdef HAVE_INTERNAL_SD
    else if (args[0] == "FileV2" || args[0] == "FILE") {
      if (args.size() != 4)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 3 arguments");
      // Token and refresh token stored in a file.
      // Automatically switches to "SAuthV2" if on a VM.
      // Args: clientID, clientSecret, filepath
      mgr->setAuthProviderFromFileV2(args[1], args[2], args[3]);
    }
    else if (args[0] == "SauthV2") {
      if (args.size() != 4)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 3 arguments");
      // "Client Credentials Grant" associated with VM.
      // Automatically switches to "FileV2" if not on a VM.
      // This is why client id and secret are needed.
      // Currently it makes no difference if the application
      // uses "FileV2" or "SauthV2". But that might change.
      // Args: clientID, clientSecret, clientCredentialsFile
      mgr->setAuthProviderSauthV2(args[1], args[2], args[3]);
    }
    else if (args[0] == "ServiceSauthV2") {
      if (args.size() != 3)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 2 arguments");
      // "Client Credentials Grant" unlocked by just the id and secret.
      // Args: clientID, clientSecret
      mgr->setAuthProviderServiceSauthV2(args[1], args[2]);
    }
    else if (args[0] == "SauthImpersonationToken") {
      if (args.size() != 5)
        throw OpenZGY::Errors::ZgyUserError(cred + "expects 4 arguments");
      // Impersonation token signed by SAuth.
      // Args: clientId, clientSecret, impToken, impTokenContext
      mgr->setAuthProviderSauthImpersonationToken(args[1], args[2], args[3], args[4]);
    }
#endif
    else {
      throw OpenZGY::Errors::ZgyUserError(cred + "not supported");
    }
  }

  static std::string tokenCallbackForwarder(const void *app_data)
  {
    auto fn = reinterpret_cast<const std::function<std::string()>*>(app_data);
    return (*fn)();
  }

  static std::string tokenCallbackInvalid(const void *)
  {
    throw OpenZGY::Errors::ZgyUserError
      ("Refreshing the token at this point is no longer possible");
  }

  /**
   * Pass initial credentials down to the SDAPI layer.
   *
   * Thread safety: This is meant to be called from the constructor
   * when opening a file. So there should not be any race condition.
   * CAVEAT: If the manager is cached and shared between open files
   * then this raises some challenges.
   */
  void authorizeManager(
       const std::string& token,
       const OpenZGY::SeismicStoreIOContext::tokencb_t& tokencb2)
  {
    std::lock_guard<std::mutex> lk(mutex_);

    if (tokencb2) {
      // Careful: the manager indirectly gets an unsafe pointer to the
      // functor. Need to control the lifetime of the functor itself,
      // as well as requiring its implementation to remain valid.
      tokenrefresh2_ = tokencb2;
      manager_->setAuthProviderCallback
        (&SDGenericDatasetWrapper::tokenCallbackForwarder, &tokenrefresh2_);
      // If the callback returned an empty, expired, or invalid token
      // then show the actual SDAPI error text. tokenmessage_ should
      // only be used if the code here is reasonably sure it knows
      // what went wrong, and wants to throw a more readable error.
      tokenmessage_ = "";
      return;
    }

    // A missing token at this point is probably an error, but in case
    // the SDAPI library has some tricks up its sleeve don't throw
    // am error unless the SDAPI does so first.
    if (token.empty())
      tokenmessage_ = "Missing access token or callback in iocontext";
    else
      tokenmessage_ = "";

    setAuthProvider(manager_.get(), token, logger_);
  }
};

/**
 * Automatically close the dataset when the last reference to it goes away.
 * This is a fallback. Please do NOT rely on this behavior. Especially if
 * the dataset is open for write. Instead, make sure xx_close() is called
 * in a timely manner. If we get here with an open dataset then:
 *
 *  - Exceptions will be logged and swallowed.
 *  - The C++ runtime might abort due to exception during unwind.
 *  - Cleanup such as making the dataset read-only might be skipped.
 */
SDGenericDatasetWrapper::~SDGenericDatasetWrapper()
{
  // It is tempting to disable the auth callback here, not below.
  // If we are called from a static destructor then the callback
  // probably isn't valid and might even crash. But if going out
  // of scope earlier this adds problems we really don't need.
  // Locks may need to be reset and the dataset might need to be
  // made readonly. That requires idtokens.
  info_.reset();
  try {
    // Explicit close so we have control of how to deal with errors. Or not.
    // Hopefully the destructors won't try to close again when they see
    // that we tried to do so ourselves. Note: currently they will try that.
    if (dataset_)
      wrapper_close(false);
  }
  catch (const std::exception& ex) {
    if (std::string(ex.what()).find("dataset is not open") == std::string::npos)
      this->logger_(0, "SDGenericDataset::close(): " + std::string(ex.what()));
  }

  // Careful: the manager indirectly gets an unsafe pointer to our
  // &tokenrefresh2_. See authorizeManager(). this->manager_ should
  // not be shared outside the SeismicStoreFile (although the wrapper
  // might in the future be shared) which means the manager should
  // never need to refresh the token after this->dataset_ has been
  // closed. But, better to be safe.
  if (tokenrefresh2_ && manager_) {
    if (!manager_.unique() && logger_(0, ""))
      logger_(0, "Manager not unique! Disable the auth callback");
    else if (logger_(1, ""))
      logger_(1, "Routinely disable the auth callback");
    manager_->setAuthProviderCallback
      (&SDGenericDatasetWrapper::tokenCallbackInvalid);
  }
  manager_.reset();
}

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> SeismicStoreFile   /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

SeismicStoreFile::SeismicStoreFile(const std::string& filename, OpenMode mode, const IOContext *iocontext)
  : FileADT()
  , _config()
{
  auto context = dynamic_cast<const OpenZGY::SeismicStoreIOContext*>(iocontext);
  _logger = ((context && context->_logger) ? context->_logger :
             LoggerBase::standardCallback
             (LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE"),
              "openzgy-cloud: ", ""));
  if (!context)
    throw OpenZGY::Errors::ZgyUserError("Opening a file from seismic store requires a SeismicStoreIOContext");
  this->_config.reset(new OpenZGY::SeismicStoreIOContext(*context));
  _rtimer.reset(new SummaryPrintingTimerEx(mode == OpenMode::ReadWrite || mode == OpenMode::Truncate ? "Cloud.reread" : "Cloud.read"));
  _wtimer.reset(new SummaryPrintingTimerEx("Cloud.write"));

  if (_logger(3, ""))
    _sslogger(3, std::stringstream() << "SeismicStoreFile("
              << "\"" << filename << "\", " << int(mode) << ", *)\n");

  std::unordered_map<std::string, std::string> extra;
  using seismicdrive::api::json::Constants;
  if (!context->_legaltag.empty())
    extra[Constants::kLegalTagLabel] = context->_legaltag;
  if (!context->_writeid.empty())
    extra[Constants::kWriteIdLabel] = context->_writeid;
  if (!context->_seismicmeta.empty())
    extra[Constants::kDDMSSeismicMetadataLabel] = context->_seismicmeta;
  //TODO-Low: if (!context->_pedantic.empty())
  //  extra[Constants::KPedantic] = context->_pedantic;

  bool sd_mgr_log = Environment::getNumericEnv("OPENZGY_SDMANAGER_LOG",0) > 0;
  bool sd_ds_log = Environment::getNumericEnv("OPENZGY_SDDATASET_LOG",0) > 0;
  auto manager = std::make_shared<seismicdrive::SDManager>
    (context->_sdurl, context->_sdapikey);
  manager->setLogStatus(sd_mgr_log);

  // TODO-Low: Cache the manager and possibly the SDUtils instance.

  auto datasetwrapper = std::make_shared<SDGenericDatasetWrapper>
    (manager, nullptr, mode, _logger);

  datasetwrapper->authorizeManager(context->_sdtoken, context->_sdtokencb2);

  std::shared_ptr<ISDGenericDataset> dataset;
  try {
    switch (mode) {
    case OpenMode::Closed:
      break;
    case OpenMode::ReadOnly:
      dataset = _open_dataset_ro(manager, filename, extra, sd_ds_log, context);
      break;
    case OpenMode::ReadWrite:
      dataset = _open_dataset_rw(manager, filename, false, extra, sd_ds_log, context);
      break;
    case OpenMode::Truncate:
      dataset = _open_dataset_rw(manager, filename, true, extra, sd_ds_log, context);
      break;
    }
  }
  catch (const std::exception& ex) {
    // Make sure this smart pointer doesn't remain alive while unwinding.
    // This interferes with the unique() test in ~SDGenericDatasetWrapper().
    manager.reset();
    datasetwrapper->throwCloudException(ex, "Open");
  }

  datasetwrapper->updateDataset(dataset);
  this->_dataset = datasetwrapper;
  // Removed this because it causes info() to be populated early,
  // thereby negating the benefit of lazy evaluation. And, worse,
  // causing different behavior when debugging is on.
  //if (_logger(1) && mode != OpenMode::Closed)
  //  _logger(1, this->_dataset->info()->toString());
}

SeismicStoreFile::~SeismicStoreFile()
{
  // The calling layer is supposed to do an explicit xx_close() so it
  // can catch and handle exceptions, and so we can be sure the token
  // callback, if used, is still valid. Do *not* try to re-authorize
  // the manager. It might not be safe to invoke the callback any
  // more. And do a blind catch of any exception because if we don't
  // the application will crash.
  if (_dataset && _dataset->dataset() && _dataset->disposition() != OpenMode::Closed) {
    try {
      _dataset->wrapper_close(_config->_set_ro_after_write);
    }
    catch (const std::exception& ex) {
      _logger(0, "EXCEPTION closing file: " + std::string(ex.what()));
    }
    catch (...) {
      _logger(0, "EXCEPTION closing file.");
    }
  }
  _dataset.reset();
}

/**
 * Convenience for invoking _logger with a simple message. Useful for
 * logging from outside this class, e.g. in the delayed-write class.
 */
bool
SeismicStoreFile::_mylogger(int priority, const std::string& message) const
{
  return _logger(priority, message);
}

/**
 * Convenience for invoking _logger with a stringstream.
 * Due to a somewhat naughty cast, the function can be caller as:
 *
 *   if(_logger(pr1))
 *    _sslogger(pri, std::stringstream() << some << data << here);
 *
 * The first line is optional. It just prevents the expression in
 * the second line from being evaluatet if debugging is disabled.
 *
 * Thread safety: Yes, because the global _logger can only be set
 * while the very first instance is being constructed.
 */
bool
SeismicStoreFile::_sslogger(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return _logger(priority, sstream ? sstream->str() : std::string());
}

/**
 * Configure the exponential backoff used by Seismic Store.
 * Only expected to be used for debugging.
 *
 * * -1 => use defaults.
 * *  0 => Turn off exponential backoff completely.
 * * >0 =? Set maximum repeat count.
 */
void
SeismicStoreFile::_set_backoff(ISDGenericDataset* sdgd)
{
  const int retries = this->_config->_retry_count;
  if (retries >= 0) {
    seismicdrive::ExponentialRetryBackoffPolicy policy;
    if (retries == 0) {
      policy.enabled = false;
    }
    else {
      policy.enabled = true;
      policy.maxRetry = retries;
      policy.initialWaitingTimeMicroSec = 500 * 1000;
      policy.maxWaitingTimeMicroSec = 32 * 1000 * 1000;
    }
    sdgd->setExponentialRetryBackoffPolicy(&policy, seismicdrive::HttpConnectionLink::ANY);
    if (_logger(2, "")) {
      std::stringstream ss;
      ss << "Backoff " << (policy.enabled ? "enabled" : "disabled")
         << " retries " << policy.maxRetry
         << " start " << (float)policy.initialWaitingTimeMicroSec*1.0e-6
         << " max "   << (float)policy.maxWaitingTimeMicroSec*1.0e-6;
      _logger(2, ss.str());
    }
  }
}

/**
 * Allocate and open a SDGenericDataset for read.
 * If dictated by the iocontext, turn on the read-only flag first.
 */
std::shared_ptr<ISDGenericDataset>
SeismicStoreFile::_open_dataset_ro(const std::shared_ptr<seismicdrive::SDManager>& manager, const std::string& filename, const std::unordered_map<std::string, std::string>& extra, bool sd_ds_log, const SeismicStoreIOContext* /*context*/)
{
  if (_logger(5, ""))
    _sslogger(5, std::stringstream()
              << "make dataset for reading using manager "
              << std::hex << (std::uint64_t)manager.get());
  auto dataset = ISDGenericDataset::createPlainInstance
    (manager.get(), filename, sd_ds_log);
  _set_backoff(dataset.get());
  dataset->open(seismicdrive::SDDatasetDisposition::READ_ONLY, extra);
  if (this->_config->_force_ro_before_read) {
    if (!dataset->getReadonlyMode()) {
      dataset->setReadonlyMode(true);
      // For robustness, assume SDAPI needs a re-open to clear the read lock.
      // For robustness, assume the dataset instance cannot be re-used.
      dataset->close();
      dataset = ISDGenericDataset::createPlainInstance
        (manager.get(), filename, sd_ds_log);
      _set_backoff(dataset.get());
      dataset->open(seismicdrive::SDDatasetDisposition::READ_ONLY, extra);
      _logger(2, "Readonly flag forced on for \"" + filename + "\"");
    }
    else {
      _logger(2, "Readonly flag already on for \"" + filename + "\"");
    }
  }
  if (_logger(5, ""))
    _sslogger(5, std::stringstream()
              << "dataset for reading is "
              << std::hex << (std::uint64_t)dataset.get());
  return dataset;
}

/**
 * Allocate and open a SDGenericDataset for write.*
 *
 * If the file is already open for write elsewhere then SDAPI will throw.
 * If the file already has a read lock this implies that the read-only
 * flag is already off, and the SDAPI will throw.
 *
 * If dictated by the iocontext, turn off the read-only flag first.
 * This change will only happen if the file is currently unlocked.
 * Otherwise the read-only flag has to be off already.
 * This might still be a bad idea. The application assumes all responsibility.
 */
std::shared_ptr<ISDGenericDataset>
SeismicStoreFile::_open_dataset_rw(const std::shared_ptr<seismicdrive::SDManager>& manager, const std::string& filename, bool truncate, const std::unordered_map<std::string, std::string>& extra, bool sd_ds_log, const SeismicStoreIOContext *context)
{
  if (_logger(5, ""))
    _sslogger(5, std::stringstream()
              << "make dataset for writing using manager "
              << std::hex << (std::uint64_t)manager.get());
  const seismicdrive::SDDatasetDisposition disp =
    truncate ?
    seismicdrive::SDDatasetDisposition::OVERWRITE :
    seismicdrive::SDDatasetDisposition::READ_WRITE;
  auto dataset = ISDGenericDataset::createPlainInstance
    (manager.get(), filename, sd_ds_log);
  if (context && context->_writethreads > 1)
    dataset = ISDGenericDataset::injectAsyncInstance
      (dataset, context->_writethreads, _logger);
  _set_backoff(dataset.get());
  if (!this->_config->_force_rw_before_write) {
    dataset->open(disp, extra);
  }
  else {
    try {
      // Unlike the open for read case, incorrect read-only state will throw.
      dataset->open(disp, extra);
      _logger(2, "Readonly flag already off for \"" + filename + "\"");
    }
    catch (const seismicdrive::SDException&) {
      // TODO-Low: A specific SDAPI exception "read-only dataset"
      // Currently a SDExceptionSDAccessorError is thrown, which is
      // more about *where* the error occured and not *what* went wrong.
      // So the catch might as well be on SDException so that fixing
      // SDAPI won't break this code. Or maybe just std::exception?
      dataset = ISDGenericDataset::createPlainInstance
        (manager.get(), filename, sd_ds_log);
      _set_backoff(dataset.get());
      // This might throw if there is a current write lock.
      dataset->open(seismicdrive::SDDatasetDisposition::READ_ONLY, extra);
      if (dataset->getReadonlyMode()) {
        dataset->setReadonlyMode(false);
        dataset->close();
        dataset = ISDGenericDataset::createPlainInstance
          (manager.get(), filename, sd_ds_log);
        _set_backoff(dataset.get());
        // Any second throw will be passed on to the caller.
        dataset->open(disp, extra);
        _logger(2, "Readonly flag forced off for \"" + filename + "\"");
      }
      else {
        _logger(2, "Readonly flag already on? for \"" + filename + "\"");
        throw;
      }
    }
  }
  if (_logger(5, ""))
    _sslogger(5, std::stringstream()
              << "dataset for writing is "
              << std::hex << (std::uint64_t)dataset.get());
  return dataset;
}

std::shared_ptr<IFileADT>
SeismicStoreFile::xx_make_instance(const std::string& filename, OpenMode mode, const IOContext *iocontext)
{
  if (filename.substr(0, 5) == "sd://" &&
      (mode != OpenMode::ReadWrite && mode != OpenMode::Truncate)) {
    auto file = std::shared_ptr<IFileADT>(new SeismicStoreFile(filename, mode, iocontext));
    // This is a no-op unless enabled by enviroment variables.
    // Note, this might have been injected after the FileParallelizer instead.
    file = FileWithPerformanceLogger::inject(file, filename);

    // Improve multi-threading of decompress and copy-out.
    auto context = dynamic_cast<const SeismicStoreIOContext*>(iocontext);
    if (context && context->_cputhreads > 1)
      file = FileParallelizer::inject(file, context->_cputhreads);

    return file;
  }
  else
    return std::shared_ptr<IFileADT>();
}

/**
 * Thread safety: Designed to be thread safe as long as the underlying
 * SDGenericDataset is. Even when data is being written in another
 * thread.
 */
void
SeismicStoreFile::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  this->_validate_read(data, offset, size, this->xx_eof(), this->_mode());
  ReadRequest request(offset, size, nullptr);
  RawList split = this->_split_by_segment(ReadList{request});
  if (this->_config->_debug_trace)
    this->_config->_debug_trace("read", /*need=*/size, /*want=*/size,/*parts*/ split.size(), this->xx_segments(true));
  for (const RawRequest& it : split) {
    // TODO-Low: port _cached_read ?
    SimpleTimerEx tt(*this->_rtimer);
    this->_dataset->dataset()->readBlock
      (static_cast<int>(it.blocknum),
       static_cast<char*>(data)+it.outpos,
       static_cast<size_t>(it.local_offset),
       static_cast<size_t>(it.local_size));
    _rtimer->addBytesRead(it.local_size);
  }
}

/**
 * Thread safety: Designed to be thread safe as long as the underlying
 * SDGenericDataset is. Even when data is being written in another
 * thread.
 *
 * TODO-Worry: even with consolidate_overlaps=false, overlapping
 * requests might cause surprises. Since this isn't supposed to
 * happen anyway, maybe just fall back to one brick at a time?
 */
void
SeismicStoreFile::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
  if (requests.size() == 1) {
    // Handle this simple case specially. There will be more cases to test
    // but the shortcut might help performance. Especially if the memory
    // allocation can be made more efficient. For testing the shortcut can be
    // made unconditional. But that will disable the consolidate-brick logic.
    // Explicitly use xx_read() in this class, not any overrides. If xx_read()
    // is overridden then whoever did that wouldn't expect xx_readv() to change.
    // The fact that one is implemented using the other is an implementation detail.
    // Note that the delivery function can retain a reference to the data.
    // This is allowed as long as the data is still short lived. If not then
    // this isn't a disaster due to the eviction code in _allocate().
    for (const ReadRequest& r : requests) {
      std::shared_ptr<void> data = _allocate(r.size);
      this->SeismicStoreFile::xx_read(data.get(), r.offset, r.size, usagehint);
      _deliver(r.delivery, data, 0, r.size, transient_ok);
    }
    return;
  }

  // Consolidate adjacent bricks before reading.
  //
  // Remember to use a nonvirtual call when xx_readv() is implemented
  // in terms of xx_read() or vice versa. Because this is an
  // implementation detail and overriding one of them should not
  // affect the other. Similarly use a nonvirtual xx_eof() because we
  // are not supposed to know that SeismicStoreFileDelayedWrite makes
  // the file look bigger.
  //
  // This implementation can issue requests in multiple
  // threads, wait for all threads to complete, and then deliver all
  // the results. For this reason it needs to allocate a buffer to
  // hold the entire data to be read.
  //
  // TODO-Performance: Allow read from cloud and copy-out/decompress
  // in parallel inside a single request from the application. Probably
  // not worth the (significant) trouble, and probably won't help
  // multi-threaded applications anyway. Theoretically it might help
  // lowres computation on write. The main caveat is that requests that
  // cross segment boundaries will then need some complicated "partial
  // delivery" mechanism. Or the entire request may need to fall back
  // to the original implementation if a boundary crossing, which is
  // likely to be very rare, is detected.

  std::int64_t current_eof = SeismicStoreFile::xx_eof(); // exclude open segment
  _validate_readv(requests, current_eof, this->_mode());

  // For debugging / logging only
  const std::int64_t asked =
    std::accumulate(requests.begin(), requests.end(), std::int64_t(0),
                    [](std::int64_t a, const ReadRequest& b) {
                      return a + b.size;
                    });

  ReadList new_requests = ConsolidateRequests::consolidate
    (requests, _config->_maxhole, _config->_maxsize, _config->_aligned,
     false/*consolidate_overlaps*/, current_eof);
  auto work = _split_by_segment(new_requests);

  // Carefully get the required buffer size. Normally it would be
  // enough to just look at work.back() but there may be some odd
  // corner cases, and I just don't want to worry about those.

  const std::int64_t realsize =
    std::accumulate(work.begin(), work.end(), std::int64_t(0),
                    [](std::int64_t a, const RawRequest& b) {
                      return std::max(a, b.local_size + b.outpos);
                    });

  // This would probably work, but the single-brick case is already
  // handled and the case for two to four int8 bricks or two int16
  // bricks are not that interesting. At least not for applications
  // that read just one brick at a time. Those apps will not get here.
  //std::shared_ptr<void> data = _allocate(r.size);
  std::shared_ptr<char> data(new char[realsize], std::default_delete<char[]>());

  if (this->_config->_debug_trace)
    this->_config->_debug_trace("readv", /*need=*/asked, /*want=*/realsize,/*parts*/ work.size(), this->xx_segments(true));

  // Do the actual reading of the consolidated chunks, possibly using
  // multiple threads.
  //
  //  * Worry: Can there be multiple requests targeting the same area
  //    of the output buffer? Probably not although there can be multiple
  //    read requests for the same area of the file.
  //
  // If parallel_ok, can I then deliver data as it is received without
  // waiting for the last bit? That allows reading and e.g. decompressing
  // in parallel. Not a big deal if application has multiple reads in
  // flight. Otherwise this might in theory double the speed.
  //
  //  * Tricky to implement. Receiver doesn't allow partial delivery.
  //    So if one request requires concatenating data from multiple
  //    cloud reads then this needs to wait until the end. Or of really
  //    fancy, keep track of when all the data has need read for each
  //    of the original requests.
  const std::int64_t worksize = work.size();
  const std::int64_t threadcount = std::max(std::min(std::min(
      worksize,
      static_cast<std::int64_t>(omp_get_max_threads())),
      _config->_iothreads),
      static_cast<std::int64_t>(1));
  MTGuard guard("cloud-read", (int)threadcount);
  //std::cerr << "Access seismic store (" << worksize << "): ";
#pragma omp parallel for num_threads((int)threadcount) schedule(dynamic,1)
  for (std::int64_t ii=0; ii<worksize; ++ii) {
    //if (!ii) std::cerr << ("[" + std::to_string(omp_get_num_threads()) + "]");
    const auto& it = work[ii];
    guard.run([&](){
      //std::cerr << "0123456789"[omp_get_thread_num() % 10];
      SimpleTimerEx tt(*_rtimer);
      this->_dataset->dataset()->readBlock
        (static_cast<int>(it.blocknum),
         data.get() + it.outpos,
         static_cast<size_t>(it.local_offset),
         static_cast<size_t>(it.local_size));
    });
    _rtimer->addBytesRead(it.local_size);
  }
  guard.finished();
  //std::cerr << "$\n";

  // Do not try to multi-thread the following loop. Instead inject a
  // FileParallelizer instance at a higher level. At this lowest level,
  // each request might be a large consolidated read. Splitting and
  // parallelizing the CPU bound tasks should be done at a finer one
  // brick granularity and FileParallelizer is designed to do just that.

  std::int64_t pos = 0;
  for (const ReadRequest& rr : new_requests) {
    std::int64_t this_size = std::max((std::int64_t)0, std::min(rr.size, current_eof - rr.offset));
    // TODO-Worry: If this_size != rr.size, can this ever happen?
    // If yes then we might have lost track of where in the buffer
    // we should copy out from. This wory also applies to the Python code.
    _deliver(rr.delivery, data, pos, this_size, transient_ok);
    pos += this_size;
  }
}

/**
 * Thread safety: NO. The method would need more information, such as
 * the size of all previous segments "in transit", in order to do this
 * correctly. On the other hand there is nothing preventing us to
 * split a segment in parts inside this method and write those parts
 * in parallel.
 */
void
SeismicStoreFile::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  SimpleTimerEx tt(*_wtimer);
  this->_validate_write(data, offset, size, this->_mode());
  std::int64_t current_eof = SeismicStoreFile::xx_eof(); // MUST be nonvirtual
  if (_logger(5, ""))
    _sslogger(5, std::stringstream()
              << "SeismicStoreFile.xx_write("
              << "offset=" << offset
              << ", size=" << size
              << ", current EOF is " << current_eof
              << ")\n");
  std::int64_t blocknum{0}, local_offset{0}, local_size{0};
  bool overwrite{false};
  if (offset == current_eof) {
    // Sequential write from current EOF. Create a new segment.
    blocknum = this->_dataset->info()->blockCount();
  }
  else if (offset < current_eof) {
    // Rewrite existing block. Resizing not allowed.
    overwrite = true;
    this->_dataset->info()->getLocalOffset
      (offset, size, &blocknum, &local_offset, &local_size);
    // Normally we only get here to overwrite blob 0, and that is ok.
    // Writing block 0 is not multi-threaded and does not resize.
    // If opening an existing file for update it depends on how that
    // is handled elsewhere. Hopefully we still won't get here
    // If that happens then there are several caveats:
    // - May need to allow resizing the last brick, which in turn
    //   invalidates some assumptions about immutable information.
    // - The test below will fail in the parallel upload case
    //   because local_offset and local_size refers to SDAPI blocks and
    //   not the larger segments that we are asked to write. local_size
    //   will usually not be larger than one SDAPI block and will thus
    //   fail the size check.
    // - Maybe check offset+N*(segsize/segsplit) (last SDAPI block)?
    // - I am unsure whether it is only the test that is wrong or whether
    //   this case needs more special handling.
    if (local_offset != 0 || local_size != size)
      throw OpenZGY::Errors::ZgyInternalError("Cannot write resized segment.");
  }
  else {
    throw OpenZGY::Errors::ZgyUserError("Cannot write segments out of order.");
  }

  // The segsize configuration belongs to SeismicFileDelayedWrite but is
  // needed here because the blobsize configuration that belongs to us
  // isn't stored explicitly. It is stored as segsize/segsplit. The same
  // configuration instance is passed to both us and SeismicFileDelayedWrite
  // so getting it here is easy even though it technically breaks isolation.
  //
  // Side issue: If some code creates a SeismicFile instance that is to be
  // used directly without a SeismicStoreFileDelayedWrite wrapper then segsplit
  // should be set to 1. I am not sure what happens if it isn't.
  if (_config->_segsize > 0 && _config->_segsize % _config->_segsplit != 0)
    throw OpenZGY::Errors::ZgyUserError("segsize must be a multiple of segsplit");
  const std::int64_t blobsize =
    (blocknum == 0 || _config->_segsize <= 0 ||_config->_segsplit <= 1) ?
    size : _config->_segsize / _config->_segsplit;
  if (size <= blobsize) {
    do_write_one(data, blocknum, size, overwrite);
  }
  else {
    do_write_many(data, blocknum, size, blobsize, overwrite);
  }
  if (this->_config->_debug_trace)
    this->_config->_debug_trace
      (offset == current_eof ? "append" : "write",
       size, size, 1, this->xx_segments(true));
}

OpenMode
SeismicStoreFile::_mode() const
{
  return _dataset ? _dataset->disposition() : OpenMode::Closed;
}

/**
 * This is the final part of xx_write, in the case where we want to
 * write just a simgle SDAPI blob.
 */
void
SeismicStoreFile::do_write_one(
     const void*  const data,
     const std::int64_t blocknum,
     const std::int64_t size,
     const bool         overwrite)
{
  if (_logger(1, "")) {
    _sslogger(1, std::stringstream()
              << "do_write_one(*, " << blocknum << ", " << size << ", "
              << std::boolalpha << overwrite << ")");
    }
  this->_dataset->info()->checkOnWrite(blocknum, size);
  this->_dataset->dataset()->writeBlock
    (static_cast<int>(blocknum), // cast checked by checkOnWrite()
     static_cast<const char*>(data),
     static_cast<std::size_t>(size),
     overwrite);
  _wtimer->addBytesWritten(size);
  this->_dataset->updateOnWrite(blocknum, size);
}

/**
 * This is the final part of xx_write, in the case where we want to
 * write two or more SDAPI blobs of no more than blobsize bytes each.
 * It should work for a single blob as well, but that would be a
 * major overkill compared to just calling do_wrire_one instead.
 *
 * Possible scenarios.
 *
 * \li We never get here for block 0; than one cannot be split and it
 *     would rarely make sense anyway because block 0 is usually small.
 *     Calling this code anyway with blobsize == size and blobcount == 1
 *     should work but is needlessly roundabout and I have not checked
 *     that it won't hit any corner cases.
 *
 * \li If size <= blobsize then we might as well use the old code,
 *     even if the other condtions for writing in parallel pass.
 *     There will be just one SDAPI block and only one thread will
 *     be used. The short cut means slightly less code coverage
 *     when running unit tests with only small files. So, use
 *     bigger ones.
 *
 * \li In most cases we get segsize bytes here (typically 1 GB) and
 *     we are expected to write out the data in segsplit (typically 8)
 *     identical SDAPI blocks each of size segsize/segsplit.
 *
 * \li If size != segsize this is the last or
 *     (if -1 <= blocknum <= segsplit) only segment with bulk data.
 *     In that case we may end up writing fewer than segsplit SDAPI
 *     blocks. Possibly just one. And the last SDAPI block can be
 *     smaller than segsize/segsplit blocks. We can even have
 *     size < segsize/segsplit meaning all headers are in SDAPI
 *     block 0 and all data in SDAPI block 1. In the latter case
 *     there will be no parallelization.
 *
 * The last two bullets seem to suggest we need to know what segsize is,
 * But no, all this code needs to know is that it shall split the
 * incoming buffer into chunks of no more than blobsize bytes.
 */
void
SeismicStoreFile::do_write_many(
     const void*  const data,
     const std::int64_t blocknum,
     const std::int64_t size,
     const std::int64_t blobsize,
     const bool         overwrite)
{
  // Only need to run consistency checks on the the first of the
  // segsplit SDAPI blocks written. So we end up checing block
  // 0, 1, segsplit+1, ... When checking block <2 the real segment
  // size hasn't been established yet so there isn't much to check.
  // Otherwise the check asserts that _dataset->info()->block1size()
  // equals blobsize. There is always a check that blocknum fits in
  // an int. To make the cast safe.
  this->_dataset->info()->checkOnWrite(blocknum, std::min(size, blobsize));

  const int blobcount = static_cast<int>((size + blobsize - 1) / blobsize);
  MTGuard guard("cloud-write", blobcount);
#pragma omp parallel for num_threads(blobcount)
  for (int ii = 0; ii < blobcount; ++ii) {
    if (/*blocknum == 1 &&*/ ii == 0 && _logger(1, "")) {
      _sslogger(1, std::stringstream()
                << "do_write_many(*, "
                << blocknum << ".." << (blocknum + blobcount - 1) << ", "
                << size << ", "
                << blobsize << ", "
                << std::boolalpha << overwrite << ")"
                << " using " << omp_get_num_threads() << " threads");
    }
    guard.run([&](){
      const int          iter_blocknum  = ii + static_cast<int>(blocknum);
      const std::int64_t iter_offset    = ii * blobsize;
      const std::int64_t iter_endoffset = std::min(size, (ii+1) * blobsize);
      std::int64_t       iter_size      = iter_endoffset - iter_offset;
      const char*        iter_data = static_cast<const char*>(data)+iter_offset;
      if (iter_size > 0) {
        this->_dataset->dataset()->writeBlock
          (iter_blocknum, iter_data, iter_size, overwrite);
        _wtimer->addBytesWritten(iter_size);
      }
    });
  }
  guard.finished();

  // The rest of the functions, including updateOnWrite(), are unaware
  // that we might have done parallel writes. Keep them in the dark.
  // Otherwise there will be errors reported about writing out of sequence.
  // Not to mention that updateOnWrite() isn't threadsafe.
  for (int ii = 0; ii < blobcount; ++ii) {
    const int          iter_blocknum  = ii + static_cast<int>(blocknum);
    const std::int64_t iter_offset    = ii * blobsize;
    const std::int64_t iter_endoffset = std::min(size, (ii+1) * blobsize);
    std::int64_t       iter_size      = iter_endoffset - iter_offset;
    if (iter_size > 0) {
      this->_dataset->updateOnWrite(iter_blocknum, iter_size);
    }
  }
}

/**
 * \details: Thread safety: No. All other operations must be completed first.
 */
void
SeismicStoreFile::xx_close()
{
  if (!_dataset || _dataset->disposition() == OpenMode::Closed) {
    // Note: I might "be nice" to the application and simply ignore a duplicate
    // close or a close on a file that was never open in the first place.
    // In that case I should probably check using an atomic_flag.
    // But if the application issues extraneous xx_close(), let alone multiple
    // concurrent calls to xx_close(), this is a bug. That may indicate there
    // is something else wrong as well.
    throw OpenZGY::Errors::ZgyUserError("Attemping to close a file twice.");
  }

  if (_dataset) {
    try {
      switch (_dataset->disposition()) {
      case OpenMode::Closed:
        break;
      case OpenMode::ReadOnly:
      case OpenMode::ReadWrite:
      case OpenMode::Truncate:
        if (_dataset->dataset()) {
          auto victim = _dataset;
          _dataset.reset();
          victim->wrapper_close(this->_config->_set_ro_after_write);
        }
        break;
      }
    }
    catch (const std::exception& ex) {
      // Too late; _dataset is null so we don't know the tokenmessage.
      //_dataset->throwCloudException(ex, "Close");
      throw OpenZGY::Errors::ZgyInternalError
        ("Close: Seismic Store: " + std::string(ex.what()));
    }
  }
  _dataset.reset();
  _rtimer.reset();
  _wtimer.reset();
}

/**
 * Thread safety: Yes.
 *
 * If info() is called on one thread while another thread is doing a
 * write operation (which is actually not allowed) then it is
 * unspecified whether the returned value is before or after the
 * write.
 *
 * CAVEAT: Note that the method is overridden in
 * SeismicStoreFileDelayedWrite, and that one might not be safe.
 */
std::int64_t
SeismicStoreFile::xx_eof() const
{
  if (!_dataset) {
    throw OpenZGY::Errors::ZgyUserError("The file is not open in xx_eof().");
    //return -1; Might be safer if this is happening in a destructor.
  }
  return _dataset->info()->totalSize();
}

/**
 * \brief Return the size of each segment of the file.
 * \details: Thread safety: Not if writes may be in progress. Could be fixed.
 *
 * If complete=false return at most 3 numbers: The first, second,
 * and last segment size. Currently all segments except the
 * first and last are required to have the same size, so by
 * combining the results of xx_segments() and xx_eof() it is
 * possible to compute the rest of the information.
 */
std::vector<std::int64_t>
SeismicStoreFile::xx_segments(bool complete) const
{
  if (!this->_dataset)
    return std::vector<std::int64_t>{};
  return this->_dataset->info()->allSizes(complete);
}

/**
 * \details: Thread safety: Yes.
 */
bool
SeismicStoreFile::xx_iscloud() const
{
  return true;
}

/**
 * \details: Thread safety: Yes.
 */
void
SeismicStoreFile::deleteFile(const std::string& filename, bool missing_ok) const
{
  if (_logger(2, ""))
    _sslogger(2, std::stringstream()
              << "SeismicStoreFile::deleteFile("
              << "\"" << filename << "\", "
              << "missing_ok=" << std::boolalpha << missing_ok
              << ")\n");
  if (!_dataset)
    throw OpenZGY::Errors::ZgyUserError("The manager is not open in deleteFile.");
  // Make sure the returned smart pointer doesn't go out of scope.
  std::shared_ptr<seismicdrive::SDManager> smart_manager = _dataset->manager();
  seismicdrive::SDUtils utils(smart_manager.get());
  try {
    utils.deleteDataset(filename);
  }
  catch (const std::exception& ex) {
    if (std::string(ex.what()).find("does not exist") != std::string::npos) {
      if (_logger(1, ""))
        _sslogger(1, std::stringstream()
                  << "Deleting already deleted"
                  << " \"" << filename << "\"\n");
      if (!missing_ok)
        _dataset->throwCloudException(ex, "Delete");
    }
    else {
      _dataset->throwCloudException(ex, "Delete");
    }
  }
  _logger(2, "SeismicStoreFile::deleteFile DONE.\n");
}

std::string
SeismicStoreFile::altUrl(const std::string& filename) const
{
  if (!_dataset)
    throw OpenZGY::Errors::ZgyUserError("The manager is not open in altUrl.");

  // Should I strip off any "?context= first? It doesn't make sense
  // to create an alturl from another alturl. Probably doesn't
  // matter much either way.

  try {
    // Make sure the returned smart pointer doesn't go out of scope.
    std::shared_ptr<seismicdrive::SDManager> smart_manager = _dataset->manager();
    std::shared_ptr<ISDGenericDataset> dataset =
      // The ugly const-cast reflects that opening for read will in this
      // case have a permanent effect on the file. There is a reason the
      // flag starts with "force".
      const_cast<SeismicStoreFile*>(this)->
      _open_dataset_ro(smart_manager, filename,
                       std::unordered_map<std::string, std::string>(), false, nullptr);
    const std::string wid = dataset->getConsistencyID();
    const std::string ctx = dataset->getSerializedContext();
    const std::string url = filename.substr(0, filename.find("?")) + "?context=" + ctx;
    try {
      dataset->close();
    }
    catch (const std::exception& ex) {
      // Workaround for what is probably a bug in SDAPI.
      // Hopefully the dataset actually got closed anyway,
      // so we don't get a resource leak here.
      if (!strstr(ex.what(), "has been locked with different ID"))
        throw;
      _logger(0, "getSerializedContext() caused bogus exception for wid "+ wid + ": " + std::string(ex.what()));
    }
    if (_logger(2, ""))
      _sslogger(2, std::stringstream()
                << "SeismicStoreFile::altUrl(\""
                << (filename.size() < 76 ? filename : filename.substr(0, 72) + "...")
                << "\")"
                << " = \"" << "REDACTED"/*url*/ << "\"" // Don't log secrets!
                << "\n");
    return url;
  }
  catch (const seismicdrive::error::dataset::context::NotReadOnly& ex) {
    throw OpenZGY::Errors::ZgyNotReadOnlyError(ex.what());
  }
  catch (const std::exception& ex) {
    _dataset->throwCloudException(ex, "altUrl");
    throw; // not reached, but compiler might not realize that.
  }
}

std::string
SeismicStoreFile::idToken() const
{
  if (!_dataset)
    throw OpenZGY::Errors::ZgyUserError("The manager is not open in idToken.");
  try {
    // Make sure the returned smart pointer doesn't go out of scope.
    std::shared_ptr<seismicdrive::SDManager> manager = _dataset->manager();
    std::string token = manager->getIDToken();
    return token.substr(0, 7) == "Bearer " ? token.substr(7) : token;
  }
  catch (const std::exception& ex) {
    _dataset->throwCloudException(ex, "altUrl");
    throw; // not reached, but compiler might not realize that.
  }
}

/**
 * Given one or more (offset, size, ...) tuples, convert these
 * to (segment_number, offset_in_seg, size_in_seg, outpos).
 * The delivery functor in the read requests is ignored.
 *
 * "outpos" is the offset to store the data that was read, if
 * it is to be stored sequentially in one large buffer.
 *
 * Request for data past EOF will throw an exception.
 *
 * The returned list might be longer than the input if any of the
 * input requests crossed segment boundaries.
 * The return might be shorter than the input or even empty if
 * any input request was for 0 bytes.
 *
 * The algorithm is O(n^2) on segment_count * request_count
 * but both numbers should be small. If this actually becomes
 * a problem then use binary search in self._cumsize to find
 * the starting segment.
 *
 * Maybe simplify: This logic could be moved inside SDAPI or the
 * SDAPI wrapper. Reads from segment "-1" imply linear access.
 * There would be a slight change in that requests split due to
 * crossing a segment boundary would not be parallelized. But
 * that is expected to be a very rare occurrence.
 *
 * Thread safety: Only if no writes may be pending.
 */
SeismicStoreFile::RawList
SeismicStoreFile::_split_by_segment(const ReadList& requests)
{
  RawList result;
  std::int64_t outpos = 0;
  for (const ReadRequest& rr : requests) {
    std::int64_t offset{rr.offset}, size{rr.size};
    std::int64_t blocknum{0}, local_offset{0}, local_size{0};
    while (size > 0) {
      // Normally just one iteration i.e. no crossing seg boundary.
      this->_dataset->info()->getLocalOffset
        (offset, size, &blocknum, &local_offset, &local_size);
      result.push_back(RawRequest
                       (blocknum, local_offset, local_size, outpos));
      offset += local_size; // file offset, starting point in read request.
      outpos += local_size; // buffer offset, all read requests concatenated.
      size -= local_size;
    }
  }
  return result;
}

void
SeismicStoreFile::_cached_read(/*TODO-Low: seg, offset, view*/)
{
  throw std::runtime_error("SeismicStoreFile::_cached_read: not implemented yet");
}

/////////////////////////////////////////////////////////////////////////////
//    FileADT -> SeismicStoreFile -> SeismicStoreFileDelayedWrite   /////////
/////////////////////////////////////////////////////////////////////////////

SeismicStoreFileDelayedWrite::SeismicStoreFileDelayedWrite(const std::string& filename, OpenMode mode, const IOContext *iocontext)
  : FileADT()
  , _config(nullptr)
  , _relay()
  , _open_segment()
  , _usage_hint(UsageHint::Unknown)
{
  this->_ctimer.reset(new SummaryPrintingTimerEx("Cloud.readcache"));
  this->_relay.reset(new SeismicStoreFile(filename, mode, iocontext));

  // The relayed file already did this so we are making another copy.
  // Not a big deal.
  auto context = dynamic_cast<const OpenZGY::SeismicStoreIOContext*>(iocontext);
  if (!context)
    throw OpenZGY::Errors::ZgyUserError("Opening a file from seismic store requires a SeismicStoreIOContext");
  this->_config.reset(new OpenZGY::SeismicStoreIOContext(*context));

  try {
    if (mode == OpenMode::ReadWrite)
      this->_reopen_last_segment();
  }
  catch (const std::exception& ex) {
    _relay->datasetwrapper()->throwCloudException(ex, "Open write");
  }
}

SeismicStoreFileDelayedWrite::~SeismicStoreFileDelayedWrite()
{
  try {
    this->_flush(true);
  }
  catch (const std::exception& ex) {
    // The calling layer is supposed to do an explicit xx_close()
    // so it can catch and handle exceptions. This blind catch is
    // just a desperate attempt to avoid an application crash.
    _relay->_mylogger(0, "EXCEPTION flushing file: " + std::string(ex.what()));
  }
  // Note: The dataset itself will be closed in _relay's destructor.
  // That should happen very shortly.
}

std::shared_ptr<IFileADT>
SeismicStoreFileDelayedWrite::xx_make_instance(const std::string& filename, OpenMode mode, const IOContext *iocontext)
{
  if (filename.substr(0, 5) == "sd://" &&
      (mode == OpenMode::ReadWrite || mode == OpenMode::Truncate))
    return std::shared_ptr<IFileADT>(new SeismicStoreFileDelayedWrite(filename, mode, iocontext));
  else
    return std::shared_ptr<IFileADT>();
}

/**
 * Reads consist of two parts: Data before the last committed segment
 * which are handled by the relay and data still in memory in the open segment.
 * Most reads will be in just one of the regions but this is not required.
 *
 * Thread safety: No writes may be pending. Designed to be thread safe
 * when only reading.
 */
void
SeismicStoreFileDelayedWrite::xx_read(void *data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  this->_validate_read(data, offset, size, this->xx_eof(), this->_relay->_mode());
  const std::int64_t closed_size =
    std::max(std::int64_t(0), std::min(size, this->_relay->xx_eof() - offset));
  const std::int64_t opened_size =
    size - closed_size;
  const std::int64_t local_offset =
    std::max(std::int64_t(0), offset - this->_relay->xx_eof());

  if (local_offset + opened_size > static_cast<std::int64_t>(this->_open_segment.size()))
    throw OpenZGY::Errors::ZgyUserError("Trying to read past EOF");
  if (closed_size > 0)
    this->_relay->xx_read(data, offset, closed_size, this->_usage_hint);
  if (opened_size > 0) {
    // Timing of memcpy is not interesting, but the number of calls
    // and the total byte count might be.
    SimpleTimerEx tt(*this->_ctimer);
    memcpy(static_cast<char*>(data) + closed_size,
           this->_open_segment.data() + local_offset,
           opened_size);
    this->_ctimer->addBytesRead(opened_size);
  }
}

/**
 * Thread safety: Inherited from relay->xx_readv and this->xx_read().
 * No writes may be pending.
 */
void
SeismicStoreFileDelayedWrite::xx_readv(const ReadList& requests, bool parallel_ok, bool immutable_ok, bool transient_ok, UsageHint usagehint)
{
  std::int64_t end{0};
  for (const ReadRequest& rr : requests) {
    end = std::max(end, rr.offset + rr.size);
  }
  if (end <= this->_relay->xx_eof()) {
    // The open segment is not involved, so just forward the request.
    this->_relay->xx_readv
      (requests, parallel_ok, immutable_ok, transient_ok, usagehint);
  }
  else {
    // Let xx_read handle the requests one at a time. If the requests
    // consist of both open and closed segments then this is not
    // efficient since SD access won't be paralellized. But that case
    // would be a lot of effort to support and it won't happen often.
    // Using a nonvirtual call when xx_readv() is implemented in terms
    // of xx_read() or vice versa. Because this is an implementation
    // detail and overriding one of them should not affect the other.
    for (const ReadRequest& r : requests) {
      std::shared_ptr<char> data(new char[r.size], std::default_delete<char[]>());
      this->SeismicStoreFileDelayedWrite::xx_read(data.get(), r.offset, r.size, usagehint);
      _deliver(r.delivery, data, 0, r.size, transient_ok);
    }
  }
}

/**
 * Write data to seismic store, buffering the writes to get larger
 * segment sizes. Writes are only allowed at offset 0 and at EOF.
 * This is less general then the parent type which lets us rewrite
 * any segment as long as its size does not change.
 *
 * The following cases are allowed:
 *
 *   - Write or overwrite segment 0. This segment will not be buffered
 *     and may have a different size than the rest. If already written
 *     then the size cannot change. Note: ZGY stores all its header
 *     data in segment 0. ZGY will need to write headers in a single
 *     operation.
 *
 *   - Overwrite an entire segment that has already been written.
 *     The size cannot change. This might be useful in some very
 *     special situations.
 *
 *   - Overwrite an arbitrary part of the data that has been written
 *     but not flushed to seismic store. This enables some fuzzy logic
 *     used to implement read/modify/write in OpenZGY.
 *
 *   - Append data at EOF. Data is flushed to disk when the open
 *     segment grows larger than the declared segment size, or on
 *     file close.
 *
 * The code does not support a combination of these cases. I.e the
 * data cannot span the boundary between flushed and open segments,
 * and cannot cover both before and after EOF.
 * TODO-Low: Support can be added if it ever turns out to be useful.
 * Possibly this may happen if/when support is added for update.
 * The only scenarios used today are overwrite entire segment 0
 * which will then always be closed. And append at EOF which will
 * obviously then not have date both before and after EOF and will
 * throw ZgySegmentIsClosed if data spans the open and closed parts.
 *
 * If segsize is zero no buffering is done and each write will either
 * create a new segment or completely rewrite an existing segment.
 *
 * Thread safety: NO. The method would need more information, such as
 * the size of all previous segments "in transit", in order to do this
 * correctly. On the other hand there is nothing preventing us from
 * splitting a segment in parts inside this method and write those parts
 */
void
SeismicStoreFileDelayedWrite::xx_write(const void* data, std::int64_t offset, std::int64_t size, UsageHint usagehint)
{
  const std::int64_t written = this->xx_eof();
  const std::int64_t committed = this->_relay->xx_eof();

  if (_relay->_mylogger(5, ""))
    _relay->_sslogger(5, std::stringstream()
            << "SeismicStoreFileDelayedWrite.xx_write("
            << "offset=" << offset
            << ", size=" << size
            << "), current EOF is " << written
            << ", segment size " << this->_config->_segsize
            << ", split into " << this->_config->_segsplit
            << "\n");

  if (offset != 0 && this->_config->_segsize != 0 && offset < committed) {
    // Updating is only allowed at offset 0 (headers) and in the memory buffer.
    // Callers may catch this specific exception and handle it by allocating
    // a new block and leaking the old. This is why there is an explicit check.
    // Another reason is that the lower levels are slightly more permissive.
    // In 99% of the cases a "trying to change the size of a segment" error
    // will be reported. But there are some odd corner cases that might get
    // thru because the segment size happened to match. While technically
    // legal this behavior would be really inconsistent.
    // TODO-Low: If caller doesn't catch ZgySegmentIsClosed then the memory
    // buffer should also be off limits to avoid heissenbugs. But as of today
    // the exception is always caught and handled so this might be academic.
    if (_relay->_mylogger(1, ""))
      _relay->_sslogger(1, std::stringstream() << std::hex
                        << "Cannot write at 0x" << offset
                        << " because data up to 0x" << committed
                        << " has already been flushed to the cloud\n");
    throw OpenZGY::Errors::ZgySegmentIsClosed("Block has already been flushed.");
  }

  if (offset == 0 || this->_config->_segsize <= 0 || offset < committed) {
    this->_relay->xx_write(data, offset, size, usagehint);
    if (this->_config->_debug_trace)
      this->_config->_debug_trace("flush", size, size, 1, this->xx_segments(true));
    return;
  }

  // The reason I deferred this validation is that if forwarding directly
  // to the relay it isn't needed. Because the relay does the same check.
  this->_validate_write(data, offset, size, this->_relay->_mode());
  if (offset > written)
    throw OpenZGY::Errors::ZgyUserError("Data must be written sequentially.");

  // TODO-Low: Generalize: If caller doesn't catch ZgySegmentIsClosed
  // then all rewrites ought to be forbidden to avoid heissenbugs. Since
  // ZGY is our only client and does in fact catch that exception then
  // this is low priority to change.
  if (offset != written) {
    if (_relay->_mylogger(1, ""))
      _relay->_sslogger(1, std::stringstream()
              << "Write at " << offset << " which is neither 0 nor EOF"
              << " when EOF is " << this->xx_eof()
              << ". Also, size is " << size
              << " and current open segment has " << this->_open_segment.size()
              << " bytes.\n");
    //throw OpenZGY::Errors::ZgyInternalError("Invalid offset passed to xx_write()");
  }
  if (offset == written) {
    // Append data to open segment
    const char *cdata = static_cast<const char*>(data);
    this->_open_segment.insert(this->_open_segment.end(), cdata, cdata + size);
  }
  else if (offset + size <= written) {
    // Update data fully inside open segment.
    // I already know it isn't in the committed segments as it was checked above.
    memcpy(this->_open_segment.data() + (offset-committed), data, size);
  }
  else {
    // Part update the open segment, part new data appended to it.
    // Not supported. See comment at start of function.
    throw OpenZGY::Errors::ZgyInternalError("Partial overwrite/append not supported");
  }
  if (this->_usage_hint == UsageHint::Unknown)
    this->_usage_hint = usagehint;
  else if (this->_usage_hint != usagehint)
    this->_usage_hint = UsageHint::Mixed;

  if (this->_config->_debug_trace)
    this->_config->_debug_trace("queue", size, size, 1, this->xx_segments(true));

  this->_flush(false);
}

/**
 * \details: Thread safety: No. All other operations must be completed first.
 */
void
SeismicStoreFileDelayedWrite::xx_close()
{
  this->_flush(true);
  this->_ctimer.reset();
  this->_relay->xx_close();
}

/**
 * \details: Thread safety: Not if writes may be in progress. Could be fixed.
 */
std::int64_t
SeismicStoreFileDelayedWrite::xx_eof() const
{
  return (this->_relay->xx_eof() +
          static_cast<std::int64_t>(this->_open_segment.size()));
}

/**
 * \brief Return the size of each segment of the file.
 * \details: Thread safety: Not if writes may be in progress. Could be fixed.
 *
 * If complete=false return at most 3 numbers: The first, second,
 * and last segment size. Currently all segments except the
 * first and last are required to have the same size, so by
 * combining the results of xx_segments() and xx_eof() it is
 * possible to compute the rest of the information.
 */
std::vector<std::int64_t>
SeismicStoreFileDelayedWrite::xx_segments(bool complete) const
{
  std::vector<std::int64_t> result = this->_relay->xx_segments(complete);
  if (!complete && result.size() >= 3)
    result.resize(2);
  result.push_back(this->_open_segment.size());
  return result;
}

/**
 * \details: Thread safety: Yes.
 */
bool
SeismicStoreFileDelayedWrite::xx_iscloud() const
{
  return this->_relay->xx_iscloud();
}

/**
 *
 * Re-open the last segment by reading it into memory. This is
 * needed if the last segment was not full. In the unlikely case
 * where the last write managed to exactly fill the segment it isn't
 * technically needed but is done anyway to avoid corner cases.
 *
 * Note that if there is just one segment there is nothing more
 * to be done. Our callers are responsible for reading and parsing
 * the headers on open and writing them back on close.
 *
 * Do some consistency checks, over and above what the ZgyWriter
 * does to make sure this file was written by OpenZGY and that
 * critical parameters such as the segment size was not changed.
 *
 * Instead of throwing an exception the code could have silently changed
 * segsize to a valid number. But there would be several caveats.
 *
 * - The underlying SeimsicStoreFile needs the segment size in xx_write()
 *   and its value might have been cached in the constructor.
 *
 * - If segsize needs to change, what should the code do with segsplit?
 *   Regardless of what is chosen it won't be intutive.
 *
 * - Some cases might end up with odd segment size. E.g. if there are
 *   two segments, original segment size 2 GB, new segment size 1 GB
 *   and the second (not full) segment is 1025 MB then that would
 *   become the new segment size.
 *
 * Called from constructor, so calling virtual methods in this class
 * won't work.
 *
 * Catching and translating exceptions from SDAPI is done by caller.
 *
 * Thread safety: No.
 */
void
SeismicStoreFileDelayedWrite::_reopen_last_segment()
{
  const std::int64_t user_segsize = this->_config->_real_segsize;
  const std::vector<std::int64_t> segments = this->_relay->xx_segments(false);
  std::shared_ptr<SDGenericDatasetWrapper> wrapper =
    this->_relay->datasetwrapper();
  const std::int64_t numseg =
    static_cast<std::int64_t>(wrapper->dataset()->getBlockNum());
  const std::int64_t numbytes = this->_relay->xx_eof();
  if (numseg >= 2) {
    //const std::int64_t file_segsize  = dataset->getBlockSize(1);
    //const std::int64_t file_lastsize = dataset->getBlockSize(numseg-1);
    // Or trust the xx_segments
    const std::int64_t file_segsize  = segments[1];
    const std::int64_t file_lastsize = segments.back();

    // CHECK the user supplied segment size did not change
    if (numseg == 2) {
      // The second and also last segment contains data but is
      // probably not full.
      if (file_segsize > user_segsize) {
        throw OpenZGY::Errors::ZgyUpdateRules("This ZGY file was written with a larger segment size, or it was not uploaded by OpenZGY");
        // If segments[1] <= segsize*segsplit this might technically have
        // worked, but would be a really obscure corner case causing the
        // file to be rewritten with a different segment size.
      }
    }
    else { // numseg > 2
      // The second segment is now full, so this must have been the
      // segment size the file was origially written with.
      if (file_segsize != user_segsize) {
        throw OpenZGY::Errors::ZgyUpdateRules("This ZGY file was written with a different segment size, or it was not uploaded by OpenZGY");
      }
    }

    // FOOL both SeismicStoreFile and SeismicStoreFileDelayedWrite to
    // believe we are still writing the file for the first time,
    // with no way except for the OpenMode to say it is wrong.
    //
    // Need to re-open the last segment by reading it into the open
    // buffer and deleting it from the cloud file.
    std::vector<char>& seg = this->_open_segment;
    seg.resize(file_lastsize);
    wrapper->dataset()->readBlock
      (static_cast<int>(numseg-1), seg.data(), 0, seg.size());
    // This is the point of no return when it comes to preserving
    // the file unchanged in case an error happens. So, maybe not
    // delete that segment quite yet, and just logically delete it
    // instead? This is risky because e.g. xx_eof in the relayed
    // instance will need to know about the subterfuge.
    // See more caveats in SeismicStoreFile::xx_write().
    wrapper->dataset()->deleteBlock(std::to_string(numseg-1));
    // Force the cached DatasetInformation to get updated to
    // reflect the new state of the file, i.e. one segment less.
    wrapper->updateDataset(wrapper->dataset());

    // CHECK that the code didn't mess things up.
    // Should really be testing this->xx_eof() but because this
    // method is called from the constructor our own virtuals
    // will not work.
    if (this->_relay->xx_eof() + std::int64_t(seg.size()) != numbytes)
      throw OpenZGY::Errors::ZgyInternalError("Bug: eof was changed by reopen.");
    if (static_cast<std::int64_t>(wrapper->dataset()->getBlockNum()) !=numseg-1)
      throw OpenZGY::Errors::ZgyInternalError("Bug: sdapi didn't delete.");
  }
}

/**
 * Flush "this_segsize" pending bytes. Leave any residual data
 * in the open segment buffer.
 *
 * Thread safety: No.
 */
void
SeismicStoreFileDelayedWrite::_flush_part(std::int64_t this_segsize)
{
  std::vector<char>& seg = this->_open_segment;
  if (seg.size() < static_cast<std::uint64_t>(this_segsize) || seg.size() <= 0)
    throw OpenZGY::Errors::ZgyInternalError("Bad segment size in _flush_part");
  this->_relay->xx_write(seg.data(), this->_relay->xx_eof(), this_segsize, this->_usage_hint);
  seg.erase(seg.begin(), seg.begin() + this_segsize);
  if (this->_config->_debug_trace)
    this->_config->_debug_trace("flush", this_segsize, this_segsize, 1, this->xx_segments(true));
}

/**
 * Flush pending writes, but only if we have enough data to fill
 * one or more complete segments or if the file is being closed.
 * The last segment is allowed to be smaller than the others.
 *
 * Thread safety: No.
 */
void
SeismicStoreFileDelayedWrite::_flush(bool final_call)
{
  if (this->_config->_segsize > 0)
    while (this->_open_segment.size() >= static_cast<std::uint64_t>(this->_config->_segsize)) {
      _relay->_mylogger(2, "_flush " + std::to_string(this->_config->_segsize) + " bytes\n");
      this->_flush_part(this->_config->_segsize);
    }
  if (final_call && this->_open_segment.size() > 0) {
    _relay->_mylogger(2, "_final " + std::to_string(this->_open_segment.size()) + " bytes\n");
    this->_flush_part(this->_open_segment.size());
  }
  if (this->_open_segment.size() == 0)
    this->_usage_hint = UsageHint::Unknown;
}

/**
 * Expose this internal extension to the SDManager class for use in unit tests.
 * CAVEAT: The function prototype may be declared locally.
 */
OPENZGY_TEST_API void
hack_setAuthProvider(seismicdrive::SDManager *mgr,
                     const std::string& token)
{
  SDGenericDatasetWrapper::setAuthProvider(mgr, token, nullptr);
}

namespace {
  /**
   * \details: Thread safety: Yes. add_factory() is synchronized.
   */
  class Register
  {
  public:
    Register()
    {
      FileFactory::instance().add_factory(SeismicStoreFile::xx_make_instance);
      FileFactory::instance().add_factory(SeismicStoreFileDelayedWrite::xx_make_instance);
    }
  } dummy;
} // anonymous namespace for registration

} // namespace

/** \endcond */

#endif // Almost entire file excluded if !HAVE_SD

namespace InternalZGY {
/**
 * Return a functor suitable for use in ctx.sdtokencb().
 * New tokens are generated from a SDManager created with a string
 * token, which also supports e.g. "ServiceSauthV2::id::secret".
 * If those credentials support refreshing then the callback does too.
 *
 * Using one manager to provide credentials from another might not
 * be technically legal, because the callback should not make
 * calls back to SDAPI. In this case it will probably work,
 * but don't use this trick in production without careful analysis.
 *
 * getIDToken() is thread safe so I shouldn't need additional
 * locking here.
 */
OPENZGY_TEST_API std::function<std::string()>
hack_getTokenCallback(const std::string sdurl,
                      const std::string sdkey,
                      const std::string token)
{
#ifdef HAVE_SD
  std::shared_ptr<seismicdrive::SDManager> mgr
    (new seismicdrive::SDManager(sdurl, sdkey));
  SDGenericDatasetWrapper::setAuthProvider(mgr.get(), token, nullptr);
  return [mgr]() {
           std::string token = mgr->getIDToken();
           static std::string bearer("Bearer ");
           return token.rfind(bearer) ? token : token.substr(bearer.size());
         };
#else
  return std::function<std::string()>();
#endif
}
} // namespace

