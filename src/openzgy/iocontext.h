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
 * \file iocontext.h
 * \brief Backend specific context.
 *
 * Class IOContext and derivatives are used to hold backend specific
 * information such as authorization tokens.
 */

#pragma once

#include "declspec.h"
#include "exception.h"

#include <cstdint>
#include <string>
#include <vector>
#include <functional>
#include <memory>

namespace Test
{
  class TestIOContext;
}
namespace InternalZGY
{
  class SeismicStoreFile;
  class SeismicStoreFileDelayedWrite;
}
namespace OpenZGY {
  class IZgyUtils;
}

namespace OpenZGY {
#if 0
}
#endif

/**
 * \brief Base class for backend specific context.
 *
 * \details Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_API IOContext
{
public:
  virtual ~IOContext();
  /** Display the context in a human readable format for debugging. */
  virtual std::string toString() const = 0;
  /** Create an identical deep copy */
  virtual std::shared_ptr<IOContext> clone() const = 0;
};

/** \cond SSTORE */

/**
 * \brief Credentials and configuration for Seismic Store.
 *
 * Define an iocontext for seismic store, doing consistency checks
 * and applying fallbacks from environment variables and hard coded
 * defaults.
 *
 * TODO-Low: Still undecided whether I should allow SeismicStoreFile
 * to use this class directly or whether I should map the contents to
 * an internal SDConfig class.
 *
 * TODO-Low: Move this class to a separate extensions/seismic_store.h
 * to be included by applications if and only if they need that access.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are only meant to be modified when created or
 * copied or assigned prior to being made available to others.
 */
class OPENZGY_API SeismicStoreIOContext : public IOContext
{
  // Needed because there are no public accessors. Because I am lazy.
  friend class Test::TestIOContext;
  friend class InternalZGY::SeismicStoreFile;
  friend class InternalZGY::SeismicStoreFileDelayedWrite;

public:
  typedef std::function<std::string()> tokencb_t;
  typedef std::function<void(const std::string&, std::int64_t, std::int64_t, std::int64_t, const std::vector<std::int64_t>&)> debugtrace_t;
  typedef std::function<bool(int, const std::string&)> logger_t;
  std::string toString() const override;
  std::shared_ptr<IOContext> clone() const override;
  SeismicStoreIOContext();

private:
  std::string  _sdurl;
  std::string  _sdapikey;
  std::string  _sdtoken;
  std::int64_t _maxsize;
  std::int64_t _maxhole;
  std::int64_t _aligned;
  std::int64_t _segsize;
  std::int64_t _segsplit;
  std::int64_t _iothreads;
  std::int64_t _cputhreads;
  std::int64_t _writethreads;
  std::string  _legaltag;
  std::string  _writeid;
  std::string  _seismicmeta;
  tokencb_t    _sdtokencb2;
  debugtrace_t _debug_trace;
  logger_t     _logger;
  bool         _set_ro_after_write;
  bool         _force_ro_before_read;
  bool         _force_rw_before_write;
  std::int32_t _retry_count;

private: // really private. Keep off.
  std::int64_t _real_segsize;

public:
  /**
   * Where to contact the seismic store service.
   * Defaults to $OPENZGY_SDURL. There is no hard coded fallback
   * in case that variable isn't found either. This is to mitigate
   * the risk of code stopping to work in the PROD environment.
   */
  SeismicStoreIOContext& sdurl(const std::string& value)
  {
    this->_sdurl = value;
    return *this;
  }

  /**
   * Authorization for application to access the seismic store API.
   * Defaults to $OPENZGY_SDAPIKEY. There is no hard coded fallback
   * in case that variable isn't found either. This is to mitigate
   * the risk of code stopping to work in the PROD environment.
   */
  SeismicStoreIOContext& sdapikey(const std::string& value)
  {
    this->_sdapikey = value;
    return *this;
  }

  /**
   * Provide the SAuth token used to validate requests to seismic store.
   * The token is associated with the open file and cannot be changed.
   * If neither sdtoken() nor sdtokenCb() are called then the environment
   * variable "OPENZGY_TOKEN" is tried.
   */
  SeismicStoreIOContext& sdtoken(const std::string& value) {
    this->_sdtoken = value;
    this->_sdtokencb2 = tokencb_t();
    return *this;
  }

  /**
   * Deprecated overload. Token type will be ignored and should be omitted.
   */
  SeismicStoreIOContext& sdtoken(const std::string& value, const std::string&) {
    return sdtoken(value);
  }

  /**
   * Register a callback that will be invoked when a fresh access token
   * is needed for seismic store. This is an alternative to sdtoken().
   * You don't need both.
   *
   * The callback should always return a plain token. Fancy alternatives
   * such as impersonation tokens, client credentials grant, etc. should
   * all refresh automatically. They don't need a callback for that purpose.
   *
   * The callback implementation should not invoke any OpenZGY methods
   * except for any calls expressly documented as safe and deadlock-free.
   * The callback implementation must also be thread safe. Set a lock
   * if needed. Finally, the callback needs to be legal to call indefinitely.
   * It is ok for the callback to throw if it can no longer obtain a token.
   */
  SeismicStoreIOContext& sdtokencb(const tokencb_t& value) {
    this->_sdtokencb2 = value;
    this->_sdtoken = std::string();
    return *this;
  }

  /**
   * Deprecated overload. Token type will be ignored and should be omitted.
   */
  SeismicStoreIOContext& sdtokencb(const tokencb_t& value, const std::string&) {
    return sdtokencb(value);
  }

  /**
   * Convenience function to share credentials with an open ZgyUtils instance.
   * This is required when dealing with SAuth impersonation tokens. And other
   * token types than can only have a single owner.
   *
   * The ZgyUtils instance will be kept alive as long as needed by the ZgyReader
   * or whatever this IOContext is given to.
   *
   * Calling context.credentialsFrom(utils) is equivalent to calling
   *    context.sdtokencb([utils](){return utils->idToken();});
   *
   * To instead share credentials with an SDManager created outside OpenZGY:
   *    context.sdtokencb([mgr](){return mgr->getIDToken()();});
   *
   * There is no convenience function for the SDManager case because that
   * would expose SDAPI in the public iterface. Just copy/paste the lambda
   * instead. Exposing IZgyUtils here is bad enough.
   *
   * the "mgr" or "utils" argument captured by the lambda should preferably
   * be a smart smart pointer. Using a raw pointer will also work. But the
   * application then takes FULL RESPONSIBILITY for having the manager or
   * sdutils instance stay in scope long enough. E.g. by leaking it.
   */
  SeismicStoreIOContext&
  credentialsFrom(std::shared_ptr<IZgyUtils> utils);
  // Moved to iocontext.cpp due to the dependency on api,h
  //{
  //  this->sdtokencb([utils](){return utils->idtoken();});
  //  return *this;
  //}

  /**
   * Applies to read() and finalize().
   *
   * Maximum size of consolidated requests, in MB. Must be
   * between 0 and 1024. Zero is taken to mean do not consolidate.
   *
   * Tell the reader to try to consolidate neighboring bricks
   * when reading from seismic store. This is usually possible
   * when the application requests full traces or at least traces
   * traces longer then 64 samples. Setting maxsize limits this
   * consolidation to the specified size. The assumption is that
   * for really large blocks the per-block overhead becomes
   * insignificant compared to the transfer time.
   *
   * Consolidating requests has higher priority than using
   * multiple threads. So, capping maxsize might allow more
   * data to be read in parallel.
   *
   * Note that currently the spitting isn't really smart. With a
   * 64 MB limit and 65 contiguous 1 MB buffers it might end up
   * reading 64+1 MB instead of e.g. 32+33 MB.
   *
   * Note that the low level reader should not assume that
   * requests are capped at this size. They might be larger
   * e.g. when reading the header information.
   *
   * Defaults to $OPENZGY_MAXSIZE_MB if not specified, or 64 MB.
   */
  SeismicStoreIOContext& maxsize(int value)
  {
    if (value < 0 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("maxsize must be between 0 and 1024 MB");
    this->_maxsize = value * (std::int64_t)(1024*1024);
    return *this;
  }

  /**
   * Applies to read() and finalize().
   *
   * Maximum memory to waste, in MB. Must be between 0 and 1024.
   *
   * This applies when consolidate neighboring bricks when
   * reading from seismic store. Setting maxhole > 0 tells the
   * reader that it is ok to also consolidate requests that are
   * almost neighbors, with a gap up to and including maxhole.
   * The data read from the gap will be discarded unless picked
   * up by some (not yet implemented) cache.
   *
   * For cloud access with high bandwidth (cloud-to-cloud) this
   * should be at least 2 MB because smaller blocks will take
   * just as long to read. For low bandwidth cloud access
   * (cloud-to-on-prem) it should be less. If a fancy cache
   * is implemented it should be more. For accessing on-prem
   * ZGY files it probably makes no difference.
   * Defaults to $OPENZGY_MAXHOLE_MB if not specified, or 2 MB.
   */
  SeismicStoreIOContext& maxhole(int value)
  {
    if (value < 0 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("maxhole must be between 0 and 1024 MB");
    this->_maxhole = value * (std::int64_t)(1024*1024);
    return *this;
  }

  /**
   * Applies to read() and finalize().
   *
   * File read alignment, in MB. Must be between 0 and 1024.
   *
   * This is similar to the maxhole parameter. If set, starting
   * and ending offsets are extended so they both align to the
   * specified value. Set this parameter if the lower levels
   * implement a cache with a fixed blocksize and when there is
   * an assumpton that most reads will be aligned anyway.
   * TODO-Worry: Handling reads past EOF may become a challenge
   * for the implementation.
   * Defaults to $OPENZGY_ALIGNED_MB if not specified, or zero.
   */
  SeismicStoreIOContext& aligned(int value)
  {
    if (value < 0 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("aligned must be between 0 and 1024 MB");
    this->_aligned = value * (std::int64_t)(1024*1024);
    return *this;
  }

  /**
   * Applies to write() and finalize().
   *
   * Segment size used when writing, in MB. Must be between 1 and
   * 16*1024 (i.e. 16 GB). Defaults to $OPENZGY_SEGSIZE_MB if not
   * specified, or 256 MB. The default should work fine in almost all
   * cases. The write buffer needs to allocate segsize * segsplit
   * bytes, so make sure to take segsplit into account when deciding
   * on the segment size. With both of those left at the default value
   * the buffer will need 256 MB * 8 = 2 GB.
   */
  SeismicStoreIOContext& segsize(int value)
  {
    if (value < 0 || value > 16*1024)
      throw OpenZGY::Errors::ZgyUserError("segsize must be between 0 and 16*1024 MB");
    this->_real_segsize = value * (std::int64_t)(1024*1024);
    // Internally, _segsize is the buffer size to be used, not the
    // smaller chunk size seen by the cloud backend.
    // _segsize is derived, recomputed when either factor changes.
    this->_segsize = this->_real_segsize * (this->_segsplit<1 ? 1 : this->_segsplit);
    return *this;
  }

  /**
   * For unit tests only. Allows setting segment size to any value.
   */
  SeismicStoreIOContext& segsizebytes(std::int64_t nbytes)
  {
    this->_real_segsize = nbytes;
    this->_segsize = this->_real_segsize * (this->_segsplit<1 ? 1 : this->_segsplit);
    return *this;
  }

  /**
   * Applies to write() and finalize().
   *
   * Maximum number of threads to be used when writing data to the cloud.
   * Set to 1 if you don't want multithreaded uploads. Default is 8.
   *
   * Multi-threading is achieved by using an internal buffer size of
   * segsize*segsplit, then later splitting that up into separate
   * SDAPI objects than can be uploaded in parallel. So, beware that
   * if you increase segsplit you might need to decrease segsize to
   * avoid running out of memory.
   *
   * Using multiple threads for uploads can often improve throughput.
   * But if the limiting factor is the bandwidth between the client
   * and the cloud server then multiple treads won't help.
   *
   * Ignored (will always be 1) if writethreads > 1.
   */
  SeismicStoreIOContext& segsplit(int value)
  {
    if (value < 1 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("segsplit must be between 0 and 1024");
    // Temporary, for testing writethreads. No application code sets
    // this yet, so the only way is to set OPENZGY_NUMTHREADS_WRITE.
    // Once that has been set, any attempts to enable the old
    // parallel-write mechanism will be ignored.
    if (_writethreads > 1)
      value = 1;
    this->_segsplit = value;
    this->_segsize = this->_real_segsize * (this->_segsplit<1 ? 1 : this->_segsplit);
    return *this;
  }

  /**
   * Applies to read() and finalize() and misaligned write().
   *
   * Typically iothreads and cputhreads will be 1 in a file opened for
   * read because most applications handle multi threading of reads
   * themselves. Typically they are set higher in a file opened for
   * write because then the parameters are used in finalize(), to read
   * back the data that is to be decimated, and finalize() does *not*
   * do its own multi threading for reads. The parameters are also
   * used by write() and writeconst() if the old contents are needed
   * due to misaligned writes or tracking of changes.
   *
   * Use up to this many parallel requests to seismic store in order
   * to speed up processing. Set between 1 and 1024. This applies to
   * individual reads in the main API. So the reads must be for a
   * large area (i.e. covering many bricks) for the setting to be
   * of any use. Set to $OPENZGY_NUMTHREADS_IO if not set here, and 1
   * (i.e. no threading) if the environment setting is also missing.
   *
   * Whether it is useful to set the variable depends on the
   * application. Apps such as Petrel/BASE generally do their own
   * multi threading, issuing multiple read requests to the high level
   * API in parallel. In that case it might not be useful to also
   * parallelize individual requests.
   * It might even be a very bad idea if the total number of pending
   * requests hits some service- or network limit.
   *
   * If the application normally asks for small regions or moderately
   * sized regions that are contiguous on the file then this setting
   * has no effect. Contrariwise, requesting an entire horizontal
   * slice can use a lot of threads if it is allowed to. And if the
   * application is single threaded it really should allow this.
   */
  SeismicStoreIOContext& iothreads(int value)
  {
    if (value < 1 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("iothreads must be between 1 and 1024");
    this->_iothreads = value;
    return *this;
  }

  /**
   * Applies to read() and finalize() and misaligned write().
   *
   * Use up to this many parallel requests for cpu-intensive
   * operations such as decompression. Set between 1 and 1024.
   * Set to $OPENZGY_NUMTHREADS_CPU if not set here, and 1
   * (i.e. no threading) if the environment setting is also missing.
   *
   * As with iothreads it might not be much use for Petrel but
   * it is worth trying. There is less risk of hitting hard limits.
   */
  SeismicStoreIOContext& cputhreads(int value)
  {
    if (value < 1 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("cputhreads must be between 1 and 1024");
    this->_cputhreads = value;
    return *this;
  }

  /**
   * Applies to write() and finalize(). Not implemented / enabled yet.
   *
   * Alternative method of doing multi-threaded write by making SDAPI
   * write operations asynchronous.
   *
   * Setting segsplit != 1 will disable writethreads and vice versa.
   * Note that writethreads is not implemented yet, so setting
   * writethreads > 1 will actually disable multi-threaded writes.
   */
  SeismicStoreIOContext& writethreads(int value)
  {
    if (value < 1 || value > 1024)
      throw OpenZGY::Errors::ZgyUserError("writethreads must be between 0 and 1024");
    this->_writethreads = value;
    if (value != 1)
      segsplit(1);
    return *this;
  }

  /**
   * The legaltag stored in the file. Used only on create.
   */
  SeismicStoreIOContext& legaltag(const std::string& value)
  {
    this->_legaltag = value;
    return *this;
  }

  /**
   * If set, re-use this lock instead of creating a new one.
   * Works both for read and write locks; the name reflects
   * what SDAPI calls it.
   */
  SeismicStoreIOContext& writeid(const std::string& value)
  {
    this->_writeid = value;
    return *this;
  }

  /**
   * a dictionary of additional information to be associated
   * with this dataset in the data ecosystem. Currently used
   * only on create, although SDAPI allows this to be set on
   * an existing file by calling {get,set}SeismicMeta().
   * This setting cannot be set from the environment.
   */
  SeismicStoreIOContext& seismicmeta(const std::string& value)
  {
    this->_seismicmeta = value;
    return *this;
  }

  /**
   * For debugging and unit tests only.
   * Callback to be invoked immediately before a read or write
   * is passed on to seismic store. Typically used to verify
   * that consolidating bricks works as expected. Can only be
   * set programmatically. Not by an environment variable.
   *
   * The callback implementation should not invoke any OpenZGY
   * methods except for any calls expressly documented as safe
   * for this purpose and deadlock-free.
   */
  SeismicStoreIOContext& debug_trace(const debugtrace_t& value)
  {
    this->_debug_trace = value;
    return *this;
  }

  /**
   * For debugging. This callback will be invoked whenever
   * something needs to be written to the console or to a file.
   * See LoggerBase::standardCallback() for an example.
   *
   * The default is to use standardCallback() configured
   * according to OPENZGY_VERBOSE, OPENZGY_VERBOSE_LOGFILE,
   * and OPENZGY_VERBOSE_DETAILS.
   *
   * Note that it might be a good idea to use a std::ref here,
   * that could avoid some surprises.
   */
  SeismicStoreIOContext& logger(const logger_t& logger)
  {
    this->_logger = logger;
    return *this;
  }

  /**
   * Internal access to read out the current logger.
   * The logger is applied to higher level parts
   * of the system and needs this slight breakage
   * of encapsulation.
   */
  logger_t getLogger() const
  {
    return this->_logger;
  }

  /**
   * Set the ZGY file to read-only when done writing it. Has no effect
   * on files opened for read. Defaults to on. Most applications will
   * want to turn in this on because most applications do not expect
   * to update ZGY files in place.
   *
   * \internal
   * It is not clear whether this belongs in IOContext next to the two
   * "force" options or if it should be a new finalize tag
   * BuildFullAndSetReadOnly.
   *
   * IOContext
   *  - Not available for files.
   *  - Cannot default to "on" because it depends on the finalise mode
   *    what makes sense.
   *
   * FinalizeAction:
   *  - The flag cannot be set inside finalize, it needs to
   *    be deferred until the last flush. So it doesn't really
   *    belong here either.
   *  - It isn't just BuildFull which needs an ...AndSetReadOnly variant
   *    although the others are less likely to be used.
   */
  SeismicStoreIOContext& setRoAfterWrite(bool value)
  {
    this->_set_ro_after_write = value;
    return *this;
  }

  /**
   * Sneak past the mandatory locking in SDAPI by forcing the
   * read-only flag to true on the ZGY file, if needed, on each open
   * for read. This allows for less overhead, more caching, and use of
   * the altUrl feature. This option is useful if the file is
   * logically immutable but was not flagged as such. E.g. the creator
   * forgot to call setRoCloseWrite(true), or the ZGY file was not
   * created by OpenZGY. The option has no effect on files opened for
   * create or update. Caveat: Allowing a read-only open to have a
   * permanent effect of the file being opened is not ideal.
   */
  SeismicStoreIOContext& forceRoBeforeRead(bool value)
  {
    this->_force_ro_before_read = value;
    return *this;
  }

  /**
   * Dangerous option. Sneak past the mandatory locking in SDAPI by
   * forcing the read-only flag to false on the ZGY file, if needed,
   * that is about to be opened for update. The application must know
   * that the file is not open for read by somebody else. There is
   * also a risk that data might exists in cache even for a closed
   * file. The application assumes all responsibility.
   */
  SeismicStoreIOContext& forceRwBeforeWrite(bool value)
  {
    this->_force_rw_before_write = value;
    return *this;
  }

  /**
   * Maximum number of retries that SDAPI will use before concluding
   * that a web service is unavailable and not just experiencing a
   * transient problem. There is an exponential backoff between
   * attempts. 1/2 second for the first one, then doubling for each
   * attempt but never more than 32 seconds. Leaving the parameter
   * unset or set to -1 will use the defaults hard coded in SDAPI.
   */
  SeismicStoreIOContext& retryCount(int value)
  {
    this->_retry_count = value;
    return *this;
  }
};

/** \endcond */

}

// namespace
