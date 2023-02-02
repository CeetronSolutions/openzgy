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

#include "file_parallelizer.h"
#include "mtguard.h"
#include "exception.h"
#include <omp.h>
#include <string.h>
#include <iostream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \brief Help parallelize the decompression and copy-out steps.
 *
 * \detailed
 * - Intercept calls to xx_readv.
 *
 * - Send a single read request list down to the layer below.
 *   At the level where the parallelized is injected, each
 *   individal request will almost always be for a single brick.
 *
 * - Wait for all the data to be delvered. In the seismic server
 *   case the accessor will in any case defer sending results until
 *   it has all the requested data. Changing that will not happen
 *   anytime soon
 *
 * - Use an OpemMP loop to return each brick in a potentially
 *   different thread.
 *
 * The net effect is similar to using parallel loops at the
 * end of SeismicStoreFile::xx_readv() and also in in the
 * ConsolidateRequests::_consolidated_delivery() method.
 * In that case no additional buffer copy would have been needed.
 * The first place is for when bricks were not consolidated and the
 * second for when that did happen. The problem is that this may lead
 * to nested parallelization.
 *
 * Hence the code in this file. Which is in any case better
 * due to being more modular. The extra buffer copies are bad
 * though.
 */
FileParallelizer::FileParallelizer(std::shared_ptr<IFileADT> relay, std::int64_t cputhreads)
  : FileRelay(relay)
  , _cputhreads(cputhreads)
{
  //std::cerr << "Parallelizer has been created with "
  //          << _cputhreads << " threads\n";
}

FileParallelizer::~FileParallelizer()
{
  //std::cerr << "Parallelizer has been destructed\n";
}

void
FileParallelizer::xx_readv(
     const ReadList& requests,
     bool parallel_ok,  // usually true
     bool immutable_ok, // usually false
     bool transient_ok, // usually true
     UsageHint hint)
{
  const bool shortcut = false; // Set false only for debugging.
  const std::int64_t requestcount = requests.size();

  // Shortcut if no parallelizing possible.
  if (!parallel_ok || (shortcut && requests.size() <= 1)) {
    //std::cerr << "Parallelizer: Nothing to do\n";
    relay().xx_readv(requests, parallel_ok, immutable_ok, transient_ok, hint);
    return;
  }

  // Future: For each request, try fulfilling it entirely from the cache
  // and remove from the queue because that part would have much lower
  // overhead.

  // Negotiating immutable_ok and transient_ok:
  // With delivery using smart pointers, SeismicStoreFile::xx_readv()
  // has little or no additional cost for a mutable and not transient
  // buffer except for some corner cases. So pass down the immutable_ok
  // from the caller and pass the transient_ok=false needed here.
  // A future caching module might make this trickier.

  // The new request list sends the data to our buffers instead of to caller.
  ReadList newrequests(requests);
  typedef std::pair<std::shared_ptr<const void>, std::int64_t> delivery_t;
  std::vector<delivery_t> buffers(requestcount);
  for (std::int64_t ii = 0; ii < requestcount; ++ii) {
    delivery_t* destination = &buffers[ii];
    newrequests[ii].delivery =
      [destination](ReadRequest::data_t data, std::int64_t len) {
        //std::cerr << "+";
        destination->first = data;
        destination->second = len;
      };
  }

  // The actual read. It will not return until all is processed.
  // We keep a reference to the delivered data. Caller must be told.
  //std::cerr << "Parallelizer: n=" << requestcount << " ";
  relay().xx_readv(newrequests, true, immutable_ok, false, hint);

  // Deliver the buffers that we cached to the caller.
  const std::int64_t threadcount = std::min(std::min(requestcount, (std::int64_t)omp_get_max_threads()), _cputhreads);
  MTGuard guard("paralellizer", (int)threadcount);
#pragma omp parallel for num_threads((int)threadcount)
  for (std::int64_t ii = 0; ii < requestcount; ++ii) {
    guard.run([&](){
      //std::cerr << "0123456789"[omp_get_thread_num() % 10];
      _deliver(requests[ii].delivery, buffers[ii].first, 0, buffers[ii].second, transient_ok);
    });
  }
  guard.finished();
  //std::cerr << "$\n";
}

/**
 * Inject a parallelizer module.
 */
std::shared_ptr<IFileADT>
FileParallelizer::inject(std::shared_ptr<IFileADT> file, std::int64_t cputhreads)
{
  if (cputhreads > 1)
    file = std::shared_ptr<IFileADT>(new FileParallelizer(file, cputhreads));
  return file;
}

} // namespace
