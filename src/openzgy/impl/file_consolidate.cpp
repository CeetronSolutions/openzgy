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

#include "file_consolidate.h"
#include "fancy_timers.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file file_consolidate.cpp
 * \brief Called from file_sd.cpp to consolidate brics for cloud access.
 *
 * \details
 * Overview in pseudo-code of what this code needs to do.
 *
 * COMPLEX TYPES:
 *
 *    class ReadRequest {offset, size, delivery_functor(buffer, size)}
 *    class RawRequest  {seg_number, offset_in_seg, size_in_seg, outpos}
 *
 * FUNCTION CALL HIERARCHY:
 *
 *    xx_readv
 *       _consolidate_requests: ReadRequest[] => ReadRequest[]
 *          _split_requests: ReadRequest[] => ReadRequest[][]
 *             _groupsize: ReadRequest[] => integer
 *          _join_requests: ReadRequest[][] => ReadRequest[]
 *             _groupsize: ReadRequest[] => integer
 *             class ConsolidatedDelivery
 *          _print_requests: ReadRequest[][] => void
 *       _split_by_segment: ReadRequest[] => RawRequest[]
 *       do_readv: RawRequest[], buffer*, threadcount => void
 *          SDGenericDataset::readBlock: segno, buffer*, offset, buffersize
 *
 *    xx_read
 *      _split_by_segment: ReadRequest[] => RawRequest[]
 *      _cached_read: segno, offset, buffer* => void
 *         SDGenericDataset::readBlock: segno, buffer*, offset, buffersize
 *
 * REMARKS:
 *
 * SeismicStoreFile._split_by_segment is already implemented in
 * file_sd.cpp, it is needed even without enabling the consolidation
 * code.
 *
 * do_readv is already in C++ in the SdGlue Python extension and needs
 * to be extracted from there; will probably need some changes. If the
 * application isn't going to enable multi-threading of individual
 * requests then this method will not be needed. On the other hand,
 * multi threading can lead to some very impressive test results. The
 * do_readv method probably belongs in file_sd.cpp and not here.
 *
 * cached_read is part of the "Rema1000" kludge; I am not sure it
 * works properly and in any case we might need a proper caching
 * module. Exclude for now.
 *
 * Most of the rest is already implemented in Python so just need a
 * translation into C++.
 *
 * In the Python version the lexical scoping is similar to the call
 * hierarchy. E.g. _split_requests and _join_requests are local
 * functions inside _consolidate_requests. This changes in the C++
 * version.
 */

/**
 * Given a list of requests as passed to xx_readv, try to reduce the
 * number of requests by consolidating adjacent or nearly adjacent
 * reads. If successful this means we will be reading with larger
 * block sizes.
 *
 * Return a new list of requests that the caller may pass on to
 * xx_readv instead of the original.
 *
 * Remember that the callback functions specified with the original
 * requests need to be called with the exact data that they expected.
 * This means that in the consolidated list the callback functions
 * need to be wrappers.
 *
 * Set consolidate_overlaps to true if you expect some of the
 * individual requests to overlap, and you are ok with risking some
 * corner cases. For example, if you request a mutable buffer then the
 * overlapping area will be delivered to more than one recipient and
 * the buffer may or may not be shared between the two. The default is
 * False which causes the code to not attempt consolidation of these.
 * Less efficient but also less surprises. In practice there should
 * never be any overlap anyway.
 *
 *     \param max_hole    See class SeismicStoreIOContext for a description.
 *     \param max_size    See class SeismicStoreIOContext for a description.
 *     \param force_align See class SeismicStoreIOContext for a description.
 *     \param consolidate_overlaps See above for a description.
 */
ReadList
ConsolidateRequests::consolidate(
      const ReadList& requests,
      std::int64_t max_hole,
      std::int64_t max_size,
      std::int64_t force_align,
      bool consolidate_overlaps,
      std::int64_t eof)
{
  ReadDoubleList all_requests =
    _split_requests(requests, max_hole, max_size, force_align, consolidate_overlaps, eof);
  ReadList new_requests = _join_requests(all_requests, force_align, eof);
  if (false && requests.size() != new_requests.size()) {
    std::cout << "Consolidated " << requests.size()
              << " into " << new_requests.size() << "\n";
    _print_requests(ReadDoubleList{requests}, "Requests:", std::cout);
    _print_requests(ReadDoubleList{new_requests}, "Consolidated:", std::cout);
  }
  const std::int64_t old_size =
    std::accumulate(requests.begin(), requests.end(), std::int64_t(0),
                    [](std::int64_t a, const ReadRequest& b) {
                      return a + b.size;
                    });
  const std::int64_t new_size =
    std::accumulate(new_requests.begin(), new_requests.end(), std::int64_t(0),
                    [](std::int64_t a, const ReadRequest& b) {
                      return a + b.size;
                    });
  if (new_size < old_size && !consolidate_overlaps)
    throw std::runtime_error("Assert failed in ConsolidateRequests::consolidate");
  return new_requests;
}

/**
 * Given a list of (offset, size, functor) return offset and size for
 * the entire group. The offset is the linear offset from the start of
 * the file; it has not yet been converted to segment and local
 * offset. The returned value includes any padding for force_align.
 *
 * TODO-High the padding is WRONG, because the alignment should be
 * done per segment. We may end up crossing segment boundaries
 * needlessly. And/or going past EOF. Going past EOF is critical
 * because in the subsequent call to _split_by_segment() we will
 * end up trying to actually read that part.
 *
 * Crossing segment boundaries is less of a problem.
 *
 *   - It will not happen if the headers are aligned at least to
 *     force_align, which is typically the cache bricksize.
 *
 *   - It will not happen if the file was uploaded with sdutil.
 *     In that case there will be just one segment.
 *
 *   - It is (alomst) not an issue if a proper cache is attached.
 *
 *   - A naive cache can align to 256 KB, this virtually guarantees
 *     the header area will be sufficiently aligned if the file
 *     was initially stored on the cloud.
 *
 *   - In other cases there will be a small performance penalty but
 *     only when reading close to a segment boundary or when reading
 *     the headers. Opening a file may see a noticeable slowdown
 *     but not I think anything dramatic.
 */
std::pair<std::int64_t, std::int64_t> // returns offset, size
ConsolidateRequests::_groupsize(
      const ReadList& g,
      std::int64_t force_align,
      std::int64_t eof)
{
  if (g.empty())
    return std::make_pair(0, 0);
  std::int64_t beg = g.front().offset;
  std::int64_t end = g.front().offset + g.front().size;
  for (const ReadRequest& e : g) {
    beg = std::min(beg, e.offset);
    end = std::max(end, e.offset + e.size);
  }
  if (beg != g.front().offset)
    throw std::runtime_error("Assert failed, input ReadList not sorted");
  // This might fail if requests overlap
  //if (end != g.back().offset + g.back().size)
  //  throw std::runtime_error("Assert failed, input ReadList looks odd");
  if (force_align != 0) {
    beg = (beg / force_align) * force_align;
    end = ((end + force_align - 1) / force_align) * force_align;
    // This check only needed if padding was added.
    if (eof != 0)
      end = std::min(end, eof);
  }
  return std::make_pair(beg, end - beg);
}

/**
 * Make a list of lists, grouping requests that should be read in a
 * single operation. Operates on linear addresses, so if any of the
 * input requests crossed a segment boundary then this will also be
 * the case for the output.
 *
 * Minor "feature": The check for small holes (to be read and discarded)
 * is done before padding is applied. This means that even of padding
 * caused two chunks to become adjacent, they might not get consolidated.
 * I believe this is an academic issue because max_hole will typically be
 * at least as large as the alignment.
 */
ReadDoubleList
ConsolidateRequests::_split_requests(
      const ReadList& requests,
      std::int64_t max_hole,
      std::int64_t max_size,
      std::int64_t force_align,
      bool consolidate_overlaps,
      std::int64_t eof)
{
  ReadDoubleList all_requests;
  std::int64_t prev_end = 0;
  ReadList sorted_requests = requests;
  std::sort(sorted_requests.begin(), sorted_requests.end(),
            [](const ReadRequest& a, const ReadRequest& b) {
              return (a.offset < b.offset ||
                      (a.offset == b.offset && a.size < b.size));
            });
  for (const ReadRequest& request : sorted_requests) {
    std::int64_t hole = request.offset - prev_end;
    if (!all_requests.empty() && hole <= max_hole &&
        (consolidate_overlaps || hole >= 0))
    {
      // Consolidate!
      all_requests.back().push_back(request);
    }
    else
    {
      all_requests.push_back(ReadList{request});
    }
    if (max_size != 0 && all_requests.back().size() > 1 &&
        _groupsize(all_requests.back(), force_align, eof).second > max_size)
    {
      // Oops, that last request became too large. Don't consolidate after all.
      all_requests.back().pop_back();
      all_requests.push_back(ReadList{request});
    }
    prev_end = request.offset + request.size;
  }
  return all_requests;
}

/**
 * Create the final result containing one entry for each consolidated group.
 */
ReadList
ConsolidateRequests::_join_requests(
      const ReadDoubleList& all_requests,
      std::int64_t force_align,
      std::int64_t eof)
{
  ReadList new_requests;
  for (const ReadList& group : all_requests) {
    // Output of groupsize tells us what to read in order to cover
    // all the requests in this group. We can then create a single
    // ReadRequest for the group.
    std::pair<std::int64_t,std::int64_t> info =
      _groupsize(group, force_align, eof);
    auto deliver = _consolidated_delivery(group, info.first);
    new_requests.push_back(ReadRequest(info.first, info.second, deliver));
  }
  return new_requests;
}

/**
 * For debugging only, print a list of list of requests.
 */
void
ConsolidateRequests::_print_requests(
      const ReadDoubleList& all_requests,
      const std::string& name,
      std::ostream& outstream)
{
  static auto printit = [&outstream]
    (const char *name, std::int64_t offset, std::int64_t size) {
      outstream << std::dec
                << "    " << name
                << " offset " << std::setw(8) << SummaryPrintingTimerEx::niceSize(offset)
                << " end " << std::setw(8) << SummaryPrintingTimerEx::niceSize(offset+size)
                << " size " << std::setw(6) << SummaryPrintingTimerEx::niceSize(size) << std::dec << "\n";
  };
  if (all_requests.empty() ||
      (all_requests.size() == 1 && all_requests.front().empty())) {
    outstream << "    (empty)" << std::endl;
    return;
  }
  outstream << name << "\n";
  for (const ReadList& group : all_requests) {
    if (all_requests.size() > 1)
      outstream << "  Group:\n";
    std::int64_t prev_offset = -1, prev_size = 0;
    for (const ReadRequest& rr : group) {
      if (prev_offset != -1) {
        std::int64_t skip_offset = prev_offset + prev_size;
        std::int64_t skip_size = rr.offset - (prev_offset + prev_size);
        if (skip_size != 0)
          printit("skip", skip_offset, skip_size);
      }
      printit("read", rr.offset, rr.size);
      prev_offset = rr.offset;
      prev_size = rr.size;
    }
  }
}

/**
 * Help distribute a single delivery from a consolidated read request
 * to all the requesters. Slice the data so each requester gets
 * exactly what they originally asked for.
 *
 * Note that if the original request had overlapping reads we might
 * want to force a buffer copy. Because we don't know whether the
 * recipient asked for a mutable buffer. It is tempting to disallow
 * overlapping reads completely.
 *
 * Caveat: Handling reads past EOF may be tricky. I need some specific
 * unit tests for that.
 */
ReadRequest::delivery_t
ConsolidateRequests::_consolidated_delivery(
      const ReadList& group,
      std::int64_t begin)
{
  // Note that "group" needs to be copied by value because the list
  // being passed in will go out of scope when consolidate() exits. I
  // might choose to use a smart pointer to avoid a deep copy every
  // time the functor is copied. Technically I could make this even
  // more performant by having other methods use smart pointers as
  // well. But the value of that is questionable.
  std::shared_ptr<const ReadList> group_ptr(new ReadList(group));
  return [group_ptr,begin](ReadRequest::data_t data, std::int64_t size) {
           for (const ReadRequest& rr : *group_ptr) {
             if (rr.delivery) {
               std::int64_t end = std::min(rr.offset + rr.size - begin, size);
               std::int64_t beg = std::min(rr.offset - begin, end);
               // Caller will check the transient_ok flag. we won't.
               FileADT::_deliver(rr.delivery, data, beg, end - beg, false);
             }
           }
         };
}

} // namespace
