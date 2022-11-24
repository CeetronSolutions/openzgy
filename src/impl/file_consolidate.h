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

#include "file.h"
#include <ostream>

namespace InternalZGY {
#if 0
}
#endif

/**
 * Static class for scoping ConsolidateRequests::consolidate()
 * and a number of private methods it needs.
 *
 * Thread safety: Safe because this is a static class with no data.
 */
class ConsolidateRequests
{
  ConsolidateRequests() = delete;
  ConsolidateRequests(const ConsolidateRequests&) = delete;
  ConsolidateRequests& operator=(const ConsolidateRequests&) = delete;

public:

  static ReadList consolidate(
      const ReadList& requests,
      std::int64_t max_hole = 2*1024*1024,
      std::int64_t max_size = 64*1024*1024,
      std::int64_t force_align = 0,
      bool consolidate_overlaps = false,
      std::int64_t eof = 0);

  static void _print_requests(
      const ReadDoubleList& all_requests,
      const std::string& name,
      std::ostream& os);

private:

  static std::pair<std::int64_t, std::int64_t> _groupsize(
      const ReadList& g,
      std::int64_t force_align,
      std::int64_t eof);

  static ReadDoubleList _split_requests(
      const ReadList& requests,
      std::int64_t max_hole,
      std::int64_t max_size,
      std::int64_t force_align,
      bool consolidate_overlaps,
      std::int64_t eof);

  static ReadList _join_requests(
      const ReadDoubleList& all_requests,
      std::int64_t force_align,
      std::int64_t eof);

  static ReadRequest::delivery_t _consolidated_delivery(
      const ReadList& group,
      std::int64_t begin);

  /*
  RawRequestList _split_by_segment(
      const ReadList& requests);

  void do_readv(const RawRequestList& in,
      void *buffer,
      int threadcount);
  */
};
} // namespace
