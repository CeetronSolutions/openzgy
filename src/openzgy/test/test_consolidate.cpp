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
#include "../exception.h"
#include "../impl/file.h"
#include "../impl/fancy_timers.h"
#include "../impl/file_consolidate.h"

#include <iostream>
#include <sstream>
#include <memory>


using namespace OpenZGY;
using namespace OpenZGY::Errors;
using namespace InternalZGY;

namespace {
#if 0
}
#endif

static int _nextseq = 0;

class MyReceiver
{
public:
  std::shared_ptr<int> firstval;
  std::shared_ptr<int> gotsize;
  std::shared_ptr<int> callcount;
  int seq;
public:
  MyReceiver() : firstval(new int), gotsize(new int), callcount(new int)
  {
    *firstval = 0;
    *gotsize = 0;
    *callcount = 0;
    seq = ++_nextseq;
    //std::cout << "MyReceiver #" << seq << " was constructed" << std::endl;
  }

  MyReceiver(const MyReceiver& other) : firstval(other.firstval), gotsize(other.gotsize), callcount(other.callcount)
  {
    seq = ++_nextseq;
    //std::cout << "MyReceiver #" << seq << " was cloned" << std::endl;
  }

  void operator()(ReadRequest::data_t data, std::int64_t size)
  {
    *callcount += 1;
    *gotsize = size;
    if (size >= 3*static_cast<std::int64_t>(sizeof(int))) {
      const int* const idata = reinterpret_cast<const int*>(data.get());
      *firstval = idata[0];
      if (false)
        std::cout << "Receiver # " << seq
                  << ": Deliver " << SummaryPrintingTimerEx::niceSize(size) << ": "
                  << idata[0] << ", " << idata[1] << ", " << idata[2] << ", ..."
                  << "\n";
    }
    else {
      //std::cout << std::hex << "Deliver " << size << " bytes.\n" << std::dec;
    }
  }

  bool check(int expect_firstval, int expect_gotsize) const
  {
    if (*callcount == 1 && *gotsize == expect_gotsize && *firstval == expect_firstval)
      return true;
    std::cout << "ERROR: Did not get expected delivery."
              << "\n  Expect "  << 1 << " delivery"
              << " with first " << expect_firstval
              << " and size "   << expect_gotsize
              << "\n  Actual "  << *callcount << " delivery"
              << " with first " << *firstval
              << " and size "   << *gotsize
              << std::endl;
    return false;
  }
};

/**
 * Trivial case: Empty list in, empty list out.
 */
static void
test_zero()
{
  ReadList list{};
  ReadList result = ConsolidateRequests::consolidate(list);
  TEST_CHECK(result.size() == 0);
}

/**
 * Trivial case: Single-block request in, equivalent request out.
 * The code might do a short cut and just return the input list,
 * or it can insert a (in this case) unneeded method to forward
 * the delivery.
 */
static void
test_one()
{
  MyReceiver receiver;
  ReadList list{ReadRequest(0x1000, 0x3000, receiver)};
  ReadList result = ConsolidateRequests::consolidate(list);
  TEST_CHECK(result.size() == 1);
  TEST_CHECK(result[0].offset == list[0].offset);
  TEST_CHECK(result[0].size == list[0].size);
  TEST_CHECK((bool)result[0].delivery);
  std::shared_ptr<int> data(new int[list[0].size], std::default_delete<int[]>());
  data.get()[0] = list[0].offset;
  data.get()[1] = list[0].size;
  for (std::int64_t ii=2; ii<list[0].size; ++ii)
    data.get()[ii] = 42;
  result[0].delivery(data, list[0].size * sizeof(int));
}

static void
run_test(const ReadList& list_in,
         const ReadList& expect,
         std::int64_t max_hole,
         std::int64_t max_size,
         std::int64_t align,
         bool do_overlaps,
         std::int64_t filesize)
{
  // Caller did not provide delivery callbacka, neither for the list to
  // process nor for the expected result. This probably only makes sense
  // in a unit test. Here I add my own callbacks for the input list.
  // The expect list is just used for offset/size; arguably this shouldn't
  // have been a ReadList at all.
  ReadList list(list_in);
  std::vector<MyReceiver> receivers(list.size());
  for (std::size_t ii = 0; ii < list.size(); ++ii)
    list[ii].delivery = receivers[ii];

  // Run the algorithm under test.
  ReadList result = ConsolidateRequests::consolidate
    (list, max_hole, max_size, align, do_overlaps, filesize);

  //ConsolidateRequests::_print_requests(ReadDoubleList{expect}, "Expect:", std::cout);

  // Compare expected and actual consolidated list.
  // This does not check whether data is delivered correctly.
  TEST_CHECK(result.size() == expect.size());
  for (std::size_t ii = 0; ii < result.size() && ii < expect.size(); ++ii)
    TEST_CHECK(result[ii].offset == expect[ii].offset &&
               result[ii].size   == expect[ii].size);

  // Pretend we read this data file a file. The raw data contains the
  // offset from the start of the buffer so we can check what was sent.
  std::shared_ptr<int> data(new int[filesize/sizeof(int)], std::default_delete<int[]>());
  for (std::size_t ii=0; ii<filesize/sizeof(int); ++ii)
    data.get()[ii] = sizeof(int)*static_cast<int>(ii);

  // Deliver my data according to the consolidated list.
  for (const ReadRequest& rr : result)
    FileADT::_deliver(rr.delivery, data, rr.offset, rr.size, false);

  // Check that all the data that was originally requested
  // was delivered to the correct place in a timely manner.
  // One delivery to a consolidated request might end up as
  // deliveries to more than one of the original requests.
  for (std::size_t ii = 0; ii < list.size(); ++ii)
    TEST_CHECK(receivers[ii].check(list[ii].offset, list[ii].size));
}

static void
test_simple()
{
  // Test data: 67 kb, may or may not be algned to 4 kb boundaries.
  // Note 67 is not divisible by 4, so eof is not aligned.
  // The test itself submits a list of 7 regions to access.
  //
  // Visualize the contents of the request as follows:
  //   - First line: Ruler to show 4 KB aligned blocks.
  //   - Second line: A digit if this data is being requested.
  //     The actual digit is the position of this request
  //     in the input; they don't need to be sorted.
  //   - Third line: Ruler for file size, one character per KB.
  //
  // [--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--
  // --------2222------------1111333---4444-----5555----------7777-6666-
  // 0000000000111111111122222222223333333333444444444455555555556666666
  // .123456789.123456789.123456789.123456789.123456789.123456789.123456
  //
  // The first 3 chunks have start aligned to 4 KB, the rest are not.
  // The third chunk has an odd size.
  // The 6th chunk triggers a corner case: If start and end are adjusted
  // to align with 4 KB then the end position will be past eof.
  //
  // This test: No reading across holes. Request part 1 and 3 get combined,
  // everything else is unchanged. When checking the result, note that the
  // consolidated list will always be sorted on offset.

  const std::int64_t KB = 1024;

  ReadList list{ReadRequest(24*KB, 4*KB, nullptr),
                ReadRequest( 8*KB, 4*KB, nullptr),
                ReadRequest(28*KB, 3*KB, nullptr),
                ReadRequest(34*KB, 4*KB, nullptr),
                ReadRequest(43*KB, 4*KB, nullptr),
                ReadRequest(62*KB, 4*KB, nullptr),
                ReadRequest(57*KB, 4*KB, nullptr)
  };

  ReadList expect{ReadRequest( 8*KB, 4*KB, nullptr),
                  ReadRequest(24*KB, 7*KB, nullptr),
                  ReadRequest(34*KB, 4*KB, nullptr),
                  ReadRequest(43*KB, 4*KB, nullptr),
                  ReadRequest(57*KB, 4*KB, nullptr),
                  ReadRequest(62*KB, 4*KB, nullptr)
  };

  run_test(list, expect,
           /*max_hole*/0, /*max_size*/0, /*align*/0,
           /*do_overlap*/false, /*eof*/67*KB);
}

static void
test_holes()
{
  // Same input data as test_simple but allow 8 KB holes.
  // The diagram below shows an asterisk for the data that will be
  // read and discarded.
  //
  // [--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--
  // --------2222------------1111333***4444*****5555----------7777*6666-
  // 0000000000111111111122222222223333333333444444444455555555556666666
  // .123456789.123456789.123456789.123456789.123456789.123456789.123456

  const std::int64_t KB = 1024;

  ReadList list{ReadRequest(24*KB, 4*KB, nullptr),
                ReadRequest( 8*KB, 4*KB, nullptr),
                ReadRequest(28*KB, 3*KB, nullptr),
                ReadRequest(34*KB, 4*KB, nullptr),
                ReadRequest(43*KB, 4*KB, nullptr),
                ReadRequest(62*KB, 4*KB, nullptr),
                ReadRequest(57*KB, 4*KB, nullptr)
  };

  ReadList expect{ReadRequest( 8*KB,  4*KB, nullptr),
                  ReadRequest(24*KB, 23*KB, nullptr),
                  ReadRequest(57*KB,  9*KB, nullptr)
  };

  run_test(list, expect,
           /*max_hole*/8*KB, /*max_size*/0, /*align*/0,
           /*do_overlap*/false, /*eof*/67*KB);
}

static void
test_align()
{
  // Same input data as test_simple but align to 4 KB boundaries.
  // test1: No holes allowed, test2: 8 KB holes allowed.
  //
  // Note: Even if padding caused two chunks to become adjacent, they
  // might not get consolidated. See _split_requests. This is an
  // academic issue because max_hole will
  //
  // Note: The code should not attempt to read past eof. This is
  // a corner case. The test data covers that case (request #6)
  //
  // Note: Even if the input has no overlaps, padding might make some.
  // The test data covers that case (request 6 and 7).
  //
  // The diagram below shows an asterisk for the data
  // that will be read and discarded
  //
  // [--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--
  // --------2222------------1111333***4444*****5555*--------*7777*6666*
  // 0000000000111111111122222222223333333333444444444455555555556666666
  // .123456789.123456789.123456789.123456789.123456789.123456789.123456

  const std::int64_t KB = 1024;

  ReadList list{ReadRequest(24*KB, 4*KB, nullptr),
                ReadRequest( 8*KB, 4*KB, nullptr),
                ReadRequest(28*KB, 3*KB, nullptr),
                ReadRequest(34*KB, 4*KB, nullptr),
                ReadRequest(43*KB, 4*KB, nullptr),
                ReadRequest(62*KB, 4*KB, nullptr),
                ReadRequest(57*KB, 4*KB, nullptr)
  };

  // Expected result given that blocks don't consolidate
  // if they became contiguous just due to padding.
  ReadList expect1{ReadRequest( 8*KB,  4*KB, nullptr),
                   ReadRequest(24*KB,  8*KB, nullptr),
                   ReadRequest(32*KB,  8*KB, nullptr),
                   ReadRequest(40*KB,  8*KB, nullptr),
                   ReadRequest(56*KB,  8*KB, nullptr),
                   ReadRequest(60*KB,  7*KB, nullptr),
  };
  run_test(list, expect1,
           /*max_hole*/0, /*max_size*/0, /*align*/4*KB,
           /*do_overlap*/false, /*eof*/67*KB);

  // Expected result given that blocks will consolidate
  // anyway, since both hole removal and align is specified.
  ReadList expect2{ReadRequest( 8*KB,  4*KB, nullptr),
                   ReadRequest(24*KB, 24*KB, nullptr),
                   ReadRequest(56*KB, 11*KB, nullptr)
  };
  run_test(list, expect2,
           /*max_hole*/8*KB, /*max_size*/0, /*align*/4*KB,
           /*do_overlap*/false, /*eof*/67*KB);
}


static void
test_overlap()
{
  // Test data:
  //   Request 1 and 3 fully overlap.
  //   Request 2 and 4 partly overlap.
  //   Request 5 and 6 are contiguous and should be merged.
  //   Request 7 is far away and should not be merged.
  //   Request 8 and 9 start at the same offset but differs in size.
  //
  // How to handle overlapping requests is still not decided.
  //
  //   - Pass do_overlap=false to allow overlaps but don't consolidate,
  //     which means the overlapping data actually gets read twice.
  //     This is simpler and also avoids issues with mutable buffers.
  //
  //   - Pass do_overlap=true to try to consolidate everything.
  //
  //   - Implement a check in the code that forbids overlapping requests.
  //     This is problematic because specifying an alignment might cause
  //     a valid list of requests to become overlapped.
  //
  // [--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--
  // --------2222--------------------66665555--------7777---------------
  // ----------44444444------------------------------888888-------------
  // ----------------------1111-----------------------------------------
  // ----------------------3333-----------------------------------------
  // 0000000000111111111122222222223333333333444444444455555555556666666
  // .123456789.123456789.123456789.123456789.123456789.123456789.123456

  const std::int64_t KB = 1024;

  ReadList list{ReadRequest(22*KB, 4*KB, nullptr),
                ReadRequest( 8*KB, 4*KB, nullptr),
                ReadRequest(22*KB, 4*KB, nullptr),
                ReadRequest(10*KB, 8*KB, nullptr),
                ReadRequest(36*KB, 4*KB, nullptr),
                ReadRequest(32*KB, 4*KB, nullptr),
                ReadRequest(48*KB, 4*KB, nullptr),
                ReadRequest(48*KB, 6*KB, nullptr)
  };

  // 5 and 6 should consolidate, the rest not because of overlap.
  // sorted order: 2,4,1,3,6+5,7,8
  ReadList expect1{ReadRequest( 8*KB,  4*KB, nullptr),
                   ReadRequest(10*KB,  8*KB, nullptr),
                   ReadRequest(22*KB,  4*KB, nullptr),
                   ReadRequest(22*KB,  4*KB, nullptr),
                   ReadRequest(32*KB,  8*KB, nullptr),
                   ReadRequest(48*KB,  4*KB, nullptr),
                   ReadRequest(48*KB,  6*KB, nullptr)
  };

  run_test(list, expect1,
           /*max_hole*/0*KB, /*max_size*/0, /*align*/0,
           /*do_overlap*/false, /*eof*/67*KB);

  // Now allow consolidating overlapping requests.
  ReadList expect2{ReadRequest( 8*KB, 10*KB, nullptr),
                   ReadRequest(22*KB,  4*KB, nullptr),
                   ReadRequest(32*KB,  8*KB, nullptr),
                   ReadRequest(48*KB,  6*KB, nullptr)
  };
  run_test(list, expect2,
           /*max_hole*/0*KB, /*max_size*/0, /*align*/0,
           /*do_overlap*/true, /*eof*/67*KB);
}

static void
test_large()
{
  // Some requests that might have been consolidated but won't be
  // because the total size got too large.
  //
  // There are still opportunities for optimizing: When a region is split
  // in two it could be split in to roughly equal parts instead of using
  // a greedy algorithm. But frankly I doubt this will make much difference.
  // With this particular test data the greedy algorithm gives the best
  // result in any case.
  //
  // [--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--][--
  // ---111222222222233333--------455555555555555555555----6------------
  // 0000000000111111111122222222223333333333444444444455555555556666666
  // .123456789.123456789.123456789.123456789.123456789.123456789.123456

  const std::int64_t KB = 1024;

  ReadList list{ReadRequest( 3*KB,  3*KB, nullptr),
                ReadRequest( 6*KB, 10*KB, nullptr),
                ReadRequest(16*KB,  5*KB, nullptr),
                ReadRequest(29*KB,  1*KB, nullptr),
                ReadRequest(30*KB, 20*KB, nullptr),
                ReadRequest(54*KB,  1*KB, nullptr)
  };

  ReadList expect1{ReadRequest( 3*KB, 13*KB, nullptr),
                   ReadRequest(16*KB,  5*KB, nullptr),
                   ReadRequest(29*KB,  1*KB, nullptr),
                   ReadRequest(30*KB, 20*KB, nullptr),
                   ReadRequest(54*KB,  1*KB, nullptr),
  };
  run_test(list, expect1,
           /*max_hole*/0, /*max_size*/16*KB, /*align*/0*KB,
           /*do_overlap*/false, /*eof*/67*KB);

}

/**
 * This is just to exercise the ConsolidateRequests::_print_requests
 * debug function.
 */
static void
test_consolidate_print()
{
  const std::int64_t KB = 1024;
  ReadList list1{ReadRequest( 8*KB, 4*KB, nullptr),
                 ReadRequest(24*KB, 7*KB, nullptr),
                 ReadRequest(34*KB, 4*KB, nullptr),
                 ReadRequest(43*KB, 4*KB, nullptr),
                 ReadRequest(57*KB, 4*KB, nullptr),
                 ReadRequest(62*KB, 4*KB, nullptr)
  };
  ReadList list2{ReadRequest(80*KB, 4*KB, nullptr),
                 ReadRequest(84*KB, 7*KB, nullptr),
  };
  ReadDoubleList groups{list1, list2};
  std::stringstream ss;
  ConsolidateRequests::_print_requests(groups, "PrintTest:", ss);
  if (verbose())
    std::cout << ss.str();
  TEST_CHECK(ss.str().find("read offset     8 KB end    12 KB size   4 KB") != std::string::npos);

  std::stringstream ss2;
  ConsolidateRequests::_print_requests(ReadDoubleList{}, "PrintTest:", ss2);
  if (verbose())
    std::cout << ss2.str();
  TEST_CHECK(ss2.str().find("(empty)") != std::string::npos);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("consolidate.zero",        test_zero);
      register_test("consolidate.one",         test_one);
      register_test("consolidate.simple",      test_simple);
      register_test("consolidate.holes",       test_holes);
      register_test("consolidate.align",       test_align);
      register_test("consolidate.overlap",     test_overlap);
      register_test("consolidate.large",       test_large);
      register_test("consolidate.print",       test_consolidate_print);
    }
  } dummy;
} // namespace for registration
