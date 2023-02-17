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
#include "test_utils.h"
#include "../impl/timer.h"
#include "../impl/fancy_timers.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <omp.h>

using namespace InternalZGY;

namespace {
#if 0
}
#endif

static void
test_timer_oneshot()
{
  Timer t(true, "oneshot");
  Test_Utils::random_delay(5);
  t.stop();
  t.start();
  Test_Utils::random_delay(50);
  TEST_CHECK(t.getRunning());
  TEST_CHECK(t.getEnabled());
  t.stop();
  TEST_CHECK(std::string(t.getName()) == std::string("oneshot"));
  TEST_CHECK(!t.getRunning());
  TEST_CHECK(t.getTotal() >= 0.020);
  TEST_CHECK(t.getTotal() <= 1.0);
  TEST_CHECK(t.getLast()  < t.getTotal());
  TEST_CHECK(t.getSkip()  == 0);
  TEST_CHECK(t.getCount() == 2);
  std::string msg = t.getValue(true, false);
  TEST_CHECK(msg.find("Time for oneshot") != std::string::npos);
  TEST_CHECK(msg.find("in 2 calls") != std::string::npos);

  Timer t2(false, "disabled");
  TEST_CHECK(!t2.getRunning());
  TEST_CHECK(!t2.getEnabled());
  t2.stop();
  TEST_CHECK(t2.getCount() == 0);
  TEST_CHECK(t2.getTotal() == 0);
  TEST_CHECK(strlen(t2.getValue(true, false))==0);
}

static void
test_timer_summary()
{
  SummaryTimer summary("MySummary");
  {
    SimpleTimer t(summary);
    Test_Utils::random_delay(50);
    t.stop();
  }

  TEST_CHECK(std::string(summary.getName()) == std::string("MySummary"));
  TEST_CHECK(summary.getTotal() >= 0.020);
  TEST_CHECK(summary.getTotal() <= 1.0);
  TEST_CHECK(summary.getLast()  == summary.getTotal()); // because count is 1.
  TEST_CHECK(summary.getCount() == 1);
  std::string msg = summary.getValue(true, false);
  TEST_CHECK(msg.find("Time for MySummary") != std::string::npos);
  msg = summary.getCSV();
  TEST_CHECK(msg.find("\"MySummary\"") != std::string::npos);
}

static void
test_timer_printing()
{
  SummaryPrintingTimerEx summary("YourSummary");
  bool parallel_ok = false;
#pragma omp parallel num_threads(2)
  {
    SimpleTimer t(summary);
    Test_Utils::random_delay(50);
    t.stop();
    if (omp_get_thread_num() == 0 && omp_get_num_threads() > 1)
      parallel_ok = true;
  }
  TEST_CHECK(summary.getTotal() >= 0.020);
  TEST_CHECK(summary.getTotal() <= 1.0);
  TEST_CHECK(summary.getCount() == 2);
  TEST_CHECK(!parallel_ok || summary.getAdjusted() < summary.getTotal());
  std::stringstream ss, sscsv;
  summary.printToFile(ss, false, false);
  summary.printToFile(sscsv, true, true);
  std::string msg = summary.getValue(true, false);
  TEST_CHECK(ss.str().find("Time for YourSummary") != std::string::npos);
  TEST_CHECK(ss.str().find("adjusted") != std::string::npos);
  TEST_CHECK(sscsv.str().find("\"YourSummary\"") != std::string::npos);
}

static void
test_timer_nicesize()
{
  const std::int64_t K{1024};
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(0)         == "0 bytes");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(-42)       == "-42 bytes");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(3*K)       == "3 KB");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(3*K+1)     == "3073 bytes");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(7*K*K)     == "7 MB");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(-1*K*K*K)  == "-1 GB");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(2*K*K*K*K) == "2 TB");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(2*K*K*K+K) == "2097153 KB");
  TEST_CHECK(SummaryPrintingTimerEx::niceSize(7*256*K)   == "1792 KB");
  // The next one used to be handles as a special case: An even number
  // of 256 KB bricks were printed as MB even though that might involve
  // a decimal point.
  //TEST_CHECK(SummaryPrintingTimerEx::niceSize(7*256*K) == "1.75 MB");
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("timer.oneshot",           test_timer_oneshot);
      register_test("timer.summary",           test_timer_summary);
      register_test("timer.printing",          test_timer_printing);
      register_test("timer.nicesize",          test_timer_nicesize);
    }
  } dummy;
} // namespace for registration
