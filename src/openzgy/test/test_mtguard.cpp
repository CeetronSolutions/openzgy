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
#include "../impl/mtguard.h"
#include "../api.h" // for ProgressWithDots

#include <exception>
#include <atomic>
#include <functional>
#include <iostream>
#include <sstream>

namespace InternalZGY {
#if 0
}
#endif


static void
test_mtguard_loop()
{
  const int loopcount = 64;
  std::atomic<int> detonate(13);
  std::atomic<int> worked(0);
  try {
    MTGuard guard;
#pragma omp parallel for schedule(dynamic,1) num_threads(8)
    for (int ii=0; ii<loopcount; ++ii) {
      guard.run([&]() {
                  Test_Utils::random_delay(200);
                  worked.fetch_add(1);
                  if (detonate.fetch_sub(1) == 1)
                    throw std::runtime_error("Boom!");
                });
    }
    guard.finished();
    TEST_CHECK(false && "Did not get expected exception");
  }
  catch (const std::runtime_error& ex) {
    TEST_CHECK(ex.what() && std::string(ex.what()) == "Boom!");
  }
  if (verbose())
    std::cout << "\nmtguard.loop did "
              << worked.load() << " chunks."
              << std::endl;
}

static void
test_mtguard_ordered()
{
  const int loopcount = 64;
  std::atomic<int> detonate(13);
  std::atomic<int> worked(0);
  std::atomic<int> next(0);
  int order[loopcount]{0};
  try {
    MTGuard guard;
#pragma omp parallel num_threads(8)
    {
#pragma omp for ordered schedule(dynamic,1)
      for (int ii=0; ii<loopcount; ++ii) {
        guard.run([&]() {
                    Test_Utils::random_delay(ii==0 ? 400 : 100);
                    worked.fetch_add(1);
                    if (detonate.fetch_sub(1) == 1)
                      throw std::runtime_error("Boom!");
                  });
#pragma omp ordered
        guard.run([&]() {
                    order[next.fetch_add(1)] = ii+1;
                    Test_Utils::random_delay(20);
                  });
      }
    }
    guard.finished();
    TEST_CHECK(false && "Did not get expected exception");
  }
  catch (const std::runtime_error& ex) {
    TEST_CHECK(ex.what() && std::string(ex.what()) == "Boom!");
  }
  if (verbose()) {
    std::cout << "\nmtguard.ordered did "
              << worked.load() << " chunks."
              << std::endl;
    std::cout << "Order:";
    for (int ii=0; ii<loopcount; ++ii)
      if (order[ii] > 0)
        std::cout << " " << order[ii];
    std::cout << std::endl;
  }
  // ii==0 took extra long time, which makes it unlikely that it
  // was delivered first unless the order clause is working.
  for (int ii=1; ii<loopcount; ++ii) {
    TEST_CHECK(order[ii] == 0 || order[ii] > order[ii-1]);
  }
}

static void
test_mtguard_progress()
{
  const int loopcount = 64;
  std::atomic<int> detonate(13);
  std::atomic<int> worked(0);
  std::stringstream ss;
  OpenZGY::ProgressWithDots p1(51, ss);
  try {
    MTGuardWithProgress guard(std::ref(p1), loopcount);
#pragma omp parallel for schedule(dynamic,1) num_threads(8)
    for (int ii=0; ii<loopcount; ++ii) {
      guard.run([&]() {
                  Test_Utils::random_delay(200);
                  worked.fetch_add(1);
                  guard.progress();
                  if (detonate.fetch_sub(1) == 1)
                    throw std::runtime_error("Boom!");
                });
    }
    guard.finished();
    TEST_CHECK(false && "Did not get expected exception");
  }
  catch (const std::runtime_error& ex) {
    TEST_CHECK(ex.what() && std::string(ex.what()) == "Boom!");
  }
  if (verbose())
    std::cout << "\nmtguard.progress did "
              << worked.load() << " chunks.\n"
              << "Progress: " << ss.str()
              << std::endl;
}

static void
test_mtguard_noexcept()
{
  const int loopcount = 13;
  std::atomic<int> worked(0);
  std::stringstream ss;
  OpenZGY::ProgressWithDots p1(51, ss);
  MTGuardWithProgress guard(std::ref(p1), loopcount);
#pragma omp parallel for schedule(dynamic,1) num_threads(8)
  for (int ii=0; ii<loopcount; ++ii) {
    guard.run([&]() {
                Test_Utils::random_delay(100);
                worked.fetch_add(1);
                guard.progress();
              });
  }
  guard.finished();
  if (verbose())
    std::cout << "\nmtguard.noexcept did "
              << worked.load() << " chunks.\n"
              << "Progress: " << ss.str()
              << std::endl;
}

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("mtguard.loop",            test_mtguard_loop);
      register_test("mtguard.ordered",         test_mtguard_ordered);
      register_test("mtguard.progress",        test_mtguard_progress);
      register_test("mtguard.noexcept",        test_mtguard_noexcept);
    }
  } dummy;
}

} // namespace
