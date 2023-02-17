// Copyright 2017-2022, SLB
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

#include "wintools.h"

#include <atomic>
#include <string>
#include <utility>
#include <iostream>
#include <thread>
#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#define NOMINMAX
#include <Windows.h>
#endif

#include "fancy_timers.h"

namespace InternalZGY {
  namespace WindowsTools {
#if 0
  }
}
#endif

/**
 * Global number of threads created and destroyed.
 * For debugging only. Depending on the ifdefs later in
 * this file the values might remain zero.
 */
static std::atomic<int>threada{}, threadd{};

std::pair<int, int>
threadReportPair()
{
  int detached = threadd.fetch_add(0);
  int attached = threada.fetch_add(0);
  threada.fetch_sub(attached);
  threadd.fetch_sub(detached);
  return std::pair<int, int>(attached, detached);
}

std::string
threadReportString()
{
  std::pair<int, int> info = threadReportPair();
  return ("+" + std::to_string(info.first) +
    " -" + std::to_string(info.second));
}

std::atomic<int> noop_called(0);
void noop()
{
  noop_called.fetch_add(1);
}

void threadTimeTest(int count)
{
  std::cerr << "Start creating " << count << " threads" << std::endl;
  InternalZGY::WindowsTools::threadReportString();
  InternalZGY::SummaryPrintingTimerEx timer("CreateThreads");
  for (int ii = 0; ii < count; ++ii) {
    InternalZGY::SimpleTimer tt(timer);
    std::thread t(noop);
    t.join();
  }
  std::cerr << "Done  creating " << count << " threads" << " reported " << threadReportString() << std::endl;
}

#if defined(_WIN32) && 1

/**
 * Get notified when any thread is created or destroyed in the app.
 * Windows-only version. Needs to be in a dll, not the main application.
 * There must not be another DllMain in the same dll.
 * Caveat: Possibly too much overhead to allow leaving it in production code.
 */
extern "C" BOOL WINAPI
DllMain(HINSTANCE hinstDLL, DWORD fdwReason, LPVOID lpvReserved)
{
  switch (fdwReason)
  {
  case DLL_THREAD_ATTACH:
    threada.fetch_add(1);
    break;
  case DLL_THREAD_DETACH:
    threadd.fetch_add(1);
    break;
  }
  return TRUE;
}

#elif 1

/**
 * Get notified when any thread is created or destroyed in the app.
 * Portable version. Relies on the compiler having good support for
 * thread local variables with a non trivial type. So that the variable
 * gets constructed and destructed properly.
 * Caveat: Possibly too much overhead to allow leaving it in production code.
 * Caveat: Even with good TLS support, this is not guaranteed to work.
 * Because the compiler is allowed to defer constructing the "hello"
 * variable until its first use in the thread.
 */
class Hello
{
public:
  Hello()
  {
    threada.fetch_add(1);
  }
  ~Hello()
  {
    threadd.fetch_add(1);
  }
};
static thread_local Hello hello;

#endif

}} // namespace
