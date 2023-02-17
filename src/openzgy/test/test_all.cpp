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

#define CUTEST_HAS_MAIN
#include "test_all.h"
#include "../impl/environment.h"

#include <iostream>
#include <functional>
#ifndef _WIN32
#include <signal.h>
#endif

static bool enable_sd_tests()
{
  static bool enabled = InternalZGY::Environment::getStringEnv("OPENZGY_TOKEN") != "";
  return enabled;
}

static void do_nothing()
{
}

struct test__ test_list__[1000]{0};
int test_registered__ = 0;
void register_test(const char* name, void (*func)(void))
{
  if (test_registered__ >= 1000)
    exit(9);
  if (!func)
    func = do_nothing;
  int pos = 0;
  while (test_list__[pos].func && strcmp(name, test_list__[pos].name) >= 0)
    ++pos;
  for (int ii=test_registered__; ii>pos; --ii) {
    test_list__[ii].name = test_list__[ii-1].name;
    test_list__[ii].func = test_list__[ii-1].func;
  }
  ++test_registered__;
  test_list__[test_registered__].func = nullptr;
  test_list__[pos].name = name;
  test_list__[pos].func = func;
}

void register_sd_test(const char* name, void (*func)(void))
{
  if (enable_sd_tests())
    register_test(name, func);
  else {
    char *noname = new char[strlen(name)+2];
    noname[0] = '~';
    // Re-implementing strcpy() because, Windows.
    char *cp = &noname[1];
    while(*name)
      *cp++ = *name++;
    *cp++ = '\0';
    register_test(noname, do_nothing);
  }
}

static int _global_verbose = 0;

/**
 * If --no-exec is not in force, verbose will always be zero.
 * If --no-exec is given but no explicit test names then verbose will also
 * be zero because .quiet and .verbose are both executed once, with .quiet
 * first because tests are sorted lexically. To run individual tests with
 * verbose output you will need --no-exec .verbose test1 test2 ...
 */
int verbose()
{
  return _global_verbose < 0 ? 0 : _global_verbose;
}

static void inc_verbose()
{
  ++_global_verbose;
}

static void dec_verbose()
{
  --_global_verbose;
}

static void sigpipe_ignore()
{
#ifndef _WIN32
  // On Linux the ssl library can cause bogus SIGPIPE signals.
  // It should be safe to ignore, because in that case the
  // library simply returns a proper error code which the caller
  // apparently handles quite well.
  ::signal(SIGPIPE, SIG_IGN);
#else
  // Assuming that Windows doesn't have the same issue.
#endif
}

namespace {
  class Register
  {
  public:
    Register()
    {
      // Note: The function to be registered is an old style function pointer.
      // If there is no lambda-capture, a closure can be implicitly converted
      // to a pointer to function with the same parameter and return types.
      // Which means that a lambda with no capture should work. What I worry
      // about is what happens when I define a lambda in the argument list.
      // Will that have indefinite scope, or am I playing with undefined
      // behavior here?
      register_test(".quiet",   dec_verbose);
      register_test(".verbose", inc_verbose);
      register_test(".sigpipe", sigpipe_ignore);
      //register_test(".dummy",   0);
      //register_test(".quiet",   [](){--_global_verbose;});
      //register_test(".verbose", [](){++_global_verbose;});
    }
  } dummy;
}
