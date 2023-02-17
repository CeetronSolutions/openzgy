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
#include "../iocontext.h"

#include <iostream>
#include <sstream>
#include <memory>

using namespace OpenZGY;
using Test_Utils::must_throw;

namespace Test {
  // Friend of class SeismicStoreIOContext.
  class TestIOContext {
  public:
    static void test_defaults();
    static void test_setters();
    static void test_errors();
  };
}

void Test::TestIOContext::test_defaults()
{
  SeismicStoreIOContext ctx;
  if (verbose())
    std::cout << ctx.toString();
  // Not true if OPENZGY_ variables set.
  //TEST_CHECK(ctx._sdurl == "");
  //TEST_CHECK(ctx._sdapikey == "");
  // Not true if OPENZGY_TOKEN is set.
  //TEST_CHECK(ctx._sdtoken == "");
  TEST_EQUAL(ctx._maxsize,   64 * 1024*1024);
  TEST_EQUAL(ctx._maxhole,    2 * 1024*1024);
  TEST_EQUAL(ctx._aligned,    0 * 1024*1024);
  TEST_EQUAL(ctx._segsize,  256 * 1024*1024 * ctx._segsplit);
  TEST_EQUAL(ctx._segsplit, 8);
  TEST_EQUAL(ctx._iothreads, 1);
  TEST_EQUAL(ctx._cputhreads, 1);
  TEST_EQUAL(ctx._writethreads, 1);
  TEST_EQUAL(ctx._legaltag, "");
  TEST_EQUAL(ctx._writeid, "");
  TEST_EQUAL(ctx._seismicmeta, "");
  ctx.writethreads(42);
  TEST_EQUAL(ctx._segsplit, 1);
  TEST_EQUAL(ctx._writethreads, 42);
  ctx.segsplit(19);
  // Not yet decided whether this will cause the code
  // to switch back to the old "segsplit" mode or
  // just ignore the call. Currently it is ignored.
  TEST_EQUAL(ctx._segsplit, 1);
  TEST_EQUAL(ctx._writethreads , 42);
}

void Test::TestIOContext::test_setters()
{
  SeismicStoreIOContext ctx;
  ctx.sdurl("https://seismicstore.com")
    .sdapikey("secret-api-key")
    .sdtoken("algo.payload.signature", "mocktoken")
    .maxsize(42)
    .maxhole(7)
    .aligned(1)
    .segsize(15)
    .segsplit(3)
    .iothreads(8)
    .cputhreads(5)
    .legaltag("illegal")
    .writeid("WID")
    .seismicmeta("{\"foo\": 42}");
  if (verbose())
    std::cout << ctx.toString();
  TEST_CHECK(ctx._sdurl == "https://seismicstore.com");
  TEST_CHECK(ctx._sdapikey.substr(0,3) == "sec");
  TEST_CHECK(ctx._sdtoken.substr(0,3) == "alg");
  TEST_CHECK(ctx._maxsize ==   42 * 1024*1024);
  TEST_CHECK(ctx._maxhole ==    7 * 1024*1024);
  TEST_CHECK(ctx._aligned ==    1 * 1024*1024);
  TEST_CHECK(ctx._segsize ==   15 * 1024*1024 * ctx._segsplit);
  TEST_CHECK(ctx._segsplit ==  3);
  TEST_CHECK(ctx._iothreads == 8);
  TEST_CHECK(ctx._cputhreads== 5);
  TEST_CHECK(ctx._legaltag == "illegal");
  TEST_CHECK(ctx._writeid == "WID");
  TEST_CHECK(ctx._seismicmeta == "{\"foo\": 42}");
  std::string str = ctx.toString();
  TEST_CHECK(str.find("maxsize:  42 MB") != std::string::npos);
}

void Test::TestIOContext::test_errors()
{
  SeismicStoreIOContext ctx;
  TEST_CHECK(must_throw("must be between", [&](){ctx.maxsize(2000);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.maxhole(2001);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.aligned(2000);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.segsize(200000);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.segsplit(0);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.iothreads(0);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.cputhreads(0);}));
  TEST_CHECK(must_throw("must be between", [&](){ctx.writethreads(0);}));
}

namespace {
  class Register
  {
  public:
    Register()
    {
#ifndef TEST_WRITETHREADS // If writethreads forced on, this will fail.
      register_test("iocontext.defaults",            Test::TestIOContext::test_defaults);
      register_test("iocontext.setters",             Test::TestIOContext::test_setters);
      register_test("iocontext.errors",              Test::TestIOContext::test_errors);
#endif
    }
  } dummy;
} // namespace for registration
