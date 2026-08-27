// Copyright 2017-2020, Schlumberger
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
#include "../impl/cornerpoints.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <cmath>

using namespace InternalZGY;

namespace {
#if 0
}
#endif

/** Arbitrary Control Point */
struct ACP
{
  float il; /** Inline */
  float xl; /** Crossline */
  double x;  /** Easting */
  double y;  /** Northing */

  ACP(float _il = 0, float _xl = 0, double _x = 0, double _y = 0)
    : il(_il), xl(_xl), x(_x), y(_y)
  {
  }
};

/** Annotation range */
struct AnnotRange
{
  float  first;
  float  inc;
  size_t count;

  AnnotRange(float _first, float _inc, size_t _count)
    : first(_first), inc(_inc), count(_count)
  {
  }

  float operator[](size_t i) const { return first + inc*i; }
};

typedef InternalZGY::OrderedCornerPoints OCP;

static OCP createOCP(const AnnotRange& ilr, const AnnotRange& xlr, const ACP* acp)
{
  return
    OCP(ilr.first, ilr.inc, ilr.count, xlr.first, xlr.inc, xlr.count,
        acp[0].il, acp[0].xl, acp[0].x, acp[0].y,
        acp[1].il, acp[1].xl, acp[1].x, acp[1].y,
        acp[2].il, acp[2].xl, acp[2].x, acp[2].y);
}

static void testInc(float ilinc, float xlinc)
{
  // TODO-Test, maybe different numbers in IL and XL?
  const float il0 = 1.0f;
  const float xl0 = 1.0f;
  const size_t ilcnt = 101;
  const size_t xlcnt = 101;

  ACP acp[3] =
  {
    ACP(il0,              xl0,              0.00, 0.00),
    ACP(il0,              xl0 + xlcnt - 1,  1.00, 0.00),
    ACP(il0 + ilcnt - 1,  xl0,              0.00, 1.00)
  };

  AnnotRange
    ilr(il0 + 0.67f*(ilcnt - 1), ilinc, (ilcnt - 1)/4 + 1),
    xlr(xl0 + 0.31f*(xlcnt - 1), xlinc, (xlcnt - 1)/4 + 1);

  OCP ocp = createOCP(ilr, xlr, acp);

  const double errtol = 1e-7;
  const double ocpexp[4][2] =
  {
    { (xlr.first - xl0)/(xlcnt - 1),                           (ilr.first - il0)/(ilcnt - 1) },
    { (xlr.first - xl0)/(xlcnt - 1),                           (ilr.first + ilr.inc*(ilr.count - 1) - il0)/(ilcnt - 1) },
    { (xlr.first + xlr.inc*(xlr.count - 1) - xl0)/(xlcnt - 1), (ilr.first - il0)/(ilcnt - 1) },
    { (xlr.first + xlr.inc*(xlr.count - 1) - xl0)/(xlcnt - 1), (ilr.first + ilr.inc*(ilr.count - 1) - il0)/(ilcnt - 1) }
  };

  for (size_t i = 0; i < 4; ++i) {
    const double dx = ocp[i].x - ocpexp[i][0];
    const double dy = ocp[i].y - ocpexp[i][1];
    const double diff = std::sqrt(dx*dx + dy*dy);
    TEST_CHECK(diff < errtol);
  }
}

static void test_1()
{
  testInc(1.0f, 1.0f);
}

static void test_2()
{
  testInc(1.0f, -1.0f);
}

static void test_3()
{
  testInc(-1.0f, 1.0f);
}

static void test_4()
{
  testInc(-1.0f, -1.0f);
}

static void testSingleOCP(const AnnotRange& ilr, const AnnotRange& xlr, const ACP* acp, const OCP::Element* exp, const double& errtol)
{
  OCP ocp = createOCP(ilr, xlr, acp);

  for (size_t i = 0; i < 4; ++i) {

    // verify bulk-data indices
    TEST_CHECK(ocp[i].i == exp[i].i);
    TEST_CHECK(ocp[i].j == exp[i].j);

    // verify annotation indices
    TEST_CHECK(ocp[i].il == exp[i].il);
    TEST_CHECK(ocp[i].xl == exp[i].xl);

    // verify map projection coordinates
    const double dx = ocp[i].x - exp[i].x;
    const double dy = ocp[i].y - exp[i].y;
    const double diff = std::sqrt(dx*dx + dy*dy);
    TEST_CHECK(diff <= errtol);
  }
}

static void testFwdRevCrosslinesOCP(const AnnotRange& ilr, const AnnotRange& xlr, const ACP* acp, const OCP::Element* exp, const double& errtol)
{
  // test forward crosslines
  testSingleOCP(ilr, xlr, acp, exp, errtol);

  // test reverse crosslines
  AnnotRange rxlr(xlr[xlr.count - 1], -xlr.inc, xlr.count);
  OCP::Element rexp[4] = { exp[2], exp[3], exp[0], exp[1] }; // expected OCPs in permuted order corresponding to reversing the bulk-data ordering along the inline/i axis
  for (size_t i = 0; i < 4; ++i) {  // restore original (i, j) order (only annotation and coordinates are to be permuted)
    rexp[i].i = exp[i].i;
    rexp[i].j = exp[i].j;
  }
  testSingleOCP(ilr, rxlr, acp, rexp, errtol);
}

static void testOCP(const AnnotRange& ilr, const AnnotRange& xlr, const ACP* acp, size_t nacp, const OCP::Element* exp, const double& errtol)
{
  for (size_t ii = 2; ii < nacp; ++ii, ++acp) {

    // test forward inlines
    testFwdRevCrosslinesOCP(ilr, xlr, acp, exp, errtol);

    // test reverse inlines
    AnnotRange rilr(ilr[ilr.count - 1], -ilr.inc, ilr.count);
    OCP::Element rexp[4] = { exp[1], exp[0], exp[3], exp[2] }; // expected OCPs in permuted order corresponding to reversing the bulk-data ordering along the crossline/j axis
    for (size_t i2 = 0; i2 < 4; ++i2) {  // restore original (i, j) order (only annotation and coordinates are to be permuted)
      rexp[i2].i = exp[i2].i;
      rexp[i2].j = exp[i2].j;
    }
    testFwdRevCrosslinesOCP(rilr, xlr, acp, rexp, errtol);
  }
}

static void test_5()
{
  const AnnotRange ilr(4001, 1, 425);
  const AnnotRange xlr(2900, 1, 601);
  const ACP acp[4] = {
    ACP(  4001.00,   2900.00,  2664280.00,  9885031.00),
    ACP(  4001.00,   3500.00,  2654069.00,  9907418.00),
    ACP(  4425.00,   2900.00,  2638967.00,  9873485.00),
    ACP(  4425.00,   3500.00,  2628756.00,  9895872.00)
  };
  OCP::Element exp[4];
  exp[OCP::Min0Min1] = OCP::Element(0,                         0, ilr[0],             xlr[0],             2664280.00, 9885031.00);
  exp[OCP::Min0Max1] = OCP::Element(0,             xlr.count - 1, ilr[0],             xlr[xlr.count - 1], 2654069.00, 9907418.00);
  exp[OCP::Max0Min1] = OCP::Element(ilr.count - 1,             0, ilr[ilr.count - 1], xlr[0],             2638967.00, 9873485.00);
  exp[OCP::Max0Max1] = OCP::Element(ilr.count - 1, xlr.count - 1, ilr[ilr.count - 1], xlr[xlr.count - 1], 2628756.00, 9895872.00);
  const double errtol(0.5);
  testOCP(ilr, xlr, acp, 4, exp, errtol);
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("corners.test1",           test_1);
      register_test("corners.test2",           test_2);
      register_test("corners.test3",           test_3);
      register_test("corners.test4",           test_4);
      register_test("corners.test5",           test_5);
    }
  } dummy;
} // namespace for registration
