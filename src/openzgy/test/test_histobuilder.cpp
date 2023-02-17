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
#include "../impl/histogrambuilder.h"

#include <iostream>
#include <memory>
#include <list>
#include <algorithm>
#include <limits>
#include <cmath>
#include <vector>

//using namespace OpenZGY;
using namespace InternalZGY;

#define CPPUNIT_ASSERT(a) TEST_CHECK((a))
#define CPPUNIT_ASSERT_EQUAL(a,b) TEST_EQUAL(b,a)
#define CPPUNIT_ASSERT_DOUBLES_EQUAL(a,b,eps) TEST_EQUAL_FLOAT(b,a,eps)

namespace {
#if 0
}
#endif

#if 0
static void
PrintHistogram(const char *msg, const ZGY_NS::HistogramData& h)
{
  printf("%srange %+g..%+g in %d bins\n", msg, h.getmin(), h.getmax(), h.getsize());
  double width = h.getsize() <= 1 ? 0 : (h.getmax() - h.getmin()) / (h.getsize() - 1);
  double begin = h.getmin();
  for (int ii=0; ii<h.getsize(); ++ii)
    if (h.getbins()[ii] != 0)
      printf("bin[%d]: %d value %+g (%+g..%+g)\n", ii, (int)h.getbins()[ii], begin+ii*width, begin+(ii-0.5)*width, begin+(ii+0.5)*width);
}
#endif

/**
 * Test the default constructor for the StatisticsBuilder.
 * Does not depend on type, so only needs to be run once.
 */
static void
TestStatisticDataDefaultConstructor()
{
  const StatisticData s;

  CPPUNIT_ASSERT_EQUAL(0LL, s.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, s.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, s.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, s.getssq(), 1e-5);
}

//      register_test("stats.histogram_ctor",    TestStatisticDataHistogramCtor);
//      register_test("stats.add",               TestStatisticDataAdd);
//      register_test("stats.arithmetic",        TestStatisticDataArithmetic);
//      register_test("histo.normal_ctor",       TestHistogramBuilderNormalCtor);
//      register_test("histo.discrete_ctor",     TestHistogramBuilderDiscreteCtor);
//      register_test("histo.add_nan",           TestHistogramBulderAddNaN);
//      register_test("histo.copy_assign",       TestHistogramBuilderCopyAssign);
//      register_test("histo.compare",           TestHistogramBuilderCompare);
//      register_test("histo.add",               TestHistogramBuilderAdd);
//      register_test("histo.arithmetic",        TestHistogramBuilderArithmetic);

/**
 * Check copy constructor and assignment operator for StatisticData.
 * Also exercise the "copy" constructor that is passed in all the
 * statistics as discrete values, and just builds the instance from that
 * without doing any calculation.
 */
static void
TestStatisticDataCopyAssign()
{
  StatisticData a(5, 1, 42.0, 9000.0, -12.5, 19.3); // copy fom unpacked values
  StatisticData b(a);                               // using normal copy constructor
  StatisticData c;                                  // using assignment operator
  c = a;

  CPPUNIT_ASSERT_EQUAL(5LL, a.getcnt());
  CPPUNIT_ASSERT_EQUAL(1LL, a.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   a.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, a.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  a.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   a.getmax(), 1e-5);

  CPPUNIT_ASSERT_EQUAL(5LL, b.getcnt());
  CPPUNIT_ASSERT_EQUAL(1LL, b.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   b.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, b.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  b.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   b.getmax(), 1e-5);

  CPPUNIT_ASSERT_EQUAL(5LL, c.getcnt());
  CPPUNIT_ASSERT_EQUAL(1LL, c.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   c.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, c.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  c.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   c.getmax(), 1e-5);
}

static void
TestStatisticDataHistogramCtor()
{
  StatisticData::count_type hist[10] = { 0, 0, 3, 7, 2, 0, 0, 1, 0, 0 };
  //                                     1  2  3  4  5  6  7  8  9 10
  StatisticData s(hist, 10, 1, 10, true);
  StatisticData::count_type expect_cnt = hist[0] + hist[1] + hist[2] + hist[3] + hist[4] + hist[5] + hist[6] + hist[7] + hist[8] + hist[9];
  double expect_sum = double(hist[0]*1 + hist[1]*2 + hist[2]*3 + hist[3]*4 + hist[4]*5 + hist[5]*6 + hist[6]*7 + hist[7]*8 + hist[8]*9 + hist[9]*10);
  double expect_ssq = double(hist[0]*1 + hist[1]*2*2 + hist[2]*3*3 + hist[3]*4*4 + hist[4]*5*5 + hist[5]*6*6 + hist[6]*7*7 + hist[7]*8*8 + hist[8]*9*9 + hist[9]*10*10);

  CPPUNIT_ASSERT_EQUAL(expect_cnt, s.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, s.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expect_sum, s.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expect_ssq, s.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3, s.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8, s.getmax(), 1e-5);
}

static void
TestStatisticDataAdd()
{
  StatisticData s;

  CPPUNIT_ASSERT_EQUAL(0LL, s.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, s.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, s.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, s.getssq(), 1e-5);

  s.add(29.0);

  CPPUNIT_ASSERT_EQUAL(1LL, s.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, s.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29.0, s.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29.0*29.0, s.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29.0, s.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(29.0, s.getmax(), 1e-5);

  s.add(-13.17);
  s.add(101.0);

  CPPUNIT_ASSERT_EQUAL(3LL, s.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, s.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(116.83, s.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11215.45, s.getssq(), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-13.17, s.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(101.0, s.getmax(), 1e-5);
}

/**
 * Test the supported operators (+=, -=, *=) of StatisticsBulder.
 */
static void
TestStatisticDataArithmetic()
{

  // iterate each combination of empty and non-empty left-hand and right-hand operands
  const double a(7.11), b(-13.17), c(19.23);        // parameters for generating values
  size_t index[2];                                  // iteration index
  for (index[0] = 0; index[0] < 2; ++index[0]) {    // empty and non-emtpy left-hand operand
    for (index[1] = 0; index[1] < 2; ++index[1]) {  // empty and non-empty right-hand operand

      // intialize locally computed stats
      StatisticData::count_type cnt(0);
      double sum(0), ssq(0);
      double min = 0, max = 0; // initialized in case cnt==0

      // create operands
      std::auto_ptr<StatisticData> s[2];
      s[0].reset(new StatisticData);
      s[1].reset(new StatisticData);
      for (size_t k = 0; k < 2; ++k) { // left-hand and right-hand operand
        if (index[k] != 0) {
          // this (left or right) operand should not be empty, it should contain a single value.
          const double dvalue = a + b*index[0] + c*index[1] + k;
          s[k]->add(dvalue);

          // update locally computed stats
          ++cnt;
          sum += dvalue;
          ssq += (dvalue*dvalue);
          if (cnt == 1) {
            min = max = dvalue;
          }
          else {
            min = std::min(min, dvalue);
            max = std::max(max, dvalue);
          }
        }
      }

      // perform addition
      *s[0] += *s[1];

      // verify
      CPPUNIT_ASSERT_EQUAL(cnt, s[0]->getcnt());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, s[0]->getsum(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(ssq, s[0]->getssq(), 1e-6);
      if (cnt > 0) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(min, s[0]->getmin(), 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(max, s[0]->getmax(), 1e-6);
      }

      // try to subtract the second value again.
      *s[0] -= *s[1];
      CPPUNIT_ASSERT_EQUAL(cnt, s[0]->getcnt()+s[1]->getcnt());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, s[0]->getsum()+s[1]->getsum(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(ssq, s[0]->getssq()+s[1]->getssq(), 1e-6);

      // Subtracting self should result in an empty range.
      StatisticData s2;
      s2 = *s[1];
      s2 -= *s[1];
      CPPUNIT_ASSERT_EQUAL(0LL, s2.getcnt());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, s2.getsum(), 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, s2.getssq(), 1e-6);
      if (s[1]->getcnt() != 0) {
        // test makes no sense if s[1] is already empty...
        // calculate -s[1].
        s2 -= *s[1];
        CPPUNIT_ASSERT_EQUAL(-s[1]->getcnt(), s2.getcnt());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-s[1]->getsum(), s2.getsum(), 1e-6);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(-s[1]->getssq(), s2.getssq(), 1e-6);
        s2 += *s[0];
        // demonstrate that (-s1 + s0) == (s0 - s1)
        CPPUNIT_ASSERT_EQUAL(cnt, s[0]->getcnt() + s[1]->getcnt());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, s[0]->getsum() + s[1]->getsum(), 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ssq, s[0]->getssq() + s[1]->getssq(), 1e-3);
      }
    }
  }
}

static void
TestHistogramBuilderNormalCtor()
{
  const HistogramBuilder h(64, 0, 100);

  // Check the statistics part
  CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getssq(), 1e-5);

  // Check the histogram part
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(100, h.gethisto().getmax(), 1e-5);
  CPPUNIT_ASSERT_EQUAL(64, h.gethisto().getsize());
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().getcount());
}

static void
TestHistogramBuilderDiscreteCtor()
{
  HistogramBuilder::count_type hist[10] = { 0, 0, 3, 7, 2, 0, 0, 1, 0, 0 };
  const HistogramBuilder h(hist, 10, 1, 10, 13, 100, 1000, 3, 8);

  // Check the statistics part. Note that sum etc. are wrong,
  // but it will simply believe the information I passed.
  CPPUNIT_ASSERT_EQUAL(13LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(100,  h.getstats().getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1000, h.getstats().getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3,    h.getstats().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8,    h.getstats().getmax(), 1e-5);

  // Check the histogram part
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, h.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10, h.gethisto().getmax(), 1e-5);
  for (int ii=0; ii<10; ++ii) {
    CPPUNIT_ASSERT_EQUAL(hist[ii], h.gethisto().getbins()[ii]);
    CPPUNIT_ASSERT_EQUAL(hist[ii], h.gethisto().get(ii+1));
  }
  CPPUNIT_ASSERT_EQUAL(13LL, h.gethisto().getcount());
}

static void
TestHistogramBulderAddNaN()
{
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  volatile float zero = 0.0;
  const float Inf = 1.0f / zero;

  static float data[3] = { NaN, Inf, -Inf };
  HistogramBuilder h(256, 0, 0); // invalid limits
  h.add(&data[0], &data[0]); // empty range
  h.add(&data[0], &data[3]); // good range but no samples

  // Check the statustics part.
  CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(3LL, h.getstats().getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getssq(), 1e-5);

  // Check the histogram part - it should have no samples.
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().getcount());

  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().get(0)); // range not established yet, will this crash?
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().get(NaN));
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().get(Inf));
}

/**
 * Check copy constructor and assignment operator for StatisticData.
 * Also exercise the "copy" constructor that is passed in all the
 * statistics as discrete values, and just builds the instance from that
 * without doing any calculation.
 */
static void
TestHistogramBuilderCopyAssign()
{
  // Value range will be set to 1..10, so the bins are:
  //                                        1  2  3  4  5  6  7  8  9 10
  HistogramBuilder::count_type data[10] = { 0, 0, 3, 7, 2, 0, 0, 1, 0, 0 };

  // For this test, pass in statistics that clearly don't match the histogram.
  // The simple copy/assign should not care, it should just store what it gets.
  HistogramBuilder ha(data, 10, -1, 100, 5, 42.0, 9000.0, -12.5, 19.3); // copy fom unpacked values
  HistogramBuilder hb(ha);                                              // using normal copy constructor
  HistogramBuilder hc(10, -1, +1);                                      // using assignment operator
  hc = ha;

  const StatisticData& sa = ha.getstats();
  const StatisticData& sb = hb.getstats();
  const StatisticData& sc = hc.getstats();

  // Make sure the implementations don't share data pointers
  CPPUNIT_ASSERT(ha.gethisto().getbins() != hb.gethisto().getbins());
  CPPUNIT_ASSERT(ha.gethisto().getbins() != hc.gethisto().getbins());

  // Histogram part of "a"
  CPPUNIT_ASSERT_EQUAL(10, ha.gethisto().getsize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, ha.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(100, ha.gethisto().getmax(), 1e-5);
  for (int ii=0; ii<ha.gethisto().getsize(); ++ii)
    CPPUNIT_ASSERT_EQUAL(data[ii], ha.gethisto().getbins()[ii]);

  // Statistics part of "a"
  CPPUNIT_ASSERT_EQUAL(5LL, sa.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, sa.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   sa.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, sa.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  sa.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   sa.getmax(), 1e-5);

  // Histogram part of "b"
  CPPUNIT_ASSERT_EQUAL(10, hb.gethisto().getsize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, hb.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(100, hb.gethisto().getmax(), 1e-5);
  for (int ii=0; ii<hb.gethisto().getsize(); ++ii)
    CPPUNIT_ASSERT_EQUAL(data[ii], hb.gethisto().getbins()[ii]);

  // Statistics part of "b"
  CPPUNIT_ASSERT_EQUAL(5LL, sb.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, sb.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   sb.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, sb.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  sb.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   sb.getmax(), 1e-5);

  // Histogram part of "c"
  CPPUNIT_ASSERT_EQUAL(10, hc.gethisto().getsize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1, hc.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(100, hc.gethisto().getmax(), 1e-5);
  for (int ii=0; ii<hc.gethisto().getsize(); ++ii)
    CPPUNIT_ASSERT_EQUAL(data[ii], hc.gethisto().getbins()[ii]);

  // Statistics part of "c"
  CPPUNIT_ASSERT_EQUAL(5LL, sc.getcnt());
  CPPUNIT_ASSERT_EQUAL(0LL, sc.getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0,   sc.getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9000.0, sc.getssq(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.5,  sc.getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.3,   sc.getmax(), 1e-5);
}

/**
 * Test compaison of two histograms.
 * Technically they don't need to be the same size in order to match.
 * But in practice the bin width must be very close, otherwise a single bin
 * in one of the histograms might map to several bins in the other.
 * It is Ok for the other histogram to have the same bin width but
 * a larger range (and more bins).
 */
static void
TestHistogramBuilderCompare()
{
  HistogramBuilder ha(128, 0, 127);
  HistogramBuilder hb(129, 0, 127);
  HistogramBuilder hc(256, -128, 127);
  HistogramBuilder hd(256, 0, 127);

  CPPUNIT_ASSERT(ha == hb);
  CPPUNIT_ASSERT(ha == hc);
  CPPUNIT_ASSERT(ha == hd);

  float data[5] = { 12.1f, 20.0f, 20.1f, 21.0f, 50.0f };
  ha.add(&data[0], &data[5]);
  hb.add(&data[0], &data[5]);
  hc.add(&data[0], &data[5]);
  hd.add(&data[0], &data[5]);

#if VERBOSE_PRINT
  PrintHistogram("\nA: ", ha);
  PrintHistogram("B: ", hb);
  PrintHistogram("C: ", hc);
  PrintHistogram("D: ", hd);
#endif

  CPPUNIT_ASSERT(ha == hb);
  CPPUNIT_ASSERT(!(ha != hb));
  CPPUNIT_ASSERT(ha == hc);
  CPPUNIT_ASSERT(!(ha != hc));
  CPPUNIT_ASSERT(ha != hd);
  CPPUNIT_ASSERT(!(ha == hd));

  // There is special case handling for empty histograms
  HistogramBuilder empty1(12, -1, +1);
  HistogramBuilder empty2(24, -1, +1);
  HistogramBuilder empty3(48, 0, 100);
  CPPUNIT_ASSERT(empty1 == empty2);
  CPPUNIT_ASSERT(empty1 == empty3);

  CPPUNIT_ASSERT(ha != empty1);
  CPPUNIT_ASSERT(empty1 != hb);
}

/**
 * Test the Add (single sample) and DoubleRange methods.
 * There are currently not exposed publicly, so the test is incomplete.
 */
static void
TestHistogramBuilderAdd()
{
  HistogramBuilder h(200, -100, 99);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(-100, h.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(99, h.gethisto().getmax(), 1e-5);
  CPPUNIT_ASSERT_EQUAL(200, h.gethisto().getsize());
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().getcount());

  double data[4] = { 29.0, -13.17, -13.34, 98.7 };
  h.add(&data[0], &data[4]);

  CPPUNIT_ASSERT_EQUAL(4LL, h.gethisto().getcount());
  CPPUNIT_ASSERT_EQUAL(2LL, h.gethisto().getbins()[87]);
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().getbins()[129]);
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().getbins()[199]);

  CPPUNIT_ASSERT_EQUAL(2LL, h.gethisto().get(-13));
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().get(29));
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().get(99));
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().get(105)); // outside histogram
}

/**
 * Test creating histograms from iterators with varying number of samples.
 */
template <typename T>
static void
TestHistogramBuilderFromRange()
{
  // create various iterator range sizes to test
  const StatisticData::count_type size[8] = { 0, 1, 2, 3, 5, 7, 11, 31 };

  // test each size
  for (size_t i = 0; i < 8; ++i) {

    // generate list of values
    const HistogramBuilder::count_type cnt(size[i]);
    std::list<T> values;
    const double a(7.11), b(-13.17), c(19.23); // parameters for generating values
    double sum(0), ssq(0);
    double min = 0, max = 0; // need to be initialized explicitly in case cnt==0.
    for (HistogramBuilder::count_type j = 0; j < cnt; ++j) {

      // generate value
      const T value(static_cast<T>(a + b*j + c*j*j));

      const double dvalue(static_cast<double>(value));

      // update locally calculated statistics
      if (j == 0) {
        min = max = dvalue;
      }
      else {
        min = std::min(min, dvalue);
        max = std::max(max, dvalue);
      }
      sum += dvalue;
      ssq += (dvalue*dvalue);

      // add to list of values
      values.push_back(value);
    }

    // create histogram and add values in list.
    // The histogram range will then be set to -128/+127 and not the conventional
    // full range of integral type. Because I have not bothered to scale the test
    // data, so if I did it that way then for <int> all my samples would end up
    // in a single bin. Note that the choice of -128/+127 is important: When type
    // is "char", this will trigger some shortcut code.
    HistogramBuilder h(256, -128, +127);
    h.add(values.begin(), values.end());

    // verify contents
    CPPUNIT_ASSERT_EQUAL(cnt, h.getstats().getcnt());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sum, h.getstats().getsum(), 1e-5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ssq, h.getstats().getssq(), 1e-3);
    if (cnt == 0) {
      CPPUNIT_ASSERT_EQUAL(cnt, h.gethisto().getcount());
    }
    else {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(min, h.getstats().getmin(), 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(max, h.getstats().getmax(), 1e-5);
      CPPUNIT_ASSERT(h.gethisto().getcount() > 0);
    }
  }
}

/**
 * Invoke all the templated unit tests for a given template type.
 * This is a convenience just to avoid having to list N tests * M types
 * in the test suite.
 */
template <typename T>
static void
TestByType()
{
  TestHistogramBuilderFromRange<T>();

  // Currently not supported for t=double.
  // fails due to inaccuracy caused by T -> float -> double cast in StatisticBuilder class
  // also class HistogramBuilder will need some template specialization for double
  // Note: might work now but yagni.
}

static void test_char()   { TestByType<char>(); }
static void test_schar()  { TestByType<signed char>(); }
static void test_uchar()  { TestByType<unsigned char>(); }
static void test_short()  { TestByType<short>(); }
static void test_ushort() { TestByType<unsigned short>(); }
static void test_int()    { TestByType<int>(); }
static void test_uint()   { TestByType<unsigned int>(); }
static void test_float()  { TestByType<float>(); }
//static void test_double() { TestByType<double>(); } // not supported

/**
 * Test +=, -=, etc.
 * None of these are templated, so we might as well test them
 * separately instead of further complicating the templated tests.
 */
static void
TestHistogramBuilderArithmetic()
{
  float adata[] = { 0.1f, 2.9f, 3.1f, 14.0f }; // bins 0, 3, 3, 14
  float bdata[] = { 4.0f, 5.0f, 6.8f,  6.6f }; // bins 4, 5, 7, 7
  long long acount[] = {1,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0};
  long long bcount[] = {0,0,0,0,1,1,0,2,0,0,0,0,0,0,0,0};

  // All histograms are fixed size, 16 bins with width 1,
  // range 0..15 so bin 0 is centered on 0.0, 1 is centered
  // on 1.0 et cetera.
  HistogramBuilder a(16, 0.0, 15.0);
  HistogramBuilder b(16, 0.0, 15.0);
  HistogramBuilder c(16, 0.0, 15.0);
  a.add(&adata[0], &adata[4]);
  b.add(&bdata[0], &bdata[4]);
  c.add(&adata[0], &adata[4]);
  c.add(&bdata[0], &bdata[4]);

  // Iterating over index in our acount[], bcount[] arrays.
  // The sample value corresponding to each index happens to
  // be the index itself, since I chose origin=0, width=1.
  // Inside the loop, in acount[ii] ii is used as an index
  // while in a.get(ii) it is used as a value. Yes this could
  // been a bit clearer - sorry.
  for (int ii=0; ii<16; ++ii) {
    CPPUNIT_ASSERT_EQUAL(acount[ii], a.gethisto().get(ii));
    CPPUNIT_ASSERT_EQUAL(bcount[ii], b.gethisto().get(ii));
    CPPUNIT_ASSERT_EQUAL(acount[ii] + bcount[ii], c.gethisto().get(ii));
  }

  CPPUNIT_ASSERT_EQUAL(4LL, a.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(4LL, b.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(8LL, c.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(20.1, a.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(22.4, b.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.5, c.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, a.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, a.getstats().getmax(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 4.0, b.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 6.8, b.getstats().getmax(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, c.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, c.getstats().getmax(), 1e-6);

  // The statistics can also be extracted from the bins only,
  // but this is less accurate. The range will be that of the
  // occupied bins. plus 1/2 bin on either side since the bin
  // itself is only the center value.
  {
    StatisticData hs = c.gethiststats();
    CPPUNIT_ASSERT_EQUAL(8LL, hs.getcnt());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, hs.getmin(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(14.5, hs.getmax(), 1e-6);
  }

  // test operator== and operator!=
  CPPUNIT_ASSERT(a == a);
  CPPUNIT_ASSERT(a != b);
  CPPUNIT_ASSERT(a != c);
  CPPUNIT_ASSERT(c != a);

  // test those operators against an empty histogram
  HistogramBuilder empty(256, -1, +1), e2(16, -1, +1);
  CPPUNIT_ASSERT(empty == e2);
  CPPUNIT_ASSERT(empty != a);

  // test copy and assign
  HistogramBuilder h(c);
  CPPUNIT_ASSERT(h == c);
  CPPUNIT_ASSERT(h != a);
  h = a;
  CPPUNIT_ASSERT(h == a);
  CPPUNIT_ASSERT(h != c);

  // test addition
  h = a;
  h += b;
  CPPUNIT_ASSERT(h == c);
  // Comparison doesn't check the range...
  CPPUNIT_ASSERT_EQUAL(8LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, h.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, h.getstats().getmax(), 1e-6);

  // test subtraction
  h -= a;
  CPPUNIT_ASSERT(h == b);
  // The new min/max renge will not be exect, only approximated
  // by using the histogram. The true range is 4.0 - 6.8
  // this means bins 4 and 7 are occupied, and all we can tell
  // from the histogram alone is that the input must have been
  // in the range 3.5 to 7.5 since those extremes would have been
  // rounded to 4 anf 7 respectively.
  CPPUNIT_ASSERT_EQUAL(4LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.5, h.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 7.5, h.getstats().getmax(), 1e-6);

  h -= b;
  CPPUNIT_ASSERT(h == empty);
  CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());

  // test multiply with a constant
  h = c;
  h *= 4;
  CPPUNIT_ASSERT_EQUAL(4*c.getstats().getcnt(), h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, h.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, h.getstats().getmax(), 1e-6);
  for (int ii=0; ii<16; ++ii) {
    CPPUNIT_ASSERT_EQUAL(4*c.gethisto().get(ii), h.gethisto().get(ii));
  }

  h -= c;
  h -= c;
  h -= c;
  CPPUNIT_ASSERT(h == c);

  // Since the += operator allows the two histograms to have different
  // sampling, it can also be used to manually resize the histogram.
  // Some information may be lost if the old and new bins don't align
  // the same way.

  // Convert the histogram to one that is coarser, has fewer bins,
  // and doesn't even cover the same range as "c".
  // This histogram will have 4 bins holding values 0,1,2-3,4,5-6,7,8-9,10,11
  // if you consider integers only. And the assignment from "c" will do just
  // integers, since its bin-centers happen to fall on whole numbers 0,1,2,...
  // Looking at the input data, expect 1,4,2,0 entries plus one outside.

  HistogramBuilder coarse(4, 1.0, 10.0);
  coarse += c;
  CPPUNIT_ASSERT_EQUAL(8LL, coarse.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(1LL, coarse.gethisto().getbins()[0]);
  CPPUNIT_ASSERT_EQUAL(4LL, coarse.gethisto().getbins()[1]);
  CPPUNIT_ASSERT_EQUAL(2LL, coarse.gethisto().getbins()[2]);
  CPPUNIT_ASSERT_EQUAL(0LL, coarse.gethisto().getbins()[3]);
  // min/max range is a bit tricky. The lower range should be preserved.
  // For the upper range, the code ignores the nabove() count and can
  // from the histogram alone deduce that no value should be above
  // 11.5 which is the end of the last bin. But since the last bin is
  // unused, the max range is reduced with another bin width.
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, coarse.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 8.5, coarse.getstats().getmax(), 1e-6);

  // Convert the histogram to a larger one with both more samples
  // and a greater range.
  HistogramBuilder fine(2001, -100.0, +100.0);
  fine += c;

  for (int ii=0; ii<16; ++ii) {
    CPPUNIT_ASSERT_EQUAL(acount[ii] + bcount[ii], fine.gethisto().get(ii));
    CPPUNIT_ASSERT_EQUAL(0LL, fine.gethisto().get(ii-0.5));
    CPPUNIT_ASSERT_EQUAL(0LL, fine.gethisto().get(ii+0.5));
  }

  // Test scaling a histogram.
  // This is normally used when the histogram was build based on
  // integral data, and we ought to have applied a transform
  // on every single sample. But we prefer to do it later.
  HistogramBuilder scaled = c;

  // Scale in the other direction as the expected usage,
  // i.e. converting to a value range suitable to an int8.
  scaled.scale(0, 15, -128, 127);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(-126, scaled.getstats().getmin(), 1.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+110, scaled.getstats().getmax(), 1.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-128, scaled.gethisto().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+127, scaled.gethisto().getmax(), 1e-6);

  // Scale it back again, see that we get back the original.
  scaled.scale(-128, 127, 0, 15);
  CPPUNIT_ASSERT(scaled == c);

  scaled.scale(0, 15, 1, 16);
  // This scaling should add 1 to every sample - since we had 8 of them,
  // it is fairly straight forward to check everything.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.gethisto().getmin()+1, scaled.gethisto().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.gethisto().getmax()+1, scaled.gethisto().getmax(), 1e-6);
  for (int ii=0; ii<16; ++ii)
    CPPUNIT_ASSERT_EQUAL(c.gethisto().get(ii), scaled.gethisto().get(ii+1));
  CPPUNIT_ASSERT_EQUAL(scaled.getstats().getcnt(), scaled.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmin()+1, scaled.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmax()+1, scaled.getstats().getmax(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getsum()+8, scaled.getstats().getsum(), 1e-6);

  // Test constructing a histogram from a "foreign" type.
  // This is handled by using a constructor that gets all the
  // information it needs as discrete arguments. So we can also
  // use it to copy one of our own instances.
  const StatisticData& s = c.getstats();
  HistogramBuilder foreign(c.gethisto().getbins(), c.gethisto().getsize(), c.gethisto().getmin(), c.gethisto().getmax(),
                           s.getcnt(), s.getsum(), s.getssq(),
                           s.getmin(), s.getmax());
  CPPUNIT_ASSERT(foreign == c);
  // Comparison doesn't check the range...
  CPPUNIT_ASSERT_EQUAL(c.getstats().getcnt(), foreign.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmin(), foreign.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmax(), foreign.getstats().getmax(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getsum(), foreign.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getssq(), foreign.getstats().getssq(), 1e-6);

  // Check that infinites do not affect statistics except for being counted.
  float zero = 0.0;
  float trouble[2] = { 1.0f/zero, std::numeric_limits<float>::quiet_NaN() };
  HistogramBuilder withnan = c;
  withnan.add(&trouble[0], &trouble[2]);

  CPPUNIT_ASSERT_EQUAL(8LL, withnan.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(2LL, withnan.getstats().getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.5, withnan.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, withnan.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, withnan.getstats().getmax(), 1e-6);

  // The statistics can also be extracted from the bins only,
  // but this is less accurate. The range will be that of the
  // occupied bins. plus 1/2 bin on either side since the bin
  // itself is only the center value.
  // Also, the histogram no longer keeps track of infinite-count.
  {
    StatisticData hs = withnan.gethiststats();
    CPPUNIT_ASSERT_EQUAL(8LL, hs.getcnt());
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, hs.getmin(), 1e-6);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(14.5, hs.getmax(), 1e-6);
  }
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("histo.stats.dflt_ctor",   TestStatisticDataDefaultConstructor);
      register_test("histo.stats.copy_assign", TestStatisticDataCopyAssign);
      register_test("histo.stats.histo_ctor",  TestStatisticDataHistogramCtor);
      register_test("histo.stats.add",         TestStatisticDataAdd);
      register_test("histo.stats.arithmetic",  TestStatisticDataArithmetic);
      register_test("histo.normal_ctor",       TestHistogramBuilderNormalCtor);
      register_test("histo.discrete_ctor",     TestHistogramBuilderDiscreteCtor);
      register_test("histo.add_nan",           TestHistogramBulderAddNaN);
      register_test("histo.copy_assign",       TestHistogramBuilderCopyAssign);
      register_test("histo.compare",           TestHistogramBuilderCompare);
      register_test("histo.add",               TestHistogramBuilderAdd);
      register_test("histo.arithmetic",        TestHistogramBuilderArithmetic);
      register_test("histo.type.char",         test_char);
      register_test("histo.type.schar",        test_schar);
      register_test("histo.type.uchar",        test_uchar);
      register_test("histo.type.short",        test_short);
      register_test("histo.type.ushort",       test_ushort);
      register_test("histo.type.int",          test_int);
      register_test("histo.type.uint",         test_uint);
      register_test("histo.type.float",        test_float);
    }
  } dummy;
} // namespace for registration
