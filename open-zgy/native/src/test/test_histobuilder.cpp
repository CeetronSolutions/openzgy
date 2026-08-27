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

#define HAVE_EXPANDABLE_BUILDER 1
#if HAVE_EXPANDABLE_BUILDER
#include "../impl/expandablebuilder.h"
#endif

//using namespace OpenZGY;
using namespace InternalZGY;

#define CPPUNIT_ASSERT(a) TEST_CHECK((a))
#define CPPUNIT_ASSERT_EQUAL(a,b) TEST_EQUAL(b,a)
#define CPPUNIT_ASSERT_DOUBLES_EQUAL(a,b,eps) TEST_EQUAL_FLOAT(b,a,eps)


/*
 * Here is a tip for comparing the old ZGY test this is based on.
 *    Salmon/Zgy/UtilityLib/Test/HistogramBuilderTest.cpp
sed -e '
  s/  *$//
  s/HistogramBuilderTest:://
  s/^void /static void\n/
  s/^void$/static void/
  s/\bAdd\b/add/g
  s/GetBins/getbins/g
  s/GetSize/getsize/g
  s/\([.>]\)Scale\b/\1scale/g
  s/Bulder/Builder/g
  ' HistogramBuilderTest.cpp > /tmp/HistogramBuilderTest.cpp
some_diff_tool /tmp/HistogramBuilderTest.cpp test_histobuilder.cpp
 *
 */
namespace {
#if 0
}
#endif

#if 0 // for ad-hoc debugging
static std::string fmt(const HistogramBase& h)
{
  const double width = (h.getmax() - h.getmin()) / (h.getsize() - 1);
  std::stringstream ss;
  ss << "HistogramBase size = " << h.getsize()
     << " range " << h.getmin() << " to " << h.getmax()
     << " extreme " << h.getmin() - width/2 << " to " << h.getmax() + width/2
     << " width " << width;
  return ss.str();
}
#endif

static void
printBuilder(const char *msg, const HistogramData& h, const StatisticData& s, bool details)
{
  if (details)
    std::cout << msg << "stats " << s.toString() << " histo " << h.toString(true);
  else
    std::cout << msg << "stats " << s.toString() << " histo " << h.toString() << "\n";
}

static void
printBuilder(const char *msg, const HistogramBuilder& builder, bool details)
{
  printBuilder(msg, builder.gethisto(), builder.getstats(), details);
}

static double
getBinNumber(const HistogramData& h, double value)
{
  double A{0}, B{0};
  h.calculateConversionFactors(&A, &B);
  return (A + B*value);
}

static bool isAlmostIntegral(double value, double eps = 1e-6)
{
  return (std::fabs(value - std::floor(value)) < eps ||
          std::fabs(value - std::ceil(value)) < eps);
}

static bool isZeroCentric(const HistogramData& h)
{
  const double bin = getBinNumber(h, 0.0);
  return isAlmostIntegral(bin);
}

static bool isAntiZeroCentric(const HistogramData& h)
{
  const double bin = getBinNumber(h, 0.0);
  return isAlmostIntegral(bin + 0.5);
}

static bool isSymmetric(const HistogramData& h)
{
  return (std::fabs(h.getmin() + h.getmax()) /
          (std::fabs(h.getmin()) + std::fabs(h.getmax())))
    < 1e-5;
}

/**
 * Moved from HistogramData, because it is only needed by tests.
 *
 * Check all our bins, to see that the other histogram has the exact
 * same count when looking up the value that is the center point of
 * our own bins. Note that the other histogram could still have more
 * samples than we do, if it has a wider range or more densly sampled
 * bins. When comparing two histograms you probably want to compare
 * in both directions. Or at least do an extra check on the total
 * sample count - but I prefer the double check as it is guaranteed
 * to be symmetric.
 */
static bool
compare1(const HistogramData& self, const HistogramData& other)
{
  // special handling needed for empty histograms.
  // If histogram is not empty, we are guaranteed
  // that width is not zero. For empty one it might.

  HistogramData::count_type self_count = self.getcount();
  HistogramData::count_type other_count = other.getcount();

  if (self_count == 0 || other_count == 0)
    return self_count == other_count; // true only if both are empty.

  double width = (self.getmax() - self.getmin()) / (self.getsize() - 1);
  double begin = self.getmin();
  for (int ii=0; ii<self.getsize(); ++ii) {
    double value = begin + ii*width;
    if (self.getbins()[ii] != other.get(value))
      return false;
  }
  return true;
}

static bool
compare(const HistogramData& a, const HistogramData& b)
{
  return compare1(a, b) && compare1(b, a);
}

static bool
compare(const HistogramBuilder& a, const HistogramBuilder& b)
{
  return compare(a.gethisto(), b.gethisto());
}

static bool
checkStats(
     const HistogramBuilder& builder,
     std::int64_t cnt,
     double min, double max, double sum)
{
  bool ok{true};
  const StatisticData& s = builder.getstats();
  ok = TEST_EQUAL(s.getcnt(), cnt) && ok;
  ok = TEST_EQUAL_FLOAT(s.getmin(), min, 1e-5) && ok;
  ok = TEST_EQUAL_FLOAT(s.getmax(), max, 1e-5) && ok;
  ok = TEST_EQUAL_FLOAT(s.getsum(), sum, 1e-5) && ok;
  return ok;
}

static bool
checkHisto(
     const HistogramBuilder& builder,
     std::int64_t cnt,
     double min, double max, double binwidth)
{
  bool ok{true};
  const HistogramData& h = builder.gethisto();
  ok = TEST_EQUAL(h.getcount(), cnt) && ok;
  // Not testing the precise min/max limits here, because of the weird
  // slightly-not-symmetric code that I might be removing anyway.
  ok = TEST_EQUAL_FLOAT(h.getmin(), min, 0.5) && ok;
  ok = TEST_EQUAL_FLOAT(h.getmax(), max, 0.5) && ok;
  ok = TEST_EQUAL_FLOAT(h.getbinwidth(), binwidth, 1e-6) && ok;
  ok = TEST_CHECK(isAntiZeroCentric(h)) && ok;
  ok = TEST_CHECK(isSymmetric(h)) && ok;
  return ok;
}

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

/*
 * NOT PORTED FROM Salmon, the statistics-only builder has been removed.
static void
TestStatisticBuilderAddNaN()
template <typename T>
static void
TestStatisticBuilderFromRange()
 */

/**
 * Test the supported operators (+=, -=, *=) of StatisticsBuilder.
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

#if HAVE_EXPANDABLE_BUILDER // Only ExpandableBuilder has ctor without value range.
/**
 * Test the default constructor for the HistogramBuilder.
 * Does not depend on type, so only needs to be run once.
 */
static void
TestHistogramBuilderDefaultConstructor()
{
  const ExpandableBuilder h(256);

  CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getssq(), 1e-5);
  CPPUNIT_ASSERT_EQUAL(h.gethisto().getsize(), (int)std::count(h.gethisto().getbins(), h.gethisto().getbins() + h.gethisto().getsize(), 0));
  TEST_CHECK(isAntiZeroCentric(h.gethisto()));
  TEST_CHECK(isSymmetric(h.gethisto()));
}
#endif

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
TestHistogramBuilderAddNaN()
{
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  volatile float zero = 0.0;
  const float Inf = 1.0f / zero;

  static float data[3] = { NaN, Inf, -Inf };
  // Note: using invalid limits min=max=0 will now throw.
  // So, just choose something.
  HistogramBuilder h(256, -99, 99);
  h.add(&data[0], &data[0]); // empty range
  h.add(&data[0], &data[3]); // good range but no samples

  // Check the statistics part.
  CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_EQUAL(3LL, h.getstats().getinf());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getsum(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getssq(), 1e-5);

  // Check the histogram part - it should have no samples.
  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().getcount());

  CPPUNIT_ASSERT_EQUAL(0LL, h.gethisto().get(0));
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
 * The test is not very useful after operator== was removed by YAGNI.
 * Comparison is only done inside this unit test.
 *
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

  CPPUNIT_ASSERT(compare(ha, hb));
  CPPUNIT_ASSERT(compare(ha, hc));
  CPPUNIT_ASSERT(compare(ha, hd));

  float data[5] = { 12.1f, 20.0f, 20.1f, 21.0f, 50.0f };
  ha.add(&data[0], &data[5]);
  hb.add(&data[0], &data[5]);
  hc.add(&data[0], &data[5]);
  hd.add(&data[0], &data[5]);

  if (verbose()) {
    printBuilder("\nA: ", ha, true);
    printBuilder("B: ", hb, true);
    printBuilder("C: ", hc, true);
    printBuilder("D: ", hd, true);
  }

  CPPUNIT_ASSERT(compare(ha, hb));
  CPPUNIT_ASSERT(compare(ha, hc));
  CPPUNIT_ASSERT(!compare(ha, hd));

#if HAVE_EXPANDABLE_BUILDER
  // There is special case handling for empty histograms
  ExpandableBuilder empty1(12);
  ExpandableBuilder empty2(24);
  //illegal ExpandableBuilder empty3(48, 0, 100, true);
  CPPUNIT_ASSERT(compare(empty1, empty2));
  //CPPUNIT_ASSERT(compare(empty1, empty3));

  CPPUNIT_ASSERT(!(compare(ha, empty1)));
  CPPUNIT_ASSERT(!(compare(empty1, hb)));
#endif
}

/**
 * Test the add (single sample) and DoubleRange methods.
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

#if HAVE_EXPANDABLE_BUILDER

/**
 * Test the add (single sample) and DoubleRange methods.
 * There are currently not exposed publicly, so the test is incomplete.
 */
static void
TestHistogramBuilderAddExpand()
{
  // TODO-WIP-BrickedAPI: A non-symmetric ExpandableBuilder is illegal!
  ExpandableBuilder h(200, -100, 99, false);

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

  double data1[1] = { 105 };
  h.add(&data1[0], &data1[1]);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(-199, h.gethisto().getmin(), 1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(199, h.gethisto().getmax(), 1e-5);
  CPPUNIT_ASSERT_EQUAL(200, h.gethisto().getsize());

  // bin [100] range is now -0.5..+1.5, center value 0.5
  // bin [99]  has range -2.5..-0.5
  // bin [93]  has range -14.5 .. -12.5
  // So, where do the (-13) samples end up now? (answer: 93)
  //for (int ii=0; ii<200; ++ii)
  //  if (h.getbins()[ii])
  //    printf("Bin [%d]: %d\n", ii, int(h.getbins()[ii]));

  CPPUNIT_ASSERT_EQUAL(5LL, h.gethisto().getcount());
  CPPUNIT_ASSERT_EQUAL(2LL, h.gethisto().getbins()[93]);
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().getbins()[114]);
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().getbins()[149]);

  CPPUNIT_ASSERT_EQUAL(2LL, h.gethisto().get(-13));
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().get(29));
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().get(99));
  CPPUNIT_ASSERT_EQUAL(1LL, h.gethisto().get(105)); // now it is inside
}

/**
 * Test histograms with different value ranges and both fixed and dynamic ranges.
 * The same number of samples (3) will be added in each case, and only one of
 * the samples lie inside the initial value range.
 */
template <typename T>
static void
TestHistogramBuilderValueRange()
{
  // iterate fixed and dynamic histogram range
  for (size_t i = 0; i < 2; ++i) {

    // set fixed range indicator
    const bool fixed(i > 0);

    for (size_t j = 0; j < 5; ++j) {

      // create value range
      const double min(13.17 + 5.71*i);
      const double max(min + 19.47 + 2.23*j);

      // construct empty histogram
      // Note that we won't exercise the special shortcut for "char" here,
      // even for fixed histograms, because in that case min/max need
      // to be set to the limits for char.
      // TODO-WIP-BrickedAPI: A non-symmetric ExpandableBuilder is illegal!
      ExpandableBuilder h(256, min, max, fixed);

      // verify
      CPPUNIT_ASSERT_EQUAL(0LL, h.getstats().getcnt());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getsum(), 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0, h.getstats().getssq(), 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(min, h.gethisto().getmin(), 1e-5);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(max, h.gethisto().getmax(), 1e-5);
      CPPUNIT_ASSERT_EQUAL(h.gethisto().getsize(), (int)std::count(h.gethisto().getbins(), h.gethisto().getbins() + h.gethisto().getsize(), 0));

      // add some values
      const T value[3] = { static_cast<T>(min - 7), static_cast<T>((min + max)/2), static_cast<T>(max + 13) };
      h.add(&value[0], &value[0] + 3);

      // verify behaviour
      CPPUNIT_ASSERT_EQUAL(3LL, h.getstats().getcnt());
      CPPUNIT_ASSERT(h.getstats().getsum() != 0);
      CPPUNIT_ASSERT(h.getstats().getssq() > 0);
      if (fixed) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(min, h.gethisto().getmin(), 1e-5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(max, h.gethisto().getmax(), 1e-5);
        CPPUNIT_ASSERT_EQUAL(h.gethisto().getsize()-1, (int)std::count(h.gethisto().getbins(), h.gethisto().getbins() + h.gethisto().getsize(), 0)); // all but 1 bin empty
      }
      else {
        CPPUNIT_ASSERT(h.gethisto().getmin() <= (min - 7));
        CPPUNIT_ASSERT(h.gethisto().getmax() >= (max + 13));
        CPPUNIT_ASSERT_EQUAL(h.gethisto().getsize()-3, (int)std::count(h.gethisto().getbins(), h.gethisto().getbins() + h.gethisto().getsize(), 0)); // all but 3 bins empty
      }
    }
  }
}
#endif

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

/*
 * NOT PORTED FROM Salmon, stats/histo are now always in storage type.
 *
template <typename T>
static void
TestHistogramBuilderFromRangeWithConversion()
 */

/**
 * Invoke all the templated unit tests for a given template type.
 * This is a convenience just to avoid having to list N tests * M types
 * in the test suite.
 */
template <typename T>
static void
TestByType()
{
#if HAVE_EXPANDABLE_BUILDER
  TestHistogramBuilderValueRange<T>();
#endif
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
  CPPUNIT_ASSERT(compare(a, a));
  CPPUNIT_ASSERT(!compare(a, b));
  CPPUNIT_ASSERT(!compare(a, c));
  CPPUNIT_ASSERT(!compare(c, a));

  // test those operators against an empty histogram
  HistogramBuilder empty_builder(256, -1, +1), e2(16, -1, +1);
  CPPUNIT_ASSERT(compare(empty_builder, e2));
  CPPUNIT_ASSERT(!compare(empty_builder, a));

  // test copy and assign
  HistogramBuilder h(c);
  CPPUNIT_ASSERT(compare(h, c));
  CPPUNIT_ASSERT(!compare(h, a));
  h = a;
  CPPUNIT_ASSERT(compare(h, a));
  CPPUNIT_ASSERT(!compare(h, c));

  // test addition
  h = a;
  h += b;
  CPPUNIT_ASSERT(compare(h, c));
  // Comparison doesn't check the range...
  CPPUNIT_ASSERT_EQUAL(8LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, h.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(14.0, h.getstats().getmax(), 1e-6);

  // test subtraction
  h -= a;
  CPPUNIT_ASSERT(compare(h, b));
  // The new min/max range will not be exect, only approximated
  // by using the histogram. The true range is 4.0 - 6.8
  // this means bins 4 and 7 are occupied, and all we can tell
  // from the histogram alone is that the input must have been
  // in the range 3.5 to 7.5 since those extremes would have been
  // rounded to 4 and 7 respectively.
  // UPDATE: Narrow-to-histogram has been removed from the builder,
  // so the range will now be wider.
  CPPUNIT_ASSERT_EQUAL(4LL, h.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, h.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 14.0, h.getstats().getmax(), 1e-6);

  h -= b;
  CPPUNIT_ASSERT(compare(h, empty_builder));
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
  CPPUNIT_ASSERT(compare(h, c));

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
  // UPDATE: Narrow-to-histogram has been removed from the builder,
  // so the range will now be wider.
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.1, coarse.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( 14.0, coarse.getstats().getmax(), 1e-6);

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
  CPPUNIT_ASSERT(compare(scaled, c));

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
  CPPUNIT_ASSERT(compare(foreign, c));
  // Comparison doesn't check the range...
  CPPUNIT_ASSERT_EQUAL(c.getstats().getcnt(), foreign.getstats().getcnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmin(), foreign.getstats().getmin(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getmax(), foreign.getstats().getmax(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getsum(), foreign.getstats().getsum(), 1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(c.getstats().getssq(), foreign.getstats().getssq(), 1e-6);

  // Check that infinites do not affect statistics except for being counted.
  const float trouble[2] = { std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()};
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

static void
test_print()
{
  InternalZGY::HistogramBuilder builder(256, 0, 256);
  const std::vector<float> data {99, 56, 23, 101, 109, 5, 55};
  builder.add(data.data(), data.data() + data.size());
  const std::string dump = builder.gethisto().toString(/*verbose=*/true);
  if (verbose())
    std::cout << dump << std::flush;
}

static void
test_graph()
{
  InternalZGY::HistogramBuilder builder(256, 0, 256);
  const std::vector<float> data {99, 56, 23, 101, 109, 5, 55};
  builder.add(data.data(), data.data() + data.size());
  std::stringstream ss;
  builder.gethisto().toGraph(ss);
  const std::string dump = ss.str();
  if (verbose())
    std::cout << dump << std::flush;
}

#if HAVE_EXPANDABLE_BUILDER
/*
 * Test adding data in different order to a histogram.
 * This will affect the final histogram range and bin width,
 * but the result should of course still be correct.
 */
static void
TestHistogramOrderOfAddedData()
{
  // With dynamic range, the order of data being added may affect the final range.
  // There is a risk that my planned tweaks will change the behavior.
  // Test using the following data blocks:
  //
  //   noValues     - block entirely consisting of NaN, Inf.
  //   singleValue1 - some NaN and Inf values, some real values but those are all the same. (1)
  //   singleValue2 - as above but the constant value being outside the other data. (1)
  //   smallRange   - data in the range -1 .. 12
  //   largeRange   - data in the range -2 .. 80
  //
  // Here is a graphical description:
  //
  //                  -2  -1  0  1  2  12  80  100
  //   singleValue1 -           [.]
  //   singleValue2 -                          [.]
  //   smallRange   -     [..............]
  //   largeRange   - [......................]
  //
  // Adding the noValues block should never change anything, regardless of when it is done.
  // The other blocks may be added in any order, i.e. 4! or 24 combinations.
  // In every case the end result should be compatible, in the sense that if there are
  // enought bins then looking up each of the original values should return the same count.
  // Also the statistics should (of course) be the same.
  // The final histogram range is expected to vary qute a bit and isn't really important
  // as lomg as it is wider than the actual data range and not so wide that the bins
  // start getting crammed together. But I can add specific asserts just to keep track
  // of whether anything changed in the implementation.

  const float NaN = std::numeric_limits<float>::quiet_NaN();
  volatile float zero = 0.0;
  const float Inf = 1.0f / zero;

  const float data[5][10] = {
    { NaN, NaN, Inf, -Inf, NaN, NaN, NaN, NaN, NaN, NaN }, // noValues
    { NaN, Inf, 1.0, 1.0, NaN, 1.0, NaN, NaN, NaN, NaN }, // singleValue1 = 3 * 1.0
    { 100, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, // singleValue2 = 1 * 100.0
    {   8,  -1,  12,  12,  11, NaN, NaN, NaN, NaN, NaN }, // smallRange (5 values)
    {  -2,  80,  40, NaN, NaN, NaN, NaN, NaN, NaN, NaN }, // largeRange (3 values)
  };

  ExpandableBuilder allValues(256);
  allValues.add(&data[0][0], &data[5][0]);

  if (verbose()) {
    printBuilder("\nAll: ", allValues, true);
  }

  for (int ii=1; ii<=4; ++ii) {
    for (int jj=1; jj<=4; ++jj) {
      for (int xx=1; xx<=4; ++xx) {
        for (int yy=1; yy<=4; ++yy) {
          if (ii!=jj && ii!=xx && ii!=yy && jj!=xx && jj!=yy && xx!=yy) {
            ExpandableBuilder h(256);
            h.add(&data[0][0], &data[0][10]); // no-op
            h.add(&data[ii][0], &data[ii][10]);
            //printf("Starting with %d -> range %+g..%+g center %+g\n", ii, h.getmin(), h.getmax(), (h.getmin()+h.getmax())/2);
            h.add(&data[0][0], &data[0][10]); // no-op
            h.add(&data[jj][0], &data[jj][10]);
            h.add(&data[0][0], &data[0][10]); // no-op
            h.add(&data[xx][0], &data[xx][10]);
            h.add(&data[0][0], &data[0][10]); // no-op
            h.add(&data[yy][0], &data[yy][10]);

            const StatisticData& s = h.getstats();
            TEST_EQUAL(s.getcnt(), 12LL);
            TEST_EQUAL(s.getinf(), 68LL); // added 10*8 values, 12 of which were real
            TEST_EQUAL_FLOAT(s.getsum(),   263, 1e-5);
            TEST_EQUAL_FLOAT(s.getssq(), 18481, 1e-5);
            TEST_EQUAL_FLOAT(s.getmin(),    -2, 1e-5);
            TEST_EQUAL_FLOAT(s.getmax(),   100, 1e-5);

            // The number of bins don't change, only the histogram range.
            TEST_CHECK(!h.gethisto().isfixed());
            TEST_EQUAL(h.gethisto().getsize(), 256);

            // Total number of real added values is the same, regardless of ordering.
            ExpandableBuilder::count_type count = 0;
            for (int binno = 0; binno < h.gethisto().getsize(); ++binno)
              count += h.gethisto().getbins()[binno];
            TEST_EQUAL(count, 12LL);

            // TODO-WIP-BrickedAPI: The results depend on rounding.
            // E.g. adding sample 1.0 might end up in the bin for +0..+1,
            // but asking how many samples have the value 1.0 might
            // look in the bin for +1..+2 and this return zero.
            // With one scenario, might need to add 0.01 to all
            // negative values and subtract 0.01 from the positive.

            // Assume that we have enough bins, so each of out discrete values
            // ended up in separate bins. Thie means we can check each of them.
            TEST_EQUAL(h.gethisto().get(-2), 1LL);
            TEST_EQUAL(h.gethisto().get(-1), 1LL);
            TEST_EQUAL(h.gethisto().get(1), 3LL);
            TEST_EQUAL(h.gethisto().get(8), 1LL);
            TEST_EQUAL(h.gethisto().get(11), 1LL);
            TEST_EQUAL(h.gethisto().get(12), 2LL);
            TEST_EQUAL(h.gethisto().get(40), 1LL);
            TEST_EQUAL(h.gethisto().get(80), 1LL);
            TEST_EQUAL(h.gethisto().get(100), 1LL);

#if 0
            // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM is much simpler.
            // In practice, the first non-empty added block will determine
            // the final histogram value range. The numbers below have simply
            // been copied from a previous run, they have not been validated.
            switch (ii) {
              case 1: // 1.0 added first, causing range to be set to 1.0-eps, 1.0+eps
                TEST_EQUAL_FLOAT(h.gethisto().getmin(), -127, 0.1);
                TEST_EQUAL_FLOAT(h.gethisto().getmax(), 129, 0.1);
                break;
              case 2: // 100.0 added first, causing range to be set to 100.0-eps, 100.0+eps
                TEST_EQUAL_FLOAT(h.gethisto().getmin(), -28, 0.1);
                TEST_EQUAL_FLOAT(h.gethisto().getmax(), 228, 0.1);
                break;
              case 3: // -1..12 added first
                TEST_EQUAL_FLOAT(h.gethisto().getmin(), -98.5, 0.1);
                TEST_EQUAL_FLOAT(h.gethisto().getmax(), 109.5, 0.1);
                break;
              case 4: // -2..80 added first
                TEST_EQUAL_FLOAT(h.gethisto().getmin(), -43, 0.1);
                TEST_EQUAL_FLOAT(h.gethisto().getmax(), 121, 0.1);
                break;
              default:
                TEST_CHECK(false && "Invalid value of loop iterator");
            }
#endif
            // NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM is much simpler,
            // the final range is supposed to be the same for each ordering.
            // Results will depend on the choice of center- or extreme range.
            // And whether we have converted to zero-centric yet.
            TEST_EQUAL(h.gethisto().getcount(), h.getstats().getcnt());
            TEST_EQUAL_FLOAT(h.gethisto().getbinwidth(), 1.0, 1e-7);
            TEST_EQUAL_FLOAT(h.gethisto().getmin(), -127.5, 0.1);
            TEST_EQUAL_FLOAT(h.gethisto().getmax(), 127.5, 0.1);
            if (verbose()) {
              std::stringstream ss;
              ss << "Order: " << ii << " " << jj << " " << xx << " " << yy << " -> ";
              printBuilder(ss.str().c_str(), h, false);
            }
          }
        }
      }
    }
  }
}

/**
 * Calling ExpandableBuilder::operator+= with corner cases.
 *
 * - empty += full
 * - empty -= full
 * - full += empty
 * - full -= empty
 * - full -= full
 * - full -= copy(full)
 * - full += full
 * - full += copy(full)
 *
 * Note, when using oerator+= to change the number of bins,
 * it is advisable to only decrease the size, and only by an
 * integral factor. Otherwise the histogram will look bad.
 * And trimRange() will be slightly wrong because many couns
 * and up in the slightly wrong place.
 */
static void
TestHistogramBuilderWidenNot()
{
  const float data[4]{5,6,7,8};

  // Adding an empty builder to self is a no-op.
  // Ditto for subtract.
  {
    ExpandableBuilder a(256);
    ExpandableBuilder b(64);
    a.add(&data[0], &data[4]);
    a += b;
    TEST_CHECK(checkStats(a, 4, 5, 8, 26));
    TEST_CHECK(checkHisto(a, 4, -16, +16, 0.125));

    a -= b;
    TEST_CHECK(checkStats(a, 4, 5, 8, 26));
    TEST_CHECK(checkHisto(a, 4, -16, +16, 0.125));
  }

  // Adding a non-empty builder to empty self is an assign.
  // First the case where self doesn't have any range yet.
  {
    ExpandableBuilder a(256);
    ExpandableBuilder b(512); // Set to 64 to trigger a bug / "feature".
    b.add(&data[0], &data[4]);
    a += b;
    TEST_CHECK(checkStats(a, 4, 5, 8, 26));
    TEST_CHECK(checkHisto(a, 4, -16, +16, 0.125));
  }

  // Adding a non-empty builder to empty self is an assign.
  // The old range in "self" will be kept in the histogram.
  // trimRange() might be able to fix the statistics.
  // Use builders of the same size, to avoid going into the
  // who-cares-anyway cases.
  {
    ExpandableBuilder a(256);
    ExpandableBuilder b(256);
    b.add(&data[0], &data[4]);

    // Establish a too-wide range in the otherwise empty "a".
    const float widedata[1]{300};
    a.add(&widedata[0], &widedata[1]);
    ExpandableBuilder tmp(a);
    a -= tmp;
    TEST_CHECK(checkStats(a, 0, 300, 300, 0));
    // TODO-WIP-BrickedAPI: ??? The slight asymmetry is confusing.
    TEST_CHECK(checkHisto(a, 0, -510, +510, 4.0));

    // TODO-WIP-BrickedAPI: Depends on trimRange, which is affected by rounding.
    // The max() can be 8 (intuitive) or 12 (nunerically unstable).
    a += b;
    TEST_CHECK(checkStats(a, 4, 5, 300, 26)); // NOT trim changed 300 to 8.
    TEST_CHECK(checkHisto(a, 4, -510, +510, 4.0));
  }

  // Adding self to self (ref. equal) is a corner case.
  {
    ExpandableBuilder a(256);
    a.add(&data[0], &data[4]);
    a += a;
    TEST_CHECK(checkStats(a, 2*4, 5, 8, 2*26));
    TEST_CHECK(checkHisto(a, 2*4, -16, +16, 0.125));
  }

  // Adding self to a copy of self is not a corner case,
  // but it should give the same result.
  {
    ExpandableBuilder a(256);
    a.add(&data[0], &data[4]);
    a += ExpandableBuilder(a);
    TEST_CHECK(checkStats(a, 2*4, 5, 8, 2*26));
    TEST_CHECK(checkHisto(a, 2*4, -16, +16, 0.125));
  }

  // Subtracting self from self (ref. equal) is a corner case.
  // The builder should retain ranges.
  // The compiler might rightfully warn about "a -= a" because (a)
  // it is sily and probably unintentional, and (b) there is a fair
  // chance that whoever implemented the class failed to handle this
  // corner case. But I did handle it, and this is my unit test.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wself-assign-overloaded"
#endif
  {
    ExpandableBuilder a(256);
    a.add(&data[0], &data[4]);
    a -= a;
    TEST_CHECK(checkStats(a, 0, 5, 8, 0));
    TEST_CHECK(checkHisto(a, 0, -16, +16, 0.125));
    TEST_CHECK(!a.getstats().empty());
  }
#ifdef __clang__
#pragma clang diagnostic pop
#endif

  // Subtracting self from copy of self is not a corner case.
  // Should expect the same result as above.
  {
    ExpandableBuilder a(256);
    a.add(&data[0], &data[4]);
    a -= ExpandableBuilder(a);
    TEST_CHECK(checkStats(a, 0, 5, 8, 0));
    TEST_CHECK(checkHisto(a, 0, -16, +16, 0.125));
    TEST_CHECK(!a.getstats().empty());
  }

  // Subtracting something from empty builder
  {
    ExpandableBuilder a(256);
    ExpandableBuilder b(256);
    b.add(&data[0], &data[4]);
    a -= b;
    TEST_CHECK(checkStats(a, -4, 5, 8, -26));
    TEST_CHECK(checkHisto(a, -4, -16, +16, 0.125));

    a += b;
    TEST_CHECK(checkStats(a, 0, 5, 8, 0));
    TEST_CHECK(checkHisto(a, 0, -16, +16, 0.125));
  }
}

/**
 * Focusing only of widen range where it actually ends up being needed.
 */
static void
TestHistogramBuilderWidenNeeded()
{
  const float data[8]{5,6,7,8,100,200,300,400};
  {
    ExpandableBuilder a(512);
    a.add(&data[0], &data[4]);
    if (verbose())
      printBuilder("\n\nAdd 8:   ", a, true);
    TEST_CHECK(isAntiZeroCentric(a.gethisto()));
    TEST_CHECK(isSymmetric(a.gethisto()));
    TEST_EQUAL(a.gethisto().getcount(), 4);
    TEST_EQUAL_FLOAT(a.gethisto().getbinwidth(), 0.0625, 0.0001);
    a.add(&data[4], &data[8]);
    if (verbose())
      printBuilder("\nAdd 400: ", a, true);
    TEST_CHECK(isAntiZeroCentric(a.gethisto()));
    TEST_CHECK(isSymmetric(a.gethisto()));
    TEST_EQUAL(a.gethisto().getcount(), 8);
    TEST_EQUAL_FLOAT(a.gethisto().getbinwidth(), 2, 0.0001);
  }
}

/**
 * Floating point data that just happens to have the range -128..+127
 * should end up with that range and not double. This is a bit tricky
 * due to the slightly assymetric nature of a symmetric *and* zero cantric
 * histogram.
 */
static void
TestHistogramBuilderWiden128()
{
  const float data[5]{42, 0, -127, +127, -128};

  ExpandableBuilder a(1024);
  a.add(&data[0], &data[5]);
  if (verbose())
    printBuilder("\n", a, true);

  TEST_EQUAL(a.getstats().getcnt(), 5);
  TEST_EQUAL_FLOAT(a.getstats().getmin(), -128.0, 0.0001);
  TEST_EQUAL_FLOAT(a.getstats().getmax(), +127.0, 0.0001);

  TEST_CHECK(isAntiZeroCentric(a.gethisto()));
  TEST_CHECK(isSymmetric(a.gethisto()));
  TEST_EQUAL(a.gethisto().getcount(), 5);
  // The exact value here isn't important. It may depend on the
  // not-quite-symmetric setting, the rounding mode, and center vs. edge.
  TEST_EQUAL_FLOAT(a.gethisto().getmin(), -255.75, 0.0001);
  TEST_EQUAL_FLOAT(a.gethisto().getmax(), +255.75, 0.0001);
  TEST_EQUAL_FLOAT(a.gethisto().getbinwidth(), 0.5, 0.0001);

  HistogramData finalh = a.finalize(256);
  TEST_CHECK(isZeroCentric(finalh));
  TEST_EQUAL(finalh.getsize(), 256);
  TEST_EQUAL(finalh.getcount(), 5);
  TEST_EQUAL_FLOAT(finalh.getmin(), -128, 0.0001);
  TEST_EQUAL_FLOAT(finalh.getmax(), +127, 0.0001);
  TEST_EQUAL_FLOAT(finalh.getbinwidth(), 1.0, 0.0001);
}

/**
 * Focusing only of widen range when the existing range is a single value
 * Calling ExpandableBuilder::operator+=
 * NOTE: SYMMETRIC_EXPANDABLE_HISTOGRAM changes the expected behavior,
 * and the tests are not that relevant.
 */
static void
TestHistogramBuilderWidenSingle()
{
  const float data[4]{42.0f, 42.0f, 42.0f, 42.0f + 42.0f/32.0f};
  const bool symmetric = true;
  const double eps{pow(2, -124.0)};

  {
    // Empty statistics, but histogram must always be sensible.
    ExpandableBuilder a(1024);
    TEST_CHECK(a.getstats().empty());
    TEST_CHECK(double(a.gethisto().getmax()-a.gethisto().getmin()) >= eps);

    // Add some samples, all with the same value. Statistics should be exact.
    // Histogram wider then needed: For symmetrical histograms this will be
    // due to the 2^N requirement, for unconstrained histograms it will be
    // due to the slop value added to avoid division by zero.
    a.add(&data[0], &data[3]);
    TEST_CHECK(a.getstats().single());
    TEST_CHECK(double(a.gethisto().getmax()-a.gethisto().getmin()) >= eps);
    TEST_CHECK(a.gethisto().getmin() <=  42 - eps); // Actually much larger.
    TEST_CHECK(a.gethisto().getmax() >= +42 + eps); // in the symmetric case.

    a.add(&data[3], &data[4]);
    TEST_CHECK(!a.getstats().single());
    TEST_CHECK(double(a.gethisto().getmax()-a.gethisto().getmin()) >= eps);
    if (!symmetric) {
      // The slop value added around "42" to avoid divide by zero
      // should have been removed. I.e. the range should have
      // been slightly narrowed.
      TEST_EQUAL_FLOAT(a.getstats().getmin(), data[0], 1.0e-5);
      TEST_EQUAL_FLOAT(a.getstats().getmax(), data[3], 1.0e-5);
    }
  }

  {
    ExpandableBuilder a(1024);
    TEST_CHECK(a.getstats().empty());
    TEST_CHECK(double(a.gethisto().getmax()-a.gethisto().getmin()) >= eps);

    // Add some samples NOT with all the same values.
    // No slop needs to be added. For the symmetric case the range is still
    // extended due to the 2^N requirement.
    a.add(&data[0], &data[4]);
    TEST_CHECK(!a.getstats().single());
    TEST_CHECK(double(a.gethisto().getmax()-a.gethisto().getmin()) >= eps);
    if (!symmetric) {
      TEST_EQUAL_FLOAT(a.getstats().getmin(), data[0], 1.0e-5);
      TEST_EQUAL_FLOAT(a.getstats().getmax(), data[3], 1.0e-5);
    }
  }
}

/**
 * StatisticData::empty() is true if no finite values have been added.
 * There might still be NaN and +/- Inf. This means that adding an
 * "empty" builder isn't quite a no-op.
 */
static void
TestHistogramBuilderWidenNaN()
{
  constexpr float nan = std::numeric_limits<float>::quiet_NaN();
  constexpr float inf = std::numeric_limits<float>::infinity();
  constexpr float messy_data[4]{nan, nan, inf, -inf};
  constexpr float data[4]{5, 6, 7, 8};

  ExpandableBuilder empty_builder(256);
  TEST_CHECK(empty_builder.getstats().empty());
  TEST_EQUAL(empty_builder.getstats().getcnt(), 0);
  TEST_EQUAL(empty_builder.getstats().getinf(), 0);

  ExpandableBuilder messy_builder(256);
  messy_builder.add(&messy_data[0], &messy_data[4]);
  TEST_CHECK(messy_builder.getstats().empty());
  TEST_EQUAL(messy_builder.getstats().getcnt(), 0);
  TEST_EQUAL(messy_builder.getstats().getinf(), 4);

  ExpandableBuilder plain_builder(256);
  plain_builder.add(&data[0], &data[4]);
  TEST_CHECK(!plain_builder.getstats().empty());
  TEST_EQUAL(plain_builder.getstats().getcnt(), 4);
  TEST_EQUAL(plain_builder.getstats().getinf(), 0);

  {
    ExpandableBuilder a(empty_builder);
    // Should get the same result as an assignment.
    a += messy_builder;
    TEST_CHECK(a.getstats().empty());
    TEST_EQUAL(a.getstats().getcnt(), 0);
    TEST_EQUAL(a.getstats().getinf(), 4);
  }

  {
    ExpandableBuilder a(plain_builder);
    // Not a no-op, the inf count shoule get updated.
    // But do we really care? Possibly, if the builder
    // is used to test for the existence of nan.
    // Should get the same result as an assignment.
    a += messy_builder;
    TEST_CHECK(!a.getstats().empty());
    TEST_EQUAL(a.getstats().getcnt(), 4);
    TEST_EQUAL(a.getstats().getinf(), 4);
  }

  {
    ExpandableBuilder a(messy_builder);
    // Should keep a's inf count.
    a += empty_builder;
    TEST_CHECK(a.getstats().empty());
    TEST_EQUAL(a.getstats().getcnt(), 0);
    TEST_EQUAL(a.getstats().getinf(), 4);
  }

  {
    ExpandableBuilder a(messy_builder);
    // Should keep a's inf count.
    a += plain_builder;
    TEST_CHECK(!a.getstats().empty());
    TEST_EQUAL(a.getstats().getcnt(), 4);
    TEST_EQUAL(a.getstats().getinf(), 4);
  }
}

static void
TestHistogramBuilderFinalize()
{
  // TODO-WIP-BrickedAPI: Depends on rounding. There is no perfect solution.
  if (verbose())
    std::cout << std::endl;

  std::vector<float> data(4567);
  for (std::int64_t ii = 0; ii < (std::int64_t)data.size(); ++ii)
    data[ii] = (float)ii - 401.0f;

  ExpandableBuilder builder(4096);
  builder.add(data.begin(), data.end());
  const HistogramData& h = builder.gethisto();
  const StatisticData& s = builder.getstats();
  // sum of an arithmetic sequence is (begin+end)/2 * count
  TEST_EQUAL_FLOAT(s.getmin(), -401.0, 1e-6);
  TEST_EQUAL_FLOAT(s.getmax(), 4165.0, 1e-6);
  TEST_EQUAL_FLOAT(((-401.0 + 4165.0) / 2) * 4567, s.getsum(), 0.01);
  // Expect range +/- 8190  (8192 edge-to-edge).
  // With most bins containing 4 samples.
  //printBuilder("\nAdded: ", builder, true);
  TEST_EQUAL_FLOAT(h.getbinwidth(), 4, 1e-7);
  // Need a generous epsilon if "noise" of 2^-30 * nbins have been added.
  // If instead the noise is set to 2^-50 this is not needed.
  TEST_EQUAL_FLOAT(h.getmin(), -8190, 1e-5);
  TEST_EQUAL_FLOAT(h.getmax(), +8190, 1e-5);
  // Make sure none of the added samples got clipped away
  TEST_EQUAL(h.getcount(), 4567);
  // All negative sample values have a small increase in factor,
  // so the bin edges get slightly more negative,
  // meaning integrals tend to be rounded upwards towards zero.
  // The same applies to integer value zero. There is a tiny
  // negative offset which means the bin number gets rounded up.
  // All positive values get rounded down because of the increased factor.
  // Bin 0 has center value -8190, range [-8192,-8188>
  // Bin 1947 add 4*1947, center -402, range [-404..-400>
  //   => my samples -401, (lowest) ends up there.
  // Bin 1948 add 4*1948, center -398, range [-400..-396>
  // Bin 2047 add 4*2047, center -2, range [-4..0>
  // Bin 2048 add 4*2048, center +2, range [0..+4] (note rounding!)
  // Bin 2049 add 4*2049, center +6, range <4..+8] (now rounds DOWN)
  // Bin 3089 add 4*3089, center 4166, range <4164..4168]
  //   => my sample 4165 (last) end up there.

  TEST_EQUAL(h.getbins()[1946], 0LL);
  TEST_EQUAL(h.getbins()[1947], 1LL);
  TEST_EQUAL(h.getbins()[1948], 4LL);
  TEST_EQUAL(h.getbins()[1949], 4LL);
  TEST_EQUAL(h.getbins()[1950], 4LL);

  // Because of rounding towards zero, the center bin might get an extra sample.
  TEST_EQUAL(h.getbins()[2047], 4LL);
  TEST_EQUAL(h.getbins()[2048], 5LL);
  TEST_EQUAL(h.getbins()[2049], 4LL);
  TEST_EQUAL(h.getbins()[3087], 4LL);
  TEST_EQUAL(h.getbins()[3088], 4LL);
  TEST_EQUAL(h.getbins()[3089], 1LL);
  TEST_EQUAL(h.getbins()[3090], 0LL);
  TEST_EQUAL(h.getbins()[3091], 0LL);
  // Looking at my calculations done by hand, 1948..3088 inclusive, but
  // excluding bin 2048, have 4560 samples total. Add 5 in bin 2048, 1
  // in 1947, 1 in 3089 for a total of 4567. Which is the correct sum.

  //std::cout << "\n";
  //for (int ii = -410; ii <= -390; ++ii)
  //  std::cout << "Real " << ii << " ==> " << h.get(ii) << "\n";
  //std::cout << std::endl;

  HistogramData finalh = builder.finalize(256);
  TEST_CHECK(isAntiZeroCentric(h));
  TEST_CHECK(isSymmetric(h));
  TEST_CHECK(isZeroCentric(finalh));
  if (verbose()) {
    printBuilder("\n", h, builder.getstats(), false);
  }
}

static void
TestHistogramBuilderOnlyZero()
{
  std::vector<float> data{0,0,0};
  ExpandableBuilder builder(512);
  builder.add(data.begin(), data.end());

  TEST_CHECK(checkStats(builder, 3, 0, 0, 0));
  // Approximate, but epsilon in compare handles it.
  TEST_CHECK(checkHisto(builder, 3, -1e-28, +1e-28, 1e-30));

  HistogramData finalh = builder.finalize(256);
  TEST_CHECK(isZeroCentric(finalh));
  TEST_CHECK(finalh.getbinwidth() > 0);
  // By convention, the single value should appear in the center bin +/- 1.
  TEST_EQUAL(finalh.getbins()[finalh.getsize()/2-1] +
             finalh.getbins()[finalh.getsize()/2] +
             finalh.getbins()[finalh.getsize()/2+1], 3LL);
  if (verbose()) {
    printBuilder("\ninternal: ", builder, true);
    printBuilder("\nfinished: ", finalh, builder.getstats(), true);
  }
}

static void
TestHistogramBuilderOnlyOne()
{
  std::vector<float> data{42,42,42};
  ExpandableBuilder builder(512);
  builder.add(data.begin(), data.end());

  TEST_CHECK(checkStats(builder, 3, 42, 42, 3*42));
  // Approximate, depends on not-quite-symmetric.
  TEST_CHECK(checkHisto(builder, 3, -64.125, +63.635, 0.25));

  HistogramData finalh = builder.finalize(256);

  // TODO-WIP-BrickedAPI: Questionable test.
  // zero-centric is not guaranteed when 0 is not included in range.
  //TEST_CHECK(isZeroCentric(finalh));
  // The actual bin width is not relevant, and how it is chosen is intricate.
  TEST_EQUAL_FLOAT(finalh.getbinwidth(), 0.5, 0.001);
  // By convention, the single value should appear in the center bin +/- 1.
  TEST_EQUAL(finalh.getbins()[finalh.getsize()/2-1] +
             finalh.getbins()[finalh.getsize()/2] +
             finalh.getbins()[finalh.getsize()/2+1], 3LL);
  if (verbose()) {
    printBuilder("\ninternal: ", builder, true);
    printBuilder("\nfinished: ", finalh, builder.getstats(), true);
  }
}

/**
 * Create a builder with no samples, having the size, limits,
 * and expandable status copied from an existing lower level
 * HistogramData.
 */
static void
TestHistogramBuilderClone()
{
  for (int ii=0; ii<2; ++ii) {
    const bool fixed = (ii==0);
    auto builder = fixed ?
      std::make_shared<ExpandableBuilder>(64, -111, +222) :
      std::make_shared<ExpandableBuilder>(64);
    std::vector<float> data{-99, 42, 16};
    builder->add(data.begin(), data.end());
    HistogramData histodata = builder->gethisto();
    ExpandableBuilder clone(histodata, StatisticData(), /*copy=*/false);
    if (verbose()) {
      printBuilder(fixed ? "\norig   fix: " : "\norig   dyn: ", *builder, true);
      printBuilder(fixed ? "\ncloned fix: " : "\ncloned dyn: ", clone, true);
    }

    TEST_EQUAL(builder->gethisto().isfixed(), fixed);
    TEST_EQUAL(builder->gethisto().getcount(), 3);
    TEST_EQUAL(builder->getstats().getcnt(), 3);
    TEST_EQUAL(clone.gethisto().isfixed(), fixed);
    TEST_EQUAL(clone.gethisto().getcount(), 0);
    TEST_EQUAL(clone.getstats().getcnt(), 0);
    TEST_EQUAL(clone.gethisto().getmin(), builder->gethisto().getmin());
    TEST_EQUAL(clone.gethisto().getmax(), builder->gethisto().getmax());
  }
}

/**
 * Example use of the histogram builder.
 *
 * Collect sample statistics and create a 4096-bin intermediate histogram.
 * When done collecting, convert the histogram to 256-bin zero-centric.
 * The intermediate histogram can also be accessed. But that will have many
 * empty bins and won't be zero-centric.
 *
 * The example also shows how to create a fixed-width histogram. The caller
 * is responsible for setting the limits so the result ends up zero-centric.
 *
 * When the builder is not fixed-width, the intermediate size a.k.a. bin count
 * must be at least twice that of the desired result. Choosing a too large
 * size might have a performance impact. Especially on the cost of operator+=.
 *
 * When the builder is fixed-width, the intermediate and final sizes
 * must be the same.
 *
 * This specific example calls add() with an iterator from a vector<float>.
 * Other forward iterators should work as well. Ditto for int8 and int16 data.
 *
 * For multi-threaded computation, use one builder for each thread and combine
 * the results at the end by using operator+=.
 */
static void
dynamicExample()
{
  const std::vector<float> data1 {-99, 56, 23};
  const std::vector<float> data2 {42};

  ExpandableBuilder builder(4096);              // Expands to fit the data.
  //ExpandableBuilder builder(256, -100, +100); // Fixed histogram limits.
  builder.add(data1.begin(), data1.end());
  builder.add(data2.begin(), data2.end());

  const StatisticData& stats = builder.getstats();
  const HistogramData& histo = builder.finalize(256);
  const HistogramData& intermediate = builder.gethisto();

  if (verbose()) {
    printBuilder("\n", intermediate, stats, false);
    printBuilder("", histo, stats, true);
  }
}

/**
 * Example of generating a fixed-width histogram, limits given by caller.
 *
 * You might also use ExpandableBuilder for this. See dynamicExample().
 * It is cleaner to use ExpandableBuilder only when really needed, and
 * the HistogramBuilder base class otherwise.
 *
 * You probably want to choose the limits carefully to make the resulting
 * histogram zero-centric.
 *
 * There is a caveat when mixing HistogramBuilder and ExpandableBuilder
 * in the same code. Do not access an ExpandableBuilder thru a pointer
 * or reference to an ExpandableBuilder. The add() method exists in
 * both classes and is not virtual. So the automatic widening won't work.
 */
static void
staticExample()
{
  const std::vector<float> data1 {-99, 56, 23};
  const std::vector<float> data2 {42};

  HistogramBuilder builder(256, -100, +100);
  builder.add(data1.begin(), data1.end());
  builder.add(data2.begin(), data2.end());

  const StatisticData& stats = builder.getstats();
  const HistogramData& histo = builder.gethisto();

  if (verbose()) {
    printBuilder("\n", histo, stats, true);
  }
}

#endif

/*
   More tests?

*) The high level add() method may in some cases bypass the "proper" add() methods
   in both the histogram and the statistics builder. Try calling both the low level
   methods and the optimized high level add() and verify that the result is the same.
   Run this test with fixed size histograms only, as the "widen" logic might only
   be available in the high level method. But be sure to run the test for all value types.
   Note that this test might not be possible until the classes have been refactored.

*) The high level add() method may invoke completely different algorithms depending on
   the value type of the data being added. So make sure the relevant tests are executed
   for all types. (May be Ok with this already).

   The high level code also makes some assumptions about the usage of the classes.
   But this is checked wherever possible.
   Specifically: When add() is being called with integral values, the histogram range
   is assumed to be identical to the coding range. If no coding range or scaling is
   given, this means that the histogram range should match the value range of the
   data being passed. A corrollary to this is that while data may be added in multiple
   chunks, the data will always be added in the same value type that was used previously.
*/

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
      register_test("histo.add_nan",           TestHistogramBuilderAddNaN);
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
      register_test("histo.print",             test_print);
      register_test("histo.graph",             test_graph);
#if HAVE_EXPANDABLE_BUILDER
      register_test("histoex.00.dflt_ctor",       TestHistogramBuilderDefaultConstructor);
      register_test("histoex.01.add_expand",      TestHistogramBuilderAddExpand);
      register_test("histoex.02.order",           TestHistogramOrderOfAddedData);
      register_test("histoex.03.widen_not",       TestHistogramBuilderWidenNot);
      register_test("histoex.04.widen_needed",    TestHistogramBuilderWidenNeeded);
      register_test("histoex.05.widen_128",       TestHistogramBuilderWiden128);
      register_test("histoex.06.widen_single",    TestHistogramBuilderWidenSingle);
      register_test("histoex.07.widen_nan",       TestHistogramBuilderWidenNaN);
      register_test("histoex.08.finalize",        TestHistogramBuilderFinalize);
      register_test("histoex.09.onlyzero",        TestHistogramBuilderOnlyZero);
      register_test("histoex.10.onlyone",         TestHistogramBuilderOnlyOne);
      register_test("histoex.11.clone",           TestHistogramBuilderClone);
      register_test("histoex.static_example",  staticExample);
      register_test("histoex.dynamic_example", dynamicExample);
#endif
    }
  } dummy;
} // namespace for registration
