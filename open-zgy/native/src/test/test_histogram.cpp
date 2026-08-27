// Copyright 2017-2023, Schlumberger
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
#include "../impl/environment.h"
#include "../api.h"
#include "../writerargs.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#define HAVE_EXPANDABLE_BUILDER 1
#if HAVE_EXPANDABLE_BUILDER

using namespace OpenZGY;

using Test_Utils::LocalFileAutoDelete;
using Test_Utils::LocalFileNoDelete;

#define CPPUNIT_ASSERT(a) TEST_CHECK((a))
#define CPPUNIT_ASSERT_EQUAL(a,b) TEST_EQUAL(b,a)
#define CPPUNIT_ASSERT_DOUBLES_EQUAL(a,b,eps) TEST_EQUAL_FLOAT(b,a,eps)

typedef std::array<std::int64_t,3> size3i_t;

namespace {
#if 0
}
#endif

/////////////////////////////////////////////////////////////////////////////
//
// Important differences from the old ZGY accessor.
//
//  - OpenZGY does not support dead traces. A whole test function has
//    been removed, and some tweaks are needed other places.
//
//  - OpenZGY alway delivers unconverted values when reading into an
//    integral buffer. Old ZGY allowed the application to request an
//    arbitrary conversion. E.g. do the regular scaling to original
//    float data and then cast the values to integer. This is silly.
//    Perhaps less so in old ZGY, where int8 data could be read as
//    int16 or vice versa. But that feature itself was also pretty
//    useless. Unfortunately the unit test from old ZGY relies on that
//    silly feature.
//
//    TODO-WIP-BrickedAPI: Didn't I change this?
//  - OpenZGY will consistently throw an error if trying to read or
//    write outside the last brick. Access between the survey edge and
//    any padding is allowed. In the old accessor the behavior is less
//    consistent. The general API will ignore any access outside the
//    survey area. On the other hand, reads even from inside the
//    survey (and possibly outside) might return BrickNotFound if no
//    bricks at all had been written in that area. If only some bricks
//    were found then there will be no error and the missing data is
//    assumed to be zero. With the bug that no conversion is done on
//    those values. TODO: OpenZGY should possibly relax its checking.
//
//  - OpenZGY includes all finite sample values in the statistics and
//    histogram. The old library was inconsistent with what it included.
//
//  - OpenZGY does not allow the application to specify a minimum
//    range for a histogram of float samples. In this test suite, that
//    makes a big difference in testUpdate(). That test assumed the
//    histogram has fixed limits even for floats.
//    TODO-WIP-BrickedAPI: I fixed this now.
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
///   SUPPORTING CODE   /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

namespace {

  // The old ZGY simply ignored any access outside the survey.
  // Some of the tests in this file relied on that. Provide
  // saferead() and safewrite() functions to skip those parts.
  // The functionality is VERY LIMITED, allowing only those
  // cases that the unit tests need. In particular, access
  // past the end is supported. Not access before the start.

  static size3i_t
  clip(const size3i_t& beg, const size3i_t& size, const size3i_t& survey)
  {
    size3i_t result;
    for (int dim = 0; dim < 3; ++dim)
      result[dim] =
        std::max(std::int64_t(0),
                 std::min(beg[dim] + size[dim], survey[dim]) - beg[dim]);
    return result;
  }

  void safecheck(
       const size3i_t& begin,
       const size3i_t& size,
       const size3i_t& newsize, int lod)
  {
    if (lod != 0)
      throw std::runtime_error("lod != 0 not implemented in this test.");
    if (begin[0] < 0 || begin[1] < 0 || begin[2] < 0)
      throw std::runtime_error("negative start not implemented in this test.");
    if (newsize[1] != size[1] || newsize[2] != size[2])
      throw std::runtime_error("Oops, only dim 0 handled for now.");
  }

  template<typename T>
  void
  safewrite(
       const std::shared_ptr<IZgyWriter>& w,
       const size3i_t& begin,
       const size3i_t& size,
       const T* data)
  {
    const size3i_t newsize = clip(begin, size, w->size());
    safecheck(begin, size, newsize, 0);
    if (newsize[0] > 0 && newsize[1] > 0 && newsize[2] > 0)
      w->write(begin, newsize, data);
  }

  template<typename T>
  void
  saferead(
       const std::shared_ptr<IZgyWriter>& w,
       const size3i_t& begin,
       const size3i_t& size,
       T* data)
  {
    const size3i_t newsize = clip(begin, size, w->size());
    safecheck(begin, size, newsize, 0);
    if (newsize != size)
      std::fill(data, data + size[0] + size[1] + size[2], static_cast<T>(0));
    if (newsize[0] > 0 && newsize[1] > 0 && newsize[2] > 0)
      w->read(begin, newsize, data);
  }

  template<typename T>
  void
  saferead(
       const std::shared_ptr<IZgyReader>& r,
       const size3i_t& begin,
       const size3i_t& size,
       T* data,
       int lod)
  {
    const size3i_t newsize = clip(begin, size, r->size());
    safecheck(begin, size, newsize, 0);
    if (newsize != size)
      std::fill(data, data + size[0] + size[1] + size[2], 0);
    if (newsize[0] > 0 && newsize[1] > 0 && newsize[2] > 0)
      r->read(begin, newsize, data, lod);
  }
}

/**
 * Hold the satistics information so we don't need to pass the discrete
 * values around the test code. Only used internally.
 * count is treated as a double, which handles 52 bits worth of integers
 * without problems. That is sufficient, but admittedly a bit odd since
 * the count is really an integer and should be modeled that way.
 * By convention, when count is zero all other params return as zero.
 */
class MyStatistics
{
public:
  std::int64_t cnt() const { return cnt_; }
  double      sum() const { return sum_; }
  double      ssq() const { return ssq_; }
  double      min() const { return cnt_ ? min_ : 0; }
  double      max() const { return cnt_ ? max_ : 0; }
  double      avg() const { return cnt_ ? sum_ / cnt_ : 0; }
  double      std() const { return cnt_ ? std::sqrt(ssq_/cnt_ - avg()*avg()) : 0; }

  MyStatistics()
    : cnt_(0), sum_(0), ssq_(0), min_(0), max_(0) { }

  MyStatistics(std::int64_t cnt, double sum, double ssq, double min, double max)
    : cnt_(cnt), sum_(sum), ssq_(ssq), min_(min), max_(max) { }

  // Implicit converters
  MyStatistics(const OpenZGY::SampleStatistics& ss) :
    cnt_(ss.cnt), sum_(ss.sum), ssq_(ss.ssq), min_(ss.min), max_(ss.max) { }

  std::string toString() const
  {
    std::stringstream ss;
    ss << std::setw(10) << min() << " to " << std::setw(10) << max()
       << " avg " << std::setw(10) << avg()
       << " std " << std::setw(10) << std()
       << " cnt " << std::setw(8) << cnt();
    return ss.str();
  }

#if 0
  bool equals(const MyStatistics& other, double eps) const
  {
    return
      cnt() == other.cnt() &&
      std::fabs(avg() - other.avg()) <= eps &&
      std::fabs(std() - other.std()) <= eps &&
      std::fabs(min() - other.min()) <= eps &&
      std::fabs(max() - other.max()) <= eps;
  }
#endif

  bool asExpected(const MyStatistics& expect, double eps) const
  {
    return (TEST_EQUAL(cnt(), expect.cnt()) &&
            TEST_EQUAL_FLOAT(avg(), expect.avg(), eps) &&
            TEST_EQUAL_FLOAT(std(), expect.std(), eps) &&
            TEST_EQUAL_FLOAT(min(), expect.min(), eps) &&
            TEST_EQUAL_FLOAT(max(), expect.max(), eps));
  }

private:
  std::int64_t cnt_;
  double sum_;
  double ssq_;
  double min_;
  double max_;
};

class MyHistogram : public OpenZGY::SampleHistogram
{
public:
  MyHistogram(const OpenZGY::SampleHistogram& other)
    : SampleHistogram(other)
  {
  }

  MyHistogram()
  {
  }

  std::string toString(int details = 1) const
  {
    const float width = (float)
      (bins.size() <= 1 ? 0.0 :
      float((maxvalue - minvalue) / (bins.size() - 1)));
    const float begin = (float)minvalue;

    std::stringstream ss;
    ss << std::setw(10) << minvalue << " to " << std::setw(10) << maxvalue
       << " cnt " << std::setw(8) << samplecount
       << " bins [";

    if (minvalue >= maxvalue && samplecount == 0) {
      ss << " (empty histogram) ]";
      return ss.str();
    }

    switch (details) {
    default:
    case 0:
      ss << " count = " << bins.size();
      break;

    case 1:
      for (std::size_t ii=0; ii<bins.size(); ++ii)
        ss << " " << bins[ii];
      break;

    case 2:
    case 3:
      for (std::size_t ii=0; ii<bins.size(); ++ii) {
        const float center = begin+ii*width;
        const float lo = (float)(begin+(ii-0.5)*width);
        const float hi = (float)(begin+(ii+0.5)*width);
        if (bins[ii] != 0 || (lo<=0 && hi>=0) || details > 2)
          ss << "\n         bin[" << std::setw(3) << ii << "]: "
             << std::setw(8) << bins[ii]
             << " * " << std::setw(8) << center
             << " (" << lo << ".." << hi << ")";
      }
    }
    ss << " ]";
    return ss.str();
  }

#if 0
  bool equals(const OpenZGY::SampleHistogram& other, double eps)
  {
    if (samplecount != other.samplecount ||
        std::fabs(minvalue - other.minvalue) > eps ||
        std::fabs(maxvalue - other.maxvalue) > eps)
      return false;
    for (std::size_t ii=0; ii<bins.size(); ++ii)
      if (bins[ii] != other.bins[ii])
        return false;
    return true;
  }
#endif

  bool asExpected(const OpenZGY::SampleHistogram& expected, double eps)
  {
    if (!TEST_EQUAL(samplecount, expected.samplecount) ||
        !TEST_EQUAL_FLOAT(minvalue, expected.minvalue, eps) ||
        !TEST_EQUAL_FLOAT(maxvalue, expected.maxvalue, eps))
      return false;
    for (std::size_t ii=0; ii<bins.size(); ++ii)
      if (bins[ii] != expected.bins[ii])
        if (!TEST_EQUAL(bins[ii], expected.bins[ii]))
          return false;
    return true;
  }

  MyStatistics stats() const
  {
    const double begin = minvalue;
    const double width = (maxvalue - minvalue) / (bins.size() - 1);

    double hcnt = 0, hsum = 0, hssq = 0, hmin = 0, hmax = 0;
    for (std::size_t binno = 0; binno < bins.size(); ++binno) {
      double count = (double)bins[binno];
      double value = begin + binno*width;
      if (count != 0) {
        if (hcnt == 0)
          hmin = value - width/2; // first non-empty bin found
        hmax = value + width/2; // last (so far) non-empty bin found.
      }
      hcnt += count;
      hsum += count * value;
      hssq += count * value * value;
    }
    // Count is "double", just to save typing.
    return MyStatistics((std::int64_t)hcnt, hsum, hssq, hmin, hmax);
  }

  int getindex(double value, bool neighbor = false)
  {
    double b = (bins.size() - 1)/(maxvalue - minvalue);
    double a = -minvalue*b;
    double index = a + b*value;
    int n = static_cast<int>(std::floor(index+0.5));
    if (neighbor) {
      if (n < index) // we rounded down
        ++n; // so, round up instead.
      else
        --n; // round down instead.
    }
    return n;
  }

  /**
   * Return the number of samples found close to this value.
   * If value is close to the boundary between two cells,
   * you might not get the bin you expected due to numerical
   * inaccuracy.
   */
  std::int64_t get(double value)
  {
    int n = getindex(value);
    return n >= 0 && (std::size_t)n < bins.size() ? bins[n] : 0;
  }

  /**
   * Return the number of samples found close to this value,
   * as the sum of the two bins closest to the value.
   * In some very specific unit tests where the histogram
   * is expected to be sparse, this can work around the
   * rounding issues.
   */
  std::int64_t getsloppy(double value)
  {
    std::int64_t result = 0;
    int n = getindex(value);
    if (n >= 0 && (std::size_t)n < bins.size())
      result += bins[n];
    n = getindex(value, true);
    if (n >= 0 && (std::size_t)n < bins.size())
      result += bins[n];
    return result;
  }
};

/**
 * Return a single sample from our synthetic volume, defined as follows:
 * The diagonal having i==j==k running across the entire cube will have
 * sample = positive or negative infinity.
 * The brick with coordinates i==0,j==1,k==0 is completely filled with NaN.
 * The brick with coordinates i==0,j==2,k==0 is completely filled with 0.5
 * All traces with j==200 will be marked as dead. This function doesn't
 * care and will return the usual test data, but the writer should mark
 * those bits in the alpha plane as unused.
 * All other sample values are real numbers in the range 0..100.
 * 0 at the cube's origin and 100 at the opposite corner.
 * This assumes a hard coded size of (117, 241, 333). If the size
 * is different, the result will not have the same range.
 * Pass isfloat=false if either the application or the storage valuetype
 * is integral. In that case we won't return any Infinite or NaN values.
 */
double
getSyntheticSample(std::int64_t i, std::int64_t j, std::int64_t k, bool isfloat)
{
  static const std::int64_t szi = 117;
  static const std::int64_t szj = 241;
  static const std::int64_t szk = 333;

  // FAILS: There are some subtle issues when zero is not part of the
  // original data. It may get added indirectly to the statistics when
  // writing partial bricks, and will not be removed again. This is
  // a "feature". To avoid failing the unit tests, I will simply ensure
  // thet zero is found somewhere...
  if (i==0 && j==0 && k==0)
    return 0;

  if (isfloat) {
    if (i==j && i==k && i%2 == 0)
      return std::numeric_limits<float>::infinity();
    if (i==j && i==k && i%2 == 1)
      return -std::numeric_limits<float>::infinity();
    if (i<64 && j>=64 && j<128 && k >= 0 && k < 64)
      return std::numeric_limits<float>::quiet_NaN();
  }
  if (i<64 && j>=128 && j<192 && k >= 0 && k < 64)
    return 0.5;
  double ri = i*1.0/(szi - 1);
  double rj = j*1.0/(szj - 1);
  double rk = k*1.0/(szk - 1);
  return 100.0 * sqrt((ri*ri + rj*rj + rk*rk)/3.0);
}

/*
 * Get synthetic data for the required brick, ordered i slowest, k fastest.
 * The cube coordinates of the start of the brick are passed in because
 * the actual samples we fill with depend on them.
 */
template <typename T>
void
getSyntheticBrick(T* result, std::int64_t oi, std::int64_t oj, std::int64_t ok, std::int64_t si, std::int64_t sj, std::int64_t sk, bool isfloat)
{
  for (std::int64_t ii=0; ii<si; ++ii)
    for (std::int64_t jj=0; jj<sj; ++jj)
      for (std::int64_t kk=0; kk<sk; ++kk)
        result[ii*sj*sk + jj*sk + kk] =
          static_cast<T>(getSyntheticSample(oi+ii, oj+jj, ok+kk, isfloat));
}

template <typename T>
MyStatistics
getSyntheticStats(std::int64_t si, std::int64_t sj, std::int64_t sk, bool isfloat)
{
  std::int64_t scnt = 0;
  double ssum = 0, sssq = 0, smin = 0, smax = 0;
  for (std::int64_t ii=0; ii<si; ++ii)
    for (std::int64_t jj=0; jj<sj; ++jj)
      if (jj != 200 || true) // skip the dead traces, but OpenZGY has no alpha.
        for (std::int64_t kk=0; kk<sk; ++kk)
          {
            T value = static_cast<T>(getSyntheticSample(ii, jj, kk, isfloat));
            double dvalue = static_cast<double>(value);
            if (std::isfinite(dvalue)) {
              if (scnt == 0)
                smin = smax = dvalue;
              else if (smin > dvalue)
                smin = dvalue;
              else if (smax < dvalue)
                smax = dvalue;
              scnt++;
              ssum += dvalue;
              sssq += dvalue*dvalue;
            }
          }
  return MyStatistics(scnt, ssum, sssq, smin, smax);
}

template<typename T> OpenZGY::SampleDataType getSampleVT()    { return OpenZGY::SampleDataType::unknown; }
template<> OpenZGY::SampleDataType getSampleVT<float>()       { return OpenZGY::SampleDataType::float32; }
template<> OpenZGY::SampleDataType getSampleVT<std::int16_t>(){ return OpenZGY::SampleDataType::int16; }
template<> OpenZGY::SampleDataType getSampleVT<std::int8_t>() { return OpenZGY::SampleDataType::int8; }

/**
 * Create a new ZGY file and return the open handle.
 */
template <typename AppT, typename StorageT>
std::shared_ptr<IZgyWriter>
createTestFile(const std::string& name, const size3i_t& size)
{
  return IZgyWriter::open
    (ZgyWriterArgs()
     .filename(name)
     .size(size[0], size[1], size[2])
     .datatype(getSampleVT<StorageT>())
     .datarange(0, 200)
     .ilstart(1234).ilinc(5)
     .xlstart(9012).xlinc(3)
     .zstart(789.012f).zinc(3.4f)
     .zunit(UnitDimension::length, "ft", 0.3048)
     .hunit(UnitDimension::length, "cm", 0.01));
}

template <typename AppT>
void checkTestFile(const std::string& name1, const std::string& name2, double eps)
{
  std::shared_ptr<IZgyReader> file1 = IZgyReader::open(name1);
  std::shared_ptr<IZgyReader> file2 = IZgyReader::open(name2);

  const size3i_t size = file1->size();
  const size3i_t size2 = file2->size();
  for (int ii=0; ii<3; ++ii)
    CPPUNIT_ASSERT_EQUAL(size[ii], size2[ii]);

  const std::int64_t nsamples = size[0] * size[1] * size[2];
  std::unique_ptr<AppT[]>brick1(new AppT[nsamples]);
  file1->read(size3i_t{0,0,0}, size, brick1.get(), /*lod=*/0);

  std::unique_ptr<AppT[]>brick2(new AppT[nsamples]);
  file2->read(size3i_t{0,0,0}, size, brick2.get(), /*lod=*/0);

  for (std::int64_t ii=0; ii<nsamples; ++ii) {
    // Explicit check first, because the macro is rather slow.
    if (std::fabs((double)brick1[ii] - (double)brick2[ii]) > eps) {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(brick1[ii], brick2[ii], eps);
      break;
    }
  }
}

/**
 * Create a file, Write synthetic test data to it.
 * The brick size to be used when writing is passed in.
 * It need not be aligned with the file's brick size.
 * If bricksize is passed as the file size, the entire
 * file is written in a single call.
 * Alpha tiles are always written this way, mostly because I am lazy.
 * If reopen=true, the file is closed and reopened
 * after every call - including the last one.
 */
template <typename AppT, typename StorageT>
void writeTestFile(
     const std::string& name,
     const size3i_t& size,
     const size3i_t& bricksize,
     bool reopen,
     MyStatistics *stat = 0,
     MyHistogram *hist = 0)
{
  //printf("writeTestFile bricked (%lld,%lld,%lld) reopen=%s\n",
  //     (long long)bricksize[0], (long long)bricksize[1], (long long)bricksize[2], reopen ? "true":"false");

  std::shared_ptr<IZgyWriter> w = createTestFile<AppT, StorageT>(name, size);

  bool isfloat =
    !std::numeric_limits<AppT>::is_integer &&
    !std::numeric_limits<StorageT>::is_integer;

  // No ALPHA support in OpenZGY. The old test flagged traces with j=0 as dead.
  // This might affect expected results.
  // The tests use size {117, 241, 333}, so 117*333 = 38961 samples should have been dead.

  std::unique_ptr<AppT[]>brick(new AppT[bricksize[0]*bricksize[1]*bricksize[2]]);
  size3i_t brickorig {0, 0, 0};
  for (brickorig[0]=0; brickorig[0]<size[0]; brickorig[0]+=bricksize[0]) {
    for (brickorig[1]=0; brickorig[1]<size[1]; brickorig[1]+=bricksize[1]) {
      for (brickorig[2]=0; brickorig[2]<size[2]; brickorig[2]+=bricksize[2]) {
        const size3i_t thissize
          {
            std::min((std::int64_t)bricksize[0], (std::int64_t)size[0]-brickorig[0]),
            std::min((std::int64_t)bricksize[1], (std::int64_t)size[1]-brickorig[1]),
            std::min((std::int64_t)bricksize[2], (std::int64_t)size[2]-brickorig[2])
          };
        getSyntheticBrick(brick.get(),
                          brickorig[0], brickorig[1], brickorig[2],
                          thissize[0], thissize[1], thissize[2],
                          isfloat);
        w->write(brickorig, thissize, brick.get());
        if (reopen) {
          w->finalize();
          w->close();
          w.reset();
          w = IZgyWriter::reopen(ZgyWriterArgs().filename(name));
        }
      }
    }
  }

  if (!reopen) {
    w->finalize();
  }

  // If caller wants them, read the statistics before the file is closed.
  if (stat) {
    *stat = MyStatistics(w->statistics());
  }
  if (hist) {
    *hist = MyHistogram(w->histogram());
  }

  // Close completely. Would happen anyway when w goes out of scope,
  // but this is simpler to debug.
  w->close();
  w.reset();
}

/**
 * As InternalZGY::HistogramBuilder::gethiststats() but with API types.
 * Calculate statistics directly from the histogram.
 * Pass is_int8 if the histogram was created from int8 data.
 */
SampleStatistics
getHistStats(const SampleHistogram& hh, const bool is_int8 = false)
{
  std::int64_t cnt{0};
  double sum{0};
  double ssq{0};
  double min{std::numeric_limits<double>::infinity()};
  double max{-std::numeric_limits<double>::infinity()};

  const std::int64_t nbins = (std::int64_t)hh.bins.size();
  const double begin = hh.minvalue;
  const double width = (hh.maxvalue - hh.minvalue) / (nbins - 1);
  const double slop = is_int8 ? 0 : width/2;

  for (std::int64_t binno = 0; binno < nbins; ++binno) {
    std::int64_t count = hh.bins[binno];
    double value = begin + binno*width;
    if (count != 0) {
      if (cnt == 0)
        min = value - slop; // first non-empty bin found
      max = value + slop; // last (so far) non-empty bin found.
    }
    cnt += count;
    sum += count * value;
    ssq += count * value * value;
  }
  return SampleStatistics(cnt, sum, ssq, min,max);
}

/**
 * Open a ZgY file, read its statistics and histogram, and close it again.
 * stats_from_histogram=true uses the file't histogram for both.
 * is_int8, only relevant when stats_from_histogram=true, statea that
 * the histogram was built from 8-bit data.
 */
void
readStatistics(const std::string& name, MyStatistics *stat, MyHistogram *hist, bool stats_from_histogram = false/*, bool is_int8 = false*/)
{
  std::shared_ptr<IZgyReader> r = IZgyReader::open(name);
  if (!stat)
    ;
  else if (stats_from_histogram)
    *stat = MyStatistics(getHistStats(r->histogram()));
  else
    *stat = MyStatistics(r->statistics());
  if (hist)
    *hist = MyHistogram(r->histogram());
}

/////////////////////////////////////////////////////////////////////////////
///   UNIT TESTS   //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Unit test to verify that statistics are correctly set when a file
 * is written incrementally, compared to when it is written in one chunk.
 *
 * The test will also test for regressions in the simple whole-file case.
 * In this case the comparison is done using the results from version 1
 * of the code, so it isn't a true unit test.
 *
 * Write the entire test file using the general api, close it.
 * Write a second test file in multiple chunks, close/reopen between each.
 * The chunks should be larger than a brick but not aligned with the brick
 * size, so we end up with writing a mix of partial and full bricks.
 *
 * Open both files and check the stored histogram and statistics.
 * both should match.
 *
 * Run this test for 4 type combinations. Either type can be float or short.
 * Can probably use aligned=false always in the automated tests, as this will
 * indirectly test some full bricks too.
 *
 * SINGLEPASS: The old premise is that both stats and histogram work
 * the same way when written incrementally. That is plain wrong in the
 * new code. Probably need to not reguister hist.multipass<f,f> at all.
 * hist.multipass<s,s> is already disabled due to changes from old ZGY.
 * hist.multipass<f,s> still works for some reason.
 */
template <typename AppT, typename StorageT>
void testMultiPass()
{
  LocalFileAutoDelete file1("test-multipass-1.zgy");
  LocalFileAutoDelete file2("test-multipass-2.zgy");
  LocalFileAutoDelete file3("test-multipass-3.zgy");
  LocalFileAutoDelete file4("test-multipass-4.zgy");

  const size3i_t size {117, 241, 333}; // 2 x 4 x 6 = 48 bricks
  const size3i_t bsalign {64, 128, 128}; // brick size aligned with file
  const size3i_t bsunalign {50, 100, 100}; // unaligned, partial writes

  bool isfloat =
    !std::numeric_limits<AppT>::is_integer &&
    !std::numeric_limits<StorageT>::is_integer;

  MyStatistics stat0, stat1, stat2, stat3, stat4;
  MyHistogram  hist0, hist1, hist2, hist3, hist4;

  //printf("\nWriting...\n"); fflush(stdout);

  // FILE WRITTEN ALL AT ONCE, FLUSHED BUT NOT CLOSED BEFORE READ BACK
  // FILE WRITTEN ALL AT ONCE, CLOSED BEFORE READ BACK STATISTICS
  writeTestFile<AppT, StorageT>(file1.name(), size, size, false, &stat0, &hist0);
  readStatistics(file1.name(), &stat1, &hist1);

  // FILE WRITTEN IN MISALIGNED BRICKS
  writeTestFile<AppT, StorageT>(file2.name(), size, bsunalign, false);
  readStatistics(file2.name(), &stat2, &hist2);

  // FILE WRITTEN IN ALIGNED BRICKS, CLOSE/REOPEN BETWEEN EACH
  writeTestFile<AppT, StorageT>(file3.name(), size, bsalign, true);
  readStatistics(file3.name(), &stat3, &hist3);

  // FILE WRITTEN IN MISALIGNED BRICKS, CLOSE/REOPEN BETWEEN EACH
  writeTestFile<AppT, StorageT>(file4.name(), size, bsunalign, true);
  readStatistics(file4.name(), &stat4, &hist4);

  //printf("Checking...\n"); fflush(stdout);

#if 0 // Uncomment if you need to look at the files using zgydymp
  unlink("file1.zgy");
  unlink("file2.zgy");
  unlink("file3.zgy");
  unlink("file4.zgy");
  link(file1.name().c_str(), "file1.zgy");
  link(file2.name().c_str(), "file2.zgy");
  link(file3.name().c_str(), "file3.zgy");
  link(file4.name().c_str(), "file4.zgy");
#endif

  // Check that the statistics are updated in memory on flush,
  // i.e. it shouldn't make a difference whether I reopen the file
  // before reading them. Only checked for the first read.
  // the assert itself will always fail due to the tacked-on "!",
  // it is writtem that way to get full control over the epsilon
  // value, while still making sure any error message will contain
  // the entire statistics and not just the one parameters that failed.

  if (!CPPUNIT_ASSERT(stat1.asExpected(stat0, 0.01)))
    CPPUNIT_ASSERT_EQUAL(stat0.toString(), stat1.toString()+"!");
  if (!CPPUNIT_ASSERT(hist1.asExpected(hist0, 0.01)))
    CPPUNIT_ASSERT_EQUAL(hist0.toString(), hist1.toString()+"!");

  // Expected value of statistics can be calculated by hand.
  MyStatistics statExpect =
    getSyntheticStats<AppT>(size[0], size[1], size[2], isfloat);

#if 0
  std::cout << "\nStats before and after close, and calculated\n"
            << stat0.toString() << "\n"
            << stat1.toString() << "\n"
            << statExpect.toString() << "\n";
#endif

  // Expected value of the histogram is not feasible to calculate exactly;
  // depending on the valuetypes used there is much rounding done.
  // So assumption is that the simple file1 case has a correct histogram.
  // But I can do *some* concrete tests. These are done later in this function.
  MyHistogram histExpect = hist1;

  // Fuzzy checks on the histogram by comparing it with the statistics,
  // since there was no proper 'expected' value for this.
  // Start with calculating the approximate statistics from the histogram.
  MyStatistics hist1stat = hist1.stats();
  CPPUNIT_ASSERT_EQUAL(stat1.cnt(), hist1.samplecount);
  CPPUNIT_ASSERT_EQUAL(stat1.cnt(), hist1stat.cnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat1.min(), hist1stat.min(), 5.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat1.max(), hist1stat.max(), 5.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat1.avg(), hist1stat.avg(), 5.0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat1.std(), hist1stat.std(), 5.0);

#if 0
  printf("\n[1]: stat %s\n     hist %s\n     data %s\n",
         stat1.toString().c_str(),
         hist1.stats().toString().c_str(),
         hist1.toString().c_str());
  printf("\n[2]: stat %s\n     hist %s\n     data %s\n",
         stat2.toString().c_str(),
         hist2.stats().toString().c_str(),
         hist2.toString().c_str());
  printf("\n[3]: stat %s\n     hist %s\n     data %s\n",
         stat3.toString().c_str(),
         hist3.stats().toString().c_str(),
         hist3.toString().c_str());
  printf("\n[4]: stat %s\n     hist %s\n     data %s\n",
         stat4.toString().c_str(),
         hist4.stats().toString().c_str(),
         hist4.toString().c_str());
#endif

  // Verify both the statistics and the histogram data.
  // Allow more slop when comparing the min/max range of the histogram.
  // It is subject to rounding and also the whims of the auto-expand feature.
  // The counts in the histogram and in the statistics should always match.

  CPPUNIT_ASSERT_EQUAL(stat1.cnt(), hist1.samplecount);
  if (!CPPUNIT_ASSERT(stat1.asExpected(statExpect, 0.1)))
    CPPUNIT_ASSERT_EQUAL(statExpect.toString(), stat1.toString()+"!");
  if (!CPPUNIT_ASSERT(hist1.asExpected(histExpect, 2.0)))
    CPPUNIT_ASSERT_EQUAL(histExpect.toString(), hist1.toString()+"!");

  CPPUNIT_ASSERT_EQUAL(stat2.cnt(), hist2.samplecount);
  if (!CPPUNIT_ASSERT(stat2.asExpected(statExpect, 0.1)))
    CPPUNIT_ASSERT_EQUAL(statExpect.toString(), stat2.toString()+"!");
  if (!CPPUNIT_ASSERT(hist2.asExpected(histExpect, 2.0)))
    CPPUNIT_ASSERT_EQUAL(histExpect.toString(), hist2.toString()+"!");

  CPPUNIT_ASSERT_EQUAL(stat3.cnt(), hist3.samplecount);
  if (!CPPUNIT_ASSERT(stat3.asExpected(statExpect, 0.1)))
    CPPUNIT_ASSERT_EQUAL(statExpect.toString(), stat3.toString()+"!");
  if (!CPPUNIT_ASSERT(hist3.asExpected(histExpect, 2.0)))
    CPPUNIT_ASSERT_EQUAL(histExpect.toString(), hist3.toString()+"!");

  CPPUNIT_ASSERT_EQUAL(stat4.cnt(), hist4.samplecount);
  if (!CPPUNIT_ASSERT(stat4.asExpected(statExpect, 0.1)))
    CPPUNIT_ASSERT_EQUAL(statExpect.toString(), stat4.toString()+"!");
  if (!CPPUNIT_ASSERT(hist4.asExpected(histExpect, 2.0)))
    CPPUNIT_ASSERT_EQUAL(histExpect.toString(), hist4.toString()+"!");

  // The contents of each of the 4 files should also be equal
  // at this point. While not related to the statistics work, it is
  // fairly cheap to read back the 4 files in the default value type
  // and compare the results. Extra credit for that...
  // Might as well use the file's storage type for reading it.
  //printf("Checking file contents...\n"); fflush(stdout);

  checkTestFile<StorageT>(file1.name(), file2.name(), 0.01);
  checkTestFile<StorageT>(file1.name(), file3.name(), 0.01);
  checkTestFile<StorageT>(file1.name(), file4.name(), 0.01);

  //printf("DONE\n"); fflush(stdout);
}

/**
 * Unit test for statistics in a partly written file.
 *
 * What are unwritten traces read back as?
 *
 *   1) when the requested trace in inside a partly written brick?
 *   2) when the requested trace was outside any written brick,
 *      but the read was using the general api and requested a
 *      region covering at least one brick?
 *   3) when the entire request was for unwritten data?
 *
 * The answer for 1,2 is "zero", but is it the zero of the storage VT?
 * No, currently real zeroes are returned. Unless the storage is integral
 * and codingmin/codingmax does not include zero. That triggers the
 * clipping issue tested elsewhere.
 * The answer for 3 is that you get a "brick not found" error.
 *
 * If the storage format is float so the histogram range is dynamic,
 * will these unwritten traces cause the histogram range to get
 * extended to cover zero - even if the real data doesn't?
 *
 * Answer is yes if the application writes partial bricks, no if
 * it doesn't. (The latter is not verified here). This is only an issue
 * for floating point values where 0 is not considered a valid data
 * point. I guess this is not a common scenario. But maybe we should
 * formalize the behavior to say that 0 will *always* be included
 * in the histogram range, no matter what. This gives more predictable
 * behavior.
 *
 * Are the unwritten traces included in the histogram?
 *
 * a) No, only the traces actually written.
 * b) Yes, but only those in partly written bricks.
 * c) Yes, all traces are counted whether written or not.
 * d) Yes, all traces including the brick padding at the end of the cube.
 *
 * The old answer was (a), current should be (b). It is not clear what
 * the correct behavior ought to be. And should I send out a heads-up
 * since I changed the behavior? Unfortunately (a) is just not doable
 * any longer since there is no way to later know exactly which of the
 * samples the application wrote. Unless the rest were filled with NaN,
 * and that would be a major change in strategy. And how to handle NaN
 * in integral data? Doable, but again this is a major design change.
 * It is easy to change the behavior to (c). Maybe that is more intuitive?
 *
 * Unrelated, but can maybe test this here:
 * If I rewrite the bricks, will the statistics get narrowed again?
 * Answer should be Yes, but the result might not be precise.
 *
 * TODO Documentation note: If the application ever write a partial
 * brick or an unaligned brick to the file using the general API,
 * the rest of the traces inside this brick are padded with zeroes.
 * Even if the application later overwrites the missing traces, some
 * side effects may remain.
 *
 * 1) If the storage type is float, the histogram has a dynamic range.
 *    It will get extended to include zero. It will not contract later.
 *    This could be a problem if the actual data range is small and
 *    far away from zero, e.g. +10000 to +10100.
 *
 * 2) The statistics will temporary have included the padding.
 *    When the tile is overwritten, statistiucs min/max might need
 *    to be shrinked. This cannot be done precisely but will be
 *    approximated using the histogram. Note that this is a very minor
 *    issue, mostly interesting to unit tests which might fail if they
 *    don't take this into account.
 *
 * 3) If the storage type is integral and the codingmin/codingmax range
 *    does not include 0 then the padding will instead end up as either
 *    the min or the max value, depending on which is closer.
 *
 * TODO, if application never writes partial bricks this ought not
 * to happen. Ideally we should have unit tests for this. But I'd prefer
 * to just state that we don't support it, and even force 0 to be part
 * of the dynamic histogram range.
 */
template <typename AppT, typename StorageT>
void testPartial()
{
  LocalFileAutoDelete file1("test-partial.zgy");
  const size3i_t size {64, 128, 150};

  std::shared_ptr<IZgyWriter> w =
    createTestFile<AppT, StorageT>(file1.name(), size);

  std::unique_ptr<AppT[]>brick(new AppT[size[0]*size[1]*size[2]]);
  for (int ii=0; ii<100; ++ii)
    brick[ii] = 42;
  brick[0] = 40;

  // Write a single trace (i=0, j=30) with length 100 (of the cube's 150).

  w->write(size3i_t{0,30,0}, size3i_t{1,1,100}, brick.get());
  w->finalize();
  w->close();
  w = IZgyWriter::reopen(ZgyWriterArgs().filename(file1.name()));

  // read back data outside the bricks we wrote
  // N/A for OpenZGY. Will not report BrickNotFound.
  //CPPUNIT_ASSERT(!a->ReadData(0, 0, 64, 0, 1, 64, 64, brick.get(), &status));
  //CPPUNIT_ASSERT_EQUAL(FileAccess::Error::BrickNotFound, status.zgyError());
  // reading back slice 1
  w->read(size3i_t{0,0,0}, size3i_t{1,size[1],size[2]}, brick.get());

  float livedata    = brick[30*size[2]+50];
  float deadinbrick = brick[30*size[2]+110];
  float deadoutside = brick[30*size[2]+140];

  MyStatistics stat0, stat1, stat2, stat3, stat4;
  MyHistogram  hist0, hist1, hist2, hist3, hist4;

  readStatistics(file1.name(), &stat1, &hist1);

  // Fill both the bricks we touched.
  for (int ii=0; ii<64*64*64; ++ii)
    brick[ii] = 30;

  w->write(size3i_t{0,0,0},  size3i_t{64,64,64}, brick.get());
  w->write(size3i_t{0,0,64}, size3i_t{64,64,64}, brick.get());
  w->finalize();
  w->close();
  w.reset();

  // Check the statistics again, see if any trace of the padding remains.

  readStatistics(file1.name(), &stat2, &hist2);

  if (verbose()) {
    std::cout << "\nAfter writing one trace\n"
              << stat1.toString() << "\n"
              << hist1.toString(2) << std::endl;
    std::cout << "After overwriting two bricks\n"
              << stat2.toString() << "\n"
              << hist2.toString(2) << std::endl;
  }

  // the expected absent value is 0 in storage, but what will a 0 storage
  // value look like after being converted by codingmin/codingmax?
  // NO, this is irrelevant. The point in the accessor where a missing
  // brick is converted to zeroes, is still in the "real" context.
  // So the only issue here is when 0 is not representable in the
  // underlying storage, i.e. if codingmin/codingmax does not include zero.
  // That is going to mess things up since statistics don't handle
  // clipping. There is a separate unit test for that case.
#if 0
  double absent = 0;
  if (std::numeric_limits<StorageT>::is_integer) {
    double oldmin = std::numeric_limits<StorageT>::min();
    double oldmax = std::numeric_limits<StorageT>::max();
    double newmin = fi.codingmin;
    double newmax = fi.codingmax;

    const double b = (newmax - newmin) / (oldmax - oldmin);
    const double a = newmin - oldmin*b;
    absent = a + b*absent;
  }
#endif

#if 0
  printf("RESULT: %g %g %g\n", livedata, deadinbrick, deadoutside);
  printf("%s\n", stat1.toString().c_str());
  printf("%s\n", hist1.stats().toString().c_str());
  printf("%s\n", stat2.toString().c_str());
  printf("%s\n", hist2.stats().toString().c_str());
#endif

  // written trace
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0, livedata, 0.1);

  // unwritten trace, in a brick we wrote
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, deadinbrick, 0.1);

  // unwritten trace, in a brick that did not exist
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, deadoutside, 0.1);

  // the unwritten traces will have extended the histogram range
  // Since the range is doubled on every expand, we cannot know
  // exactly what it ended up as - a big negative number is ok.
  // OpenZGY might have a smaller expand or none. Test is still good.
  if (hist1.minvalue > 0)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0, hist1.minvalue, 0.1);

  // the unwritten traces will have been included in the histogram,
  // but only for the bricks we touched part of. Ditto for stats.
  // CHANGED in OpenZGY: All samples in the brick are included.
  //Old: CPPUNIT_ASSERT_EQUAL(2*64*64*64, stat1.cnt());
  //Old: CPPUNIT_ASSERT_EQUAL(2*64*64*64LL, hist1.samplecount);
  CPPUNIT_ASSERT_EQUAL(std::int64_t(64*128*150), stat1.cnt());
  CPPUNIT_ASSERT_EQUAL(std::int64_t(64*128*150), hist1.samplecount);

  // After rewriting the bricks there should be no zeroes,
  // so the statistical min/max should have contracted again.
  // But the new range will not be accurate since it is based
  // on the histogram values.
  // CHANGED in OpenZGY: All samples in the brick are included.
  // CHANGED in OpenZGY: In two-pass mode, results will be exact.
  //Old: CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, stat2.min(), 2.0);
  //Old: CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, stat2.max(), 2.0);
  //Old: CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0, stat2.avg(), 0.01);
  const double filled = 2*64*64*64 / double(64*128*150);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, stat2.min(), 0.01);
  // SINGLEPASS: There is no full rebuild, so "42" still seen in stats.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42.0 /*was 30.0*/, stat2.max(), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(30.0*filled, stat2.avg(), 0.01);
}

/**
 * Unit test to check statistics and histogram after updating parts of a file.
 * Templates on application and storage value type, and on whether to perform
 * aligned (brick at a time) updates.
 *
 * Write the entire test file using the general api, close it.
 * The synthetic data used should include NaNs (if both valuetypes are float)
 * and some dead traces. Will use same test data generator as TestMultiPass.
 * Set the coding range on create to that of the synthetic data
 * plus the shift we are going to add later.
 *
 * Update a couple of bricks to contain a higher sample value.
 * Bricks are always 64x64x64, but for the unaligned case these
 * will actually overlap multiple bricks in the file.
 *
 * Considerations when choosing what to update:
 *   One brick set entirely to a constant high value.
 *   Other bricks get a constant shift up, large enough to ensure
 *   at least some of the contents end up higher than the original
 *   max value, but still within the file's coding range.
 *   At least one of the incremented bricks should cover data
 *   containing NaNs, and at least one should cover dead traces.
 *
 *     -> Set the following to constant 195:
 *        brick (0,1,1), unaligned: (0+delta, 1+delta, 1-delta)
 *        In the unaligned case it will then overlap both brick
 *        (0,1,0) originally containing all-NaN, and (0,2,0)
 *        originally containing all constant values.
 *
 *     -> Increment by constant 95:
 *        brick (1,1,1), unaligned: (1+delta, 1+delta, 1+delta)
 *        This covers at least some NaN cells.
 *        brick (1,1,3), unaligned: (1-delta, 1+delta, 3+delta)
 *        This covers the j=200 dead trace.
 *
 * After updating, close and reopen the file and check statistics.
 * The count should be the same, the sum should be "a bit" higher.
 *
 * Update again, this time restoring the bricks to their old contents.
 *
 * After updating, close and reopen the file and check statistics.
 * Everything should be back to the original except the max value
 * which can be up to 1/2 bin too wide. And if the storage is float,
 * the histogram range will have changed.
 *
 * Run this test for 4 type combinations. Either type can be float or short.
 */
template <typename AppT, typename StorageT, bool aligned>
void testUpdate()
{
  LocalFileAutoDelete file("test-update.zgy");
  const size3i_t size {117, 241, 333}; // 2 x 4 x 6 = 48 bricks
  const int delta = aligned ? 0 : 3;
  const bool isfloat =
    !std::numeric_limits<AppT>::is_integer &&
    !std::numeric_limits<StorageT>::is_integer;

  // The start point of the three 64x64x64 bricks we want to change.
  const size3i_t pos1 { 0*64+delta, 1*64+delta, 1*64-delta };
  const size3i_t pos2 { 1*64+delta, 1*64+delta, 1*64+delta };
  const size3i_t pos3 { 1*64-delta, 1*64+delta, 3*64+delta };
  const size3i_t chunk{ 64, 64, 64 };

  MyStatistics stat1, stat2, stat3;
  MyHistogram  hist1, hist2, hist3;

  writeTestFile<AppT, StorageT>(file.name(), size, size, false);
  readStatistics(file.name(), &stat1, &hist1);

  const std::int64_t nsamples = 64*64*64;
  std::unique_ptr<AppT[]>brick(new AppT[nsamples]);

  auto w = IZgyWriter::reopen(ZgyWriterArgs().filename(file.name()));

  std::int64_t wasnan = 0;

  // First brick, set to a constant value
  // Keep track of how many NaN values we clobber.
  // In the other two bricks we don't clobber any,
  // since we only increment data so NaNs remain.

  saferead(w, pos1, chunk, brick.get());
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    if (!std::isfinite(brick[ii]))
      ++wasnan;
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] = 195;
  safewrite(w, pos1, chunk, brick.get());

  // Second brick, increment value
  saferead(w, pos2, chunk, brick.get());
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] += 95;
  safewrite(w, pos2, chunk, brick.get());

  // Third brick, increment value
  saferead(w, pos3, chunk, brick.get());
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] += 95;
  safewrite(w, pos3, chunk, brick.get());

  // Flush and completely close the file.
  w->finalize();
  w->close();
  w.reset();

  readStatistics(file.name(), &stat2, &hist2);

  // Now, put things back the way we found them.
  w = IZgyWriter::reopen(ZgyWriterArgs().filename(file.name()));
  getSyntheticBrick(brick.get(), pos1[0], pos1[1], pos1[2], 64, 64, 64, isfloat);
  safewrite(w, pos1, chunk, brick.get());

  getSyntheticBrick(brick.get(), pos2[0], pos2[1], pos2[2], 64, 64, 64, isfloat);
  safewrite(w, pos2, chunk, brick.get());

  getSyntheticBrick(brick.get(), pos3[0], pos3[1], pos3[2], 64, 64, 64, isfloat);
  safewrite(w, pos3, chunk, brick.get());

  // Flush and completely close the file.
  w->finalize();
  w->close();
  w.reset();

  readStatistics(file.name(), &stat3, &hist3);

  if (verbose()) {
    std::cout << "\nAfter writing standard test file\n"
              << stat1.toString() << "\n"
              << hist1.toString(0) << std::endl;
    std::cout << "After writing larger values\n"
              << stat2.toString() << "\n"
              << hist2.toString(0) << std::endl;
    std::cout << "After putting things back\n"
              << stat3.toString() << "\n"
              << hist3.toString(0) << std::endl;
  }

#if 0
  printf("Clobbered %lld NaN\n", (long long)wasnan);
  printf("\n  beg: %s\n  upd: %s\n  end: %s\n",
         stat1.toString().c_str(),
         stat2.toString().c_str(),
         stat3.toString().c_str());
  printf("\n  beg: %s\n  upd: %s\n  end: %s\n",
         hist1.toString().c_str(),
         hist2.toString().c_str(),
         hist3.toString().c_str());
#endif

  // Ready to verify the changes in statistics.
  // Since we changed some NaNs to real values (in the float/float case),
  // the count for the second statistics should be greater.
  // The wey the bricks were chosen, NaNs were only overwritten
  // in the float,float,unaligned case.
  CPPUNIT_ASSERT_EQUAL(stat1.cnt()+wasnan, stat2.cnt());

  // Similarly, the second statistics ought to be higher because all the
  // real samples I changed ended up with a larger value.
  // Exactly how much depends on whether NaNs were involved, and also
  // rounding if the storage type was integral. The latter should
  // be taken care of by a larger slop value.
  // Again, I admit I didn't calculate the correct values by hand.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat1.avg()+(isfloat?9.0:9.0), stat2.avg(), 0.5);

  // The range of the second statistics should reflect the test data plus the
  // highest sample (195) that I changed. Some leeway dus to rounding
  // in case the storage value type is integral.
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, stat2.min(), 0.5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(195.0, stat2.max(), 0.5);

  // The final statistics after we put things back should be similar
  // to what we started with. But the max value might be a bit too high
  // (1/2 bin, i.e. around 0.5) because it had to be approximated
  // from the histogram.
  // TODO paal, why do I need a larger slop value on windows?
  // For some reason the rounding appears to work differently there.
  // stat1.max() is 100.002 on win, 99.9985 on linux. This minor difference
  // is amplified in the stat3.max, probably a coincidence since the
  // value happened to cross the magical 100.
  // SINGLEPASS: The min/max range is less reliable. No full rebuild.
  if (stat1.cnt() != stat3.cnt() ||
      std::fabs(stat1.avg() - stat3.avg()) > 1e-3 ||
      std::fabs(stat1.std() - stat3.std()) > 1e-3 /*||
      std::fabs(stat1.min() - stat3.min()) > 1e-3 ||  // SINGLEPASS:
      std::fabs(stat1.max() - stat3.max()) > 1.0 */)  // SINGLEPASS:
    CPPUNIT_ASSERT_EQUAL(stat1.toString(), stat3.toString()+"!");

  // Now check the statistics. We set the coding range large enough,
  // so even in the float case the histogram range should not change.
  // CHANGED in OpenZGY: data range is ignored for float.
  // CHANGED in OpenZGY: if float, the histigram is no longer fixed.
  const bool fixed_histogram = std::numeric_limits<StorageT>::is_integer;
  if (fixed_histogram) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.0, hist1.minvalue, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(200.0, hist1.maxvalue, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.0, hist2.minvalue, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(200.0, hist2.maxvalue, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(  0.0, hist3.minvalue, 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(200.0, hist3.maxvalue, 1e-3);

    // With 256 bins and range 0..200, the bin number
    // for value 195 should be around 250.
    // Rounding due to value types might affect this.
    CPPUNIT_ASSERT(hist1.bins[249] == 0);
    CPPUNIT_ASSERT(hist2.bins[249] > 0);
    CPPUNIT_ASSERT(hist3.bins[249] == 0);

    // Finally, the end histogram should match what we started with.
    if (!CPPUNIT_ASSERT(hist3.asExpected(hist1, 0.5)))
      CPPUNIT_ASSERT_EQUAL(hist1.toString(), hist3.toString());
  }
  else {
    // SINGLEPASS: The exact range of the histogram is an implementation
    // detail; it is allowed to be wider than necessary.
    // And it might also have shrunk back to the original range.
    TEST_CHECK(hist1.minvalue <= 0.0);
    TEST_CHECK(hist1.maxvalue >= 100.0);
    TEST_CHECK(hist2.minvalue <= 0.0);
    TEST_CHECK(hist2.maxvalue >= 100.0);
    TEST_CHECK(hist3.minvalue <= 0.0);
    TEST_CHECK(hist3.maxvalue >= 100.0);

    // with range 0..100, bin for 195 does not exist.
    // with range 0..195, bin for 195 is 255
    // SINGLEPASS: Implementation dependent.
    //CPPUNIT_ASSERT(hist2.bins[255] > 0);
    //CPPUNIT_ASSERT(hist3.bins[255] == 0);

    // Finally, the end histogram will NOT match what we started with.
    // It will have been epanded, unless a full rebuild was forced.
  }
}

/**
 * Unit test to verify correct handling of statistics when the application
 * tries to write data outside the codingrange set when the file was created.
 * This test is run with the application't value type always float.
 * If the storage type is float then the codingrange shouldn't have any
 * effect. If it is integral then the application's data will get clipped
 * and the question is whether the statistics will reflect the clipped
 * or the original values. To handle updates correctly, they ought to
 * reflect the clipped data.
 *
 * This test could have been combined with the main test above by having
 * a few samples updated to exceed the codingrange. But with the current
 * implementation this is known to fail. And when the application value
 * type is integral, you need to jump thru hoops (or at least tweak the
 * access scale range) to get these too-large values passed down.
 * A separate test is easier.
 *
 * There is a very similar issue where a sample can end up in the neighbor
 * bin due to numerical inaccuracy. Which means a later subtraction won't
 * find it. This is tricky to test for, but a fix for the clipping case
 * is likely to cover this one also.
 *
 * We only need to update a single brick for this test.
 */
template <typename StorageT>
void testClipped()
{
  LocalFileAutoDelete file("test-clipped.zgy");
  bool isfloat = !std::numeric_limits<StorageT>::is_integer;

  const size3i_t size { 64, 64, 64 };
  const size3i_t zero { 0, 0, 0 };
  auto w = createTestFile<float, StorageT>(file.name(), size);

  const std::int64_t nsamples = size[0]*size[1]*size[2];
  std::unique_ptr<float[]>brick(new float[nsamples]);

  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] = (ii%2==0) ? 50.0f : 333.0f;

  w->write(zero, size, brick.get());
  w->close();
  w.reset();

  MyStatistics stat;
  MyStatistics histstat;
  MyHistogram  hist;
  readStatistics(file.name(), &stat, &hist);
  readStatistics(file.name(), &histstat, nullptr, true);

  if (verbose()) {
    std::cout << "\n";
    std::cout << "  stat: " << stat.toString() << "\n";
    std::cout << "  stat: " << histstat.toString() << "\n";
    std::cout << "  hist: " << hist.toString() << std::endl;
  }

  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);

  if (isfloat) {
    // The codingrange of 0..200 is only used as the initial range
    // of the histogram, so both the statistics and the histogram
    // should reflect the real data.
    // If expect_dynamic_histo() the codingrange is ignored completely.
    // SINGLEPASS: the defaultvalue (0) apparently gets included in statistics.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0/*was 50*/, stat.min(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(333, stat.max(), 1e-3);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((50+333)/2.0, stat.avg(), 1e-3);

    // Do the same test extracting the statistics from the histogram.
    // This gives a less accurate version of the statistics.
    // In this case I am also verifying that the float histogram
    // is wide enough to hold both 50 and 333.
    // If expect_dynamic_histo() the codingrange is ignored,
    // and the histogram limits will be even more approximate.
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 50, histstat.min(), 1.95);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(333, histstat.max(), 0.95);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((50+333)/2.0, histstat.avg(), 0.95);
  } else {
    // The codingrange of 0..200 is the only range that can be stored.
    // So if the statistics include the illegal values, it is impossible
    // to subtract them again if this brick is changed.
    // Slop factors are a bit higher since the storage is integral.
    // SINGLEPASS: the defaultvalue (0) apparently gets included in statistics.
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0/*was 50*/, stat.min(), 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(200, stat.max(), 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((50+200)/2.0, stat.avg(), 0.5);

    // Do the same test extracting the statistics from the histogram.
    // This gives a less accurate version of the statistics.
    CPPUNIT_ASSERT_DOUBLES_EQUAL( 50, histstat.min(), 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(200, histstat.max(), 0.5);
    CPPUNIT_ASSERT_DOUBLES_EQUAL((50+200)/2.0, histstat.avg(), 0.5);
  }
}

/**
 * Verify that the min/max range reported by the Zgy statistics is in fact the
 * first and last bin center, not the full range of the histogram.
 * Also verify the choices made when the histogram is generated.
 *
 * TestEdgesSByte, TestEdgesShort, etc. could probably be made into a single
 * function templated on type. But there are enough tricky issues that
 * i prefer maintaining them separately for now.
 *
 * Background:
 *
 * Note the following caveat with histogram data. The histogram is described by
 * number of bins and a pair of min/max values. It is not obvious whether those
 * values describe the center value of the first and last bin, or whether they
 * are the open ended range representing the whole histogram. I.e. the first value
 * (inclusive) of the first bin and the last value (exclusive) of the last bin.
 * The Zgy library assumes the former definition.
 *
 * A closely related issue is how the min/max values are chosen, when the range and
 * type of input data is known. The Zgy library will set things up so the first bin
 * has its center value the data min, and the last bin has its center value the data max.
 *
 * For unsigned 8-bit data and a 256-bin histogram this is probably optimal.
 * Each bin will have a range of exactly 1: First bin [-0.5..+0.5>, Second [+0.5..+1.5>,
 * last [+254.5..+255.5>. There is no risk of rounding errors even if the values are
 * processed as float at some point. And if the input data is uniformly distributed,
 * the histogram will be uniform as well.
 *
 * You would get the same result by choosing bins [0..1>, [1..2>,...,[255..256> but in this
 * case there would be a real risk of rounding errors if the data was processad as float.
 *
 * For floating point data, the choice zgy makes is somewhat less ideal. It isn't incorrect,
 * but an evenly distributed input will typically result in the first and last bins having
 * less entries than the rest. Example: data range [0..100], 101 bins in histogram. Center
 * of first & last bin 0 and 100 respectively, which means first bin covers [-0.5..+0.5>
 * and last bin covers [+99.5..+100.5>. But this means that half the value range of the
 * first and last bin does not exist in the input.
 *
 * For 16-bit data the problem is worse. Logically if you have 256 bins and 16-bit data
 * you would expect each bin to have a value range of exactly 256. But since Zgy adds the
 * half-bin slop on each end, you end up with a bin range of ~257. And you can get
 * quantization errors in addition to the first and last bin being smaller than the rest.
 * The produced histogram is still technically correct but will not be as you expected it.
 *
 * Does this really matter? I suspect that if you have 16 or 32 bits data from something
 * that is really continuous data, then the histogram we produce is close enough.
 * For 8-bit seismic and even more importantly for discrete data then the details become
 * important. But as already explained, I believe the zgy library does the right thing here.
 */
void testEdgesSByte()
{
  LocalFileAutoDelete file1("test-edges-sbyte.zgy");
  const size3i_t size {128, 128, 128};
  const std::int64_t nsamples = size[0] * size[1] * size[2];

  const float codingmin = -128;
  const float codingmax = +127;
  auto w = IZgyWriter::open
    (ZgyWriterArgs()
     .filename(file1.name())
     .size(size[0], size[1], size[2])
     .datatype(getSampleVT<std::int8_t>())
     .datarange(codingmin, codingmax)
     .ilstart(1234).ilinc(5)
     .xlstart(9012).xlinc(3)
     .zstart(789.012f).zinc(3.4f)
     .zunit(UnitDimension::length, "ft", 0.3048)
     .hunit(UnitDimension::length, "cm", 0.01));

  // Write the whole file in one go
  std::unique_ptr<std::int8_t[]>brick(new std::int8_t[nsamples]);
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] = std::int8_t((ii%256) - 128);
  w->write(size3i_t{0,0,0}, size, brick.get());
  w->finalize();

  MyStatistics stat(w->statistics());
  MyHistogram hist(w->histogram());

  if (verbose()) {
    std::cout << "\n";
    std::cout << "  stat: " << stat.toString() << "\n";
    std::cout << "  hist: " << hist.toString() << std::endl;
  }

  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmin, stat.min(), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmax, stat.max(), 0.01);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL((codingmin+codingmax)/2.0, stat.avg(), 0.001);
  // Expect even distribution for this type
  for (std::size_t ii=0; ii<hist.bins.size(); ++ii) {
    CPPUNIT_ASSERT_EQUAL(std::int64_t(nsamples / 256), std::int64_t(hist.bins[ii]));
  }
}

void testEdgesShort()
{
  LocalFileAutoDelete file1("test-edges-short.zgy");
  const size3i_t size {128, 128, 128};
  const std::int64_t nsamples = size[0] * size[1] * size[2];

  const float codingmin = -128;
  const float codingmax = +127;
  auto w = IZgyWriter::open
    (ZgyWriterArgs()
     .filename(file1.name())
     .size(size[0], size[1], size[2])
     .datatype(getSampleVT<std::int16_t>())
     .datarange(codingmin, codingmax)
     .ilstart(1234).ilinc(5)
     .xlstart(9012).xlinc(3)
     .zstart(789.012f).zinc(3.4f)
     .zunit(UnitDimension::length, "ft", 0.3048)
     .hunit(UnitDimension::length, "cm", 0.01));

  // Write the whole file in one go
  std::unique_ptr<std::int16_t[]>brick(new std::int16_t[nsamples]);
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] = std::int16_t((ii%65536) - 32768);
  w->write(size3i_t{0,0,0}, size, brick.get());
  w->finalize();

  MyStatistics stat(w->statistics());
  MyHistogram hist(w->histogram());

  if (verbose())
    std::cout << "\n"
              << "  stat: " << stat.toString() << "\n"
              << "  hist: " << hist.toString() << std::endl;

  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmin, stat.min(), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmax, stat.max(), 0.01);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL((codingmin+codingmax)/2.0, stat.avg(), 0.001);
  // There are 256 bins. Roughly speaking, the first and last bin will get only half as many samples.
  // Quantization errors causes this to be not entirely exact.
  std::int64_t expected_bin_count = nsamples / 255;
  std::int64_t expected_bin_count_0 = (nsamples - (expected_bin_count*254)) / 2;
  for (std::size_t ii=0; ii<hist.bins.size(); ++ii) {
    CPPUNIT_ASSERT_EQUAL(ii==0 || ii + 1 == hist.bins.size() ? expected_bin_count_0 : expected_bin_count, hist.bins[ii]);
  }
}

void testEdgesFloat()
{
  LocalFileAutoDelete file1("test-edges-float.zgy");
  const size3i_t size {128, 128, 128};
  const std::int64_t nsamples = size[0] * size[1] * size[2];

  const float codingmin = 0;
  const float codingmax = 100;
  auto w = IZgyWriter::open
    (ZgyWriterArgs()
     .filename(file1.name())
     .size(size[0], size[1], size[2])
     .datatype(getSampleVT<std::int16_t>())
     .datarange(codingmin, codingmax)
     .ilstart(1234).ilinc(5)
     .xlstart(9012).xlinc(3)
     .zstart(789.012f).zinc(3.4f)
     .zunit(UnitDimension::length, "ft", 0.3048)
     .hunit(UnitDimension::length, "cm", 0.01));

  // Write the whole file in one go
  std::unique_ptr<float[]>brick(new float[nsamples]);
  for (std::int64_t ii=0; ii<nsamples; ++ii)
    brick[ii] = float((ii*((codingmax-codingmin) / nsamples)) + codingmin);
  w->write(size3i_t{0,0,0}, size, brick.get());
  w->finalize();

  MyStatistics stat(w->statistics());
  MyHistogram hist(w->histogram());

  if (verbose())
    std::cout << "\n"
              << "  stat: " << stat.toString() << "\n"
              << "  hist: " << hist.toString() << std::endl;

  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmin, stat.min(), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)codingmax, stat.max(), 0.01);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_DOUBLES_EQUAL((codingmin+codingmax)/2.0, stat.avg(), 0.001);
  // There are 256 bins. Roughly speaking, the first and last bin will get only half as many samples.
  // Quantization errors causes this to be not entirely exact.
  std::int64_t expected_bin_count = nsamples / 255;
  for (std::size_t ii=1; ii<hist.bins.size()-1; ++ii) {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bin_count, hist.bins[ii], 3.0);
  }
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bin_count/2, hist.bins[0], expected_bin_count/100);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected_bin_count/2, hist.bins[hist.bins.size()-1], expected_bin_count/100);
}

// This test new in OpenZGY.
static void testNaN()
{
  LocalFileAutoDelete lad("test-hist-nan.zgy");
  const size3i_t zero {0, 0, 0};
  const size3i_t size {117, 241, 333}; // 2 x 4 x 6 = 48 bricks

  MyStatistics stat0, stat1, stat2, stat3;
  MyHistogram  hist0, hist1, hist2, hist3;

  auto w = createTestFile<float,float>(lad.name(), size);
  const float nan{std::numeric_limits<float>::quiet_NaN()};
  if (verbose())
    std::cout << "\n% writeconst all nan" << std::endl;
  w->writeconst(zero, size, &nan);
  w->close();
  if (verbose())
    std::cout << "% closed after create all nan\n\n" << std::flush;
  readStatistics(lad.name(), &stat0, &hist0);

  w = IZgyWriter::reopen(ZgyWriterArgs().filename(lad.name()));
  const float sixteen{16};
  if (verbose())
    std::cout << "% write sample value 16" << std::endl;
  w->write(size3i_t{0,0,0}, size3i_t{1,1,1}, &sixteen);
  w->close();
  if (verbose())
    std::cout << "% closed after written sample value 16\n\n" << std::flush;
  readStatistics(lad.name(), &stat1, &hist1);

  w = IZgyWriter::reopen(ZgyWriterArgs().filename(lad.name()));
  const float fortytwo{42};
  if (verbose())
    std::cout << "% write sample value 42" << std::endl;
  w->write(size3i_t{0,0,1}, size3i_t{1,1,1}, &fortytwo);
  w->close();
  if (verbose())
    std::cout << "% closed after written sample value 42\n\n" << std::flush;
  readStatistics(lad.name(), &stat2, &hist2);

  w = IZgyWriter::reopen(ZgyWriterArgs().filename(lad.name()));
  const float bignumber{1000};
  if (verbose())
    std::cout << "% write sample value 1000 in another brick" << std::endl;
  w->write(size3i_t{64,0,0}, size3i_t{1,1,1}, &bignumber);
  w->close();
  if (verbose())
    std::cout << "% closed after written sample value 1000 in another brick\n\n" << std::flush;
  readStatistics(lad.name(), &stat3, &hist3);

  if (verbose()) {
    std::cout << "\nAfter writing all NaN\n"
              << stat0.toString() << "\n"
              << hist0.toString(0) << std::endl;
    std::cout << "After writing 16\n"
              << stat1.toString() << "\n"
              << hist1.toString(2) << std::endl;
    std::cout << "After writing 42\n"
              << stat2.toString() << "\n"
              << hist2.toString(2) << std::endl;
    std::cout << "After writing 1000\n"
              << stat3.toString() << "\n"
              << hist3.toString(2) << std::endl;
  }

  // The all-NaN file gets a bogus histogram range -128..+127 with
  // zero samples. TODO-WIP-BrickedAPI: That particular histogram
  // should be recognized, and completely ignored once a real sample
  // is seen. If statistical count is 0, ignore existing bogus min/max
  // in histogram. (and assert samplecount is 0). Else if statistical
  // min == max (and assert all histo entries in one bin), histogram
  // range set to place all samples in the first bin. I.e. set
  // histogram min to stat min. Histogram max set to 1 if value==0, or
  // 0 if value<0, or 2*value. Also in this case, histo range is to be
  // reset completely when more data arrives.

  CPPUNIT_ASSERT_EQUAL(0, stat0.cnt());
  CPPUNIT_ASSERT_EQUAL(0, stat0.min());
  CPPUNIT_ASSERT_EQUAL(0, stat0.max());
  CPPUNIT_ASSERT_EQUAL(0, hist0.samplecount);
  // SINGLEPASS: With no samples, and the rule that a dynamic histogram
  // must be symmetric which means it must contain zero, the range
  // should be +/- "tiny"
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, hist0.minvalue, 1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, hist0.maxvalue, 1e-10);
  CPPUNIT_ASSERT(hist0.minvalue < 0);
  CPPUNIT_ASSERT(hist0.maxvalue > 0);

  CPPUNIT_ASSERT_EQUAL(1, stat1.cnt());
  CPPUNIT_ASSERT_EQUAL(sixteen, stat1.min());
  CPPUNIT_ASSERT_EQUAL(sixteen, stat1.max());
  CPPUNIT_ASSERT_EQUAL(1, hist1.samplecount);
  // SINGLEPASS: Histogram range very roughly 16..16,
  // assume no wider than 8..24. Might as well test for
  // actual values from a test run. If fail, and still
  // within the expected range, feel free to update "expect".
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+12, hist1.minvalue, 0.2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+20, hist1.maxvalue, 0.2);

  CPPUNIT_ASSERT_EQUAL(2, stat2.cnt());
  CPPUNIT_ASSERT_EQUAL(sixteen, stat2.min());
  CPPUNIT_ASSERT_EQUAL(fortytwo, stat2.max());
  // SINGLEPASS: After reopen, histo range is fixed as +12..+20
  // so the "42" value gets clipped out of the histogram but
  // will be included in the statistics and the bulk data.
  CPPUNIT_ASSERT_EQUAL(1/*was 2*/, hist2.samplecount);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+12, hist2.minvalue, 0.2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+20, hist2.maxvalue, 0.2);

  // SINGLEPASS: See above.
  CPPUNIT_ASSERT_EQUAL(3, stat3.cnt());
  CPPUNIT_ASSERT_EQUAL(sixteen, stat3.min());
  CPPUNIT_ASSERT_EQUAL(bignumber, stat3.max());
  CPPUNIT_ASSERT_EQUAL(1/*was 3*/, hist3.samplecount);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+12, hist3.minvalue, 0.2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(+20, hist3.maxvalue, 0.2);

  // Note that it is unclear whether the defaultvalue gets included
  // in the statistical and/or histogram range even when no samples
  // have that value. In this particular test it isn't. If the file
  // had not been nan-filled then it probably would have been.
}

// This test new in OpenZGY. Some overlap with testNaN.
// The focus here is on resetting the histogram when
// transitioning from zero or one values to multiple.
// Since NaN handling isn't up to spec yet, disable
// some of those checks.
//
//   EMPTY: No data
//   SINGLE: Single value
//   NONZERO: Multiple values, all positive or all negative
//   ONESIDED: Multiple values, at least one zero, the rest share the same sign.
//   MIXED: At least one positive and one negative.
template<bool reuse>
static void testSpecial()
{
  // TODO, both close/reopen and just finalize.
  LocalFileAutoDelete lad("test-hist-special.zgy");
  const size3i_t zero {0, 0, 0};
  const size3i_t size {117, 241, 333}; // 2 x 4 x 6 = 48 bricks
  const std::int64_t nsamples = size[0] * size[1] * size[2];

  std::shared_ptr<IZgyWriter> w;
  float value{0};
  MyStatistics stat;
  MyHistogram  hist;

  auto next = [&w, &stat, &hist, &lad](const std::string& message) {
                w->finalize();
                if (reuse) {
                  stat = MyStatistics(w->statistics());
                  hist = MyHistogram(w->histogram());
                }
                else {
                  w->close();
                  readStatistics(lad.name(), &stat, &hist);
                  w = IZgyWriter::reopen(ZgyWriterArgs().filename(lad.name()));
                }
                if (verbose())
                  std::cout << message << "\n"
                            << stat.toString() << "\n"
                            << hist.toString(2) << "\n";
              };

  if (verbose())
    std::cout << "\n";

  // SINGLE (unwritten)
  w = createTestFile<float,float>(lad.name(), size);
  next("SINGLE (unwritten)");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat.min(), stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);

  // SINGLE (written, nonzero)
  // Writing the entire survey should reset histogram and stats.
  // Getting rid of the bogus -128..+127.
  value = 42;
  w->writeconst(zero, size, &value);
  next("SINGLE (written, nonzero)");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat.min(), stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);
  CPPUNIT_ASSERT(hist.minvalue > -64.001);
  CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  // NONZERO
  value = 15;
  w->writeconst(zero, size3i_t{64,64,64}, &value);
  next("NONZERO");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);
  CPPUNIT_ASSERT(hist.minvalue < 15.001);
  CPPUNIT_ASSERT(hist.maxvalue > 41.999);
  CPPUNIT_ASSERT(hist.minvalue > -64.001);
  CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  // ONESIDED
  value = 0;
  w->writeconst(size3i_t{0,0,64}, size3i_t{64,64,64}, &value);
  next("ONESIDED");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  // SINGLEPASS: Too few entries in histogram, because the range
  // could not be extended in this case.
  //CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);
  //CPPUNIT_ASSERT(hist.minvalue < 0.001);
  CPPUNIT_ASSERT(hist.maxvalue > 41.999);
  CPPUNIT_ASSERT(hist.minvalue > -64.001);
  CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  // MIXED
  value = -2;
  w->writeconst(size3i_t{0,64,64}, size3i_t{64,64,64}, &value);
  next("MIXED");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(nsamples, stat.cnt());
  // SINGLEPASS: Too few entries in histogram, because the range
  // could not be extended in this case.
  //CPPUNIT_ASSERT_EQUAL(nsamples, hist.samplecount);
  //CPPUNIT_ASSERT(hist.minvalue < -1.999);
  CPPUNIT_ASSERT(hist.maxvalue > 41.999);
  CPPUNIT_ASSERT(hist.minvalue > -64.001);
  CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  // EMPTY
  // NaN isn't officially supported, but only way of getting EMPTY.
  // Writing the entire file should reset everything.
  value = std::numeric_limits<float>::quiet_NaN();
  w->writeconst(zero, size, &value);
  next("EMPTY");
  //FAILS! CPPUNIT_ASSERT(stat.min() > stat.max());
  CPPUNIT_ASSERT(hist.minvalue < hist.maxvalue); // Always sane
  CPPUNIT_ASSERT_EQUAL(0, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(0, hist.samplecount);

  // SINGLE
  value = 42;
  w->writeconst(size3i_t{0,0,0}, size3i_t{64,64,64}, &value);
  next("SINGLE");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(stat.min(), stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(64*64*64, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(64*64*64, hist.samplecount);
  //FAILS! CPPUNIT_ASSERT(hist.minvalue > -64.001);
  //FAILS! CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  // EMPTY->MIXED
  value = std::numeric_limits<float>::quiet_NaN();
  w->writeconst(zero, size, &value);
  value = 42;
  w->writeconst(size3i_t{0,0,0}, size3i_t{64,64,64}, &value);
  value = 15;
  w->writeconst(size3i_t{0,0,64}, size3i_t{64,64,64}, &value);
  value = -2;
  w->writeconst(size3i_t{0,64,0}, size3i_t{64,64,64}, &value);
  next("EMPTY->MIXED");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-2, stat.min(), 0.001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(42, stat.max(), 0.001);
  CPPUNIT_ASSERT_EQUAL(3*64*64*64, stat.cnt());
  CPPUNIT_ASSERT_EQUAL(3*64*64*64, hist.samplecount);
  CPPUNIT_ASSERT(hist.minvalue < -1.999);
  CPPUNIT_ASSERT(hist.maxvalue > 41.999);
  CPPUNIT_ASSERT(hist.minvalue > -64.001);
  CPPUNIT_ASSERT(hist.maxvalue < +64.001);

  w->close();
}

/**
 * Request a fixed-range histogram in spite of this being a float file.
 * Exclude 0, the default value, from the histogram.
 */
static void
testForceRangeInit()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("test-force-range-init.zgy");
  std::shared_ptr<IZgyWriter> writer = IZgyWriter::open
    (ZgyWriterArgsV2()
     .filename(lad.name())
     .historange(10, 100)
     .size(128, 128, 128)
     .datatype(OpenZGY::SampleDataType::float32));
  const float sixteen{16};
  const float fortytwo{42};
  const float onegross{144};
  writer->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &sixteen);
  writer->writeconst(index3_t{0,0,60}, index3_t{1,1,20}, &fortytwo);
  writer->writeconst(index3_t{0,0,127}, index3_t{1,1,1}, &onegross);
  writer->finalize();
  writer->close();

  MyStatistics stat;
  MyHistogram  hist;
  readStatistics(lad.name(), &stat, &hist);
  if (verbose())
    std::cout << "\nchosen  histo: " << hist.toString(2) << std::endl;
  TEST_EQUAL(stat.cnt(), 128*128*128);
  TEST_EQUAL(stat.sum(), sixteen*10 + fortytwo*20 + onegross*1);
  TEST_EQUAL(stat.min(), 0);
  TEST_EQUAL(stat.max(), onegross);
  TEST_EQUAL(hist.minvalue, 10);
  TEST_EQUAL(hist.maxvalue, 100);
  // Both the default value (0) and one sample (onegross) are clipped,
  // so the histogram only sees the 10+20 in-range samples.
  TEST_EQUAL(hist.samplecount, 30);
}

static void
testForceRangeUpdate()
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("test-force-range-init.zgy");
  const float zero{0};
  const float sixteen{16};
  const float fortytwo{42};
  const float onemill{1000000};

  {
    // Dynamic histogram range, somewhat wider than 16..42
    std::shared_ptr<IZgyWriter> writer = IZgyWriter::open
      (ZgyWriterArgsV2()
       .filename(lad.name())
       .size(128, 128, 128)
       .datatype(OpenZGY::SampleDataType::float32));
    writer->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &sixteen);
    writer->writeconst(index3_t{0,0,60}, index3_t{1,1,20}, &fortytwo);
    writer->close();

    MyStatistics stat0;
    MyHistogram  hist0;
    readStatistics(lad.name(), &stat0, &hist0);
    if (verbose())
      std::cout << "\ndynamic histo: " << hist0.toString(0) << std::endl;
    TEST_EQUAL(stat0.cnt(), 128*128*128);
    TEST_EQUAL(stat0.sum(), sixteen*10 + fortytwo*20);
    TEST_EQUAL(stat0.min(), 0);
    TEST_EQUAL(stat0.max(), fortytwo);
    TEST_CHECK(hist0.minvalue <= 16);
    TEST_CHECK(hist0.maxvalue >= 42);
    TEST_EQUAL(hist0.samplecount, stat0.cnt());
  }

  {
    // Update the file, try to force a wider range won't work when data exists.
    std::shared_ptr<IZgyWriter> writer = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .filename(lad.name())
       .historange(10, 1000000));
    writer->writeconst(index3_t{0,0,127}, index3_t{1,1,1}, &onemill);
    writer->close();
    MyStatistics stat1;
    MyHistogram  hist1;
    readStatistics(lad.name(), &stat1, &hist1);
    if (verbose())
      std::cout << "Force ignored: " << hist1.toString(0) << std::endl;
    TEST_EQUAL(stat1.cnt(), 128*128*128);
    TEST_EQUAL(stat1.sum(), sixteen*10 + fortytwo*20 + onemill*1);
    TEST_EQUAL(stat1.min(), 0);
    TEST_EQUAL(stat1.max(), onemill);
    TEST_CHECK(hist1.minvalue <= 10);
    TEST_CHECK(hist1.maxvalue >= 42);
    TEST_EQUAL(hist1.samplecount, stat1.cnt() - 1); // One value clipped
  }

  {
    // Update the file, clear and start over, should allow changing the range
    // in this case to the range remembered from the call to reopen().
    std::shared_ptr<IZgyWriter> writer = IZgyWriter::reopen
      (ZgyWriterArgsV2()
       .filename(lad.name())
       .historange(10, 2000000));
    writer->writeconst(index3_t{0,0,0}, writer->size(), &zero);
    writer->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &sixteen);
    writer->writeconst(index3_t{0,0,60}, index3_t{1,1,20}, &fortytwo);
    writer->writeconst(index3_t{0,0,127}, index3_t{1,1,1}, &onemill);
    writer->close();
    MyStatistics stat2;
    MyHistogram  hist2;
    readStatistics(lad.name(), &stat2, &hist2);
    if (verbose())
      std::cout << "Force applied: " << hist2.toString(0) << std::endl;
    TEST_EQUAL(stat2.cnt(), 128*128*128);
    TEST_EQUAL(stat2.sum(), sixteen*10 + fortytwo*20 + onemill*1);
    TEST_EQUAL(stat2.min(), 0);
    TEST_EQUAL(stat2.max(), onemill);
    TEST_EQUAL(hist2.minvalue, 10);
    TEST_EQUAL(hist2.maxvalue, 2000000);
    TEST_EQUAL(hist2.samplecount, stat2.cnt());
  }
}

/**
 * Test whether the current statistics and histogram can be read while
 * a file is being written, from the same file handle. Spoiler alert:
 * Yes, but only after a call to finalize().
 */
static void
testStatHistFromWriter(bool finalize)
{
  typedef std::array<std::int64_t,3> index3_t;
  LocalFileAutoDelete lad("test-stat-hist-from-writer.zgy");
  constexpr float sixteen{16};
  constexpr float fortytwo{42};
  constexpr float minusfive{-5};
  constexpr index3_t size{128, 128, 256};
  constexpr std::int64_t total_samples{size[0] * size[1] * size[2]};

  static auto dump = [](const std::string& msg,
                        const SampleStatistics& s,
                        const SampleHistogram& h)
      {
        if (verbose()) {
          std::cout << msg << ":\n"
                    << "  stat: " << MyStatistics(s).toString() << "\n"
                    << "  hist: " << MyHistogram(h).toString(2) << std::endl;
        }
      };

  // Dynamic histogram range, somewhat wider than 16..42
  std::shared_ptr<IZgyWriter> writer = IZgyWriter::open
    (ZgyWriterArgsV2()
     .filename(lad.name())
     .size(128, 128, 256)
     .datatype(OpenZGY::SampleDataType::float32));

  if (finalize)
    writer->finalize();

  {
    SampleStatistics stat0 = writer->statistics();
    SampleHistogram  hist0 = writer->histogram();
    dump("Empty file", stat0, hist0);
    TEST_EQUAL_FLOAT(stat0.sum, 0, 0.0001);
    TEST_EQUAL_FLOAT(stat0.min, 0, 0.0001);
    TEST_EQUAL_FLOAT(stat0.max, 0, 0.0001);
    TEST_EQUAL(hist0.bins[127], finalize ? total_samples : 0);
  }

  // Write 30 samples.
  writer->writeconst(index3_t{0,0,0}, index3_t{1,1,10}, &sixteen);
  writer->writeconst(index3_t{0,0,60}, index3_t{1,1,20}, &fortytwo);
  if (finalize)
    writer->finalize();

  {
    SampleStatistics stat1 = writer->statistics();
    SampleHistogram  hist1 = writer->histogram();
    dump("Written 30", stat1, hist1);
    if (finalize) {
      TEST_EQUAL_FLOAT(stat1.sum, sixteen*10 + fortytwo*20, 0.001);
      TEST_EQUAL_FLOAT(stat1.min, 0, 0.001); // defaultvalue included
      TEST_EQUAL_FLOAT(stat1.max, 42, 0.001);
      TEST_EQUAL(hist1.bins[15], total_samples - 10 - 20);
    }
    else {
      TEST_EQUAL_FLOAT(stat1.sum, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat1.min, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat1.max, 0, 0.0001);
      TEST_EQUAL(hist1.bins[127], 0);
    }
  }

  // Write 30 samples more. // histo range will not grow.
  writer->writeconst(index3_t{0,0,150}, index3_t{1,1,30}, &minusfive);
  if (finalize)
    writer->finalize();

  {
    SampleStatistics stat2 = writer->statistics();
    SampleHistogram  hist2 = writer->histogram();
    dump("Written 30 more", stat2, hist2);
    if (finalize) {
      TEST_EQUAL_FLOAT(stat2.sum, sixteen*10 + fortytwo*20 + minusfive*30, 0.001);
      TEST_EQUAL_FLOAT(stat2.min, -5, 0.001);
      TEST_EQUAL_FLOAT(stat2.max, 42, 0.001);
      // Note that if close/reopen was done, the histogram range would
      // not change and the "zero" bin would have remained 15. Using
      // just a finalize() it isn't too late. The "zero" bin will move.
      TEST_EQUAL(hist2.bins[29], total_samples - 10 - 20 - 30);
    }
    else {
      TEST_EQUAL_FLOAT(stat2.sum, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat2.min, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat2.max, 0, 0.0001);
      TEST_EQUAL(hist2.bins[127], 0);
    }
  }

  // Updating an existing brick.
  writer->writeconst(index3_t{0,0,150}, index3_t{1,1,30}, &sixteen);
  if (finalize)
    writer->finalize();

  {
    SampleStatistics stat3 = writer->statistics();
    SampleHistogram  hist3 = writer->histogram();
    dump("Updated brick", stat3, hist3);
    if (finalize) {
      TEST_EQUAL_FLOAT(stat3.sum, sixteen*10 + fortytwo*20 + sixteen*30, 0.001);
      TEST_EQUAL_FLOAT(stat3.min, -5, 0.001);
      TEST_EQUAL_FLOAT(stat3.max, 42, 0.001); // will not shrink
      TEST_EQUAL(hist3.bins[15], total_samples - 10 - 20 - 30);
    }
    else {
      TEST_EQUAL_FLOAT(stat3.sum, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat3.min, 0, 0.0001);
      TEST_EQUAL_FLOAT(stat3.max, 0, 0.0001);
      TEST_EQUAL(hist3.bins[127], 0);
    }
  }

  writer->finalize();

  {
    SampleStatistics stat4 = writer->statistics();
    SampleHistogram  hist4 = writer->histogram();
    dump("Explicit finalize", stat4, hist4);
    TEST_EQUAL_FLOAT(stat4.sum, sixteen*10 + fortytwo*20 + sixteen*30, 0.001);
    TEST_EQUAL_FLOAT(stat4.min, -5, 0.001);
    TEST_EQUAL_FLOAT(stat4.max, 42, 0.001);
    TEST_EQUAL(hist4.bins[15], total_samples - 10 - 20 - 30);
  }

  writer->close();
}

/**
 * Test that annotation and coordinates can be read from a file handle
 * open for write, without doing a finalize() or close()/reopen() first.
 *
 * The test is placed here because it is very similar to the test for
 * reading statistics and histogram. It doesn't have anything to do
 * with statistics and histogram.
 */
static void
testAnnotFromWriter()
{
  LocalFileAutoDelete lad("test-annot-from-writer.zgy");

  std::shared_ptr<IZgyWriter> writer = IZgyWriter::open
    (ZgyWriterArgsV2()
     .filename(lad.name())
     .size(128, 128, 256)
     .datatype(OpenZGY::SampleDataType::float32)
     .zunit(UnitDimension::time, "ms", 0.001)
     .hunit(UnitDimension::length, "ft", 0.3048)
     .ilstart(1)
     .ilinc(2)
     .xlstart(500)
     .xlinc(5)
     .zstart(100)
     .zinc(4)
     .corners(ZgyWriterArgs::corners_t{{{1,2},{3,4},{5,6},{7,8}}}));

  TEST_EQUAL((int)writer->zunitdim(), (int)UnitDimension::time);
  TEST_EQUAL((int)writer->hunitdim(), (int)UnitDimension::length);
  TEST_EQUAL(writer->zunitname(), std::string("ms"));
  TEST_EQUAL(writer->hunitname(), std::string("ft"));
  TEST_EQUAL_FLOAT(writer->zunitfactor(), 0.001, 0.00001);
  TEST_EQUAL_FLOAT(writer->hunitfactor(), 0.3048, 0.00001);


  TEST_EQUAL_FLOAT(writer->annotstart()[0], 1.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotstart()[1], 500.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotinc()[0], 2.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotinc()[1], 5.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->zstart(), 100.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->zinc(), 4.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->corners()[0][0], 1.0f, 0.001);

  writer->close();

  // Same test, but with a file that has been reopened
  // in order to update some metadata.

  writer = IZgyWriter::reopen
    (ZgyWriterArgsV2()
     .filename(lad.name())
     .ilstart(10)
     .ilinc(20)
     .xlstart(5000)
     .xlinc(50)
     .zstart(1000)
     .zinc(40)
     .corners(ZgyWriterArgs::corners_t{{{9,2},{3,4},{5,6},{7,8}}}));

  TEST_EQUAL_FLOAT(writer->annotstart()[0], 10.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotstart()[1], 5000.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotinc()[0], 20.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->annotinc()[1], 50.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->zstart(), 1000.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->zinc(), 40.0f, 0.001);
  TEST_EQUAL_FLOAT(writer->corners()[0][0], 9.0f, 0.001);

  writer->close();
}

} // namespace for tests

namespace {
  class Register
  {
  public:
    Register()
    {
      // AppT, StorageT
      //FAILS, see top of file.
      //register_test("hist.multipass<s,s>",     [](){testMultiPass<short,short>();});
      register_test("hist.multipass<f,s>",     [](){testMultiPass<float,short>();});
      // SINGLEPASS: Fails, see function header.
      //register_test("hist.multipass<f,f>",     [](){testMultiPass<float,float>();});
      // AppT, StorageT
      //register_test("hist.partial<s,s>",       [](){testPartial<short,short>();});
      register_test("hist.partial<f,s>",       [](){testPartial<float,short>();});
      register_test("hist.partial<f,f>",       [](){testPartial<float,float>();});
      // AppT, StorageT, Aligned
      //register_test("hist.update<s,s,T>",      [](){testUpdate<short,short,true>();});
      register_test("hist.update<f,s,T>",      [](){testUpdate<float,short,true>();});
      register_test("hist.update<f,f,T>",      [](){testUpdate<float,float,true>();});
      //register_test("hist.update<s,s,F>",      [](){testUpdate<short,short,false>();});
      register_test("hist.update<f,s,F>",      [](){testUpdate<float,short,false>();});
      register_test("hist.update<f,f,F>",      [](){testUpdate<float,float,false>();});
      // StorageT
      register_test("hist.clipped<s>",         [](){testClipped<short>();});
      register_test("hist.clipped<f>",         [](){testClipped<float>();});
      // StorageT
      register_test("hist.edges<b>",           testEdgesSByte);
      register_test("hist.edges<s>",           testEdgesShort);
      register_test("hist.edges<f>",           testEdgesFloat);
      // Miscellaneous
      register_test("hist.allnan",             testNaN);
      register_test("hist.allspecial.reuse",   [](){testSpecial<true>();});
      register_test("hist.allspecial.reopen",  [](){testSpecial<false>();});
      register_test("hist.forcerange.init",    testForceRangeInit);
      register_test("hist.forcerange.update",  testForceRangeUpdate);
      register_test("hist.shfromwriter_f",    [](){testStatHistFromWriter(true);});
      register_test("hist.shfromwriter",      [](){testStatHistFromWriter(false);});
      register_test("hist.annotfromwriter",   testAnnotFromWriter);
    }
  } dummy;
} // namespace for registration

#endif // HAVE_EXPANDABLE_BUILDER
