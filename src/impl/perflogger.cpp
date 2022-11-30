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

#include "file_performance.h"
#include "perflogger.h"
#include "timer.h"
#include "environment.h"
#include "exception.h"
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <atomic>
#include <algorithm>

namespace InternalZGY {
#if 0
}
#endif

std::atomic<int> PerformanceLogger::_last_id{0};

PerformanceLogger::PerformanceLogger(const std::string& outname, std::int64_t chunksize, int hist_bincount, double hist_min, double hist_max, int interval, const std::string& srcname)
  : _outfile()
  , _chunksize(chunksize)
  , _mutex()
  , _nsamples(0)
  , _histmin(hist_min)
  , _histmax(hist_max)
  , _statmin(std::numeric_limits<double>::max())
  , _statmax(std::numeric_limits<double>::lowest())
  , _statsum(0)
  , _statssq(0)
  , _histbins(hist_bincount, 0)
  , _sumtimer(new SummaryTimer("summary"))
  , _sumtimerbeg(std::numeric_limits<double>::max())
  , _sumtimerend(std::numeric_limits<double>::lowest())
  , _suminterval(interval)
  , _sumbytes(0)
  , _first(true)
  , _id(0)
  , _srcname(srcname)
{
  _id = 1 + _last_id.fetch_add(1);
  std::string name(outname);
  if (!name.empty()) {
    std::size_t pos = name.find("{}");
    if (pos != std::string::npos)
      name = name.substr(0, pos) + std::to_string(_id) + name.substr(pos+2);
    _outfile = std::make_shared<std::ofstream>(name, std::ofstream::app);
  }
  else {
    _outfile = std::shared_ptr<std::ostream>(&std::cerr, [](std::ostream*){});
  }
}

PerformanceLogger::~PerformanceLogger()
{
}

bool
PerformanceLogger::logThisSize(std::int64_t size)
{
  return size == _chunksize || _chunksize < 0;
}

void
PerformanceLogger::add(const Timer& timer, std::int64_t blocksize)
{
  std::unique_lock<std::mutex> lk(_mutex);
  // Optional periodic reporting.
  const double now = timer.getLastStop() / timer.getFrequency();
  if (_sumtimerbeg < _sumtimerend &&
      _suminterval > 0 &&
      now >= _sumtimerbeg + _suminterval)
  {
    // The timer passed in belongs to the next interval.
    lk.unlock();
    std::string msg = dumpThroughput(true);
    lk.lock();
    if (_outfile && !msg.empty())
      *_outfile << msg << std::flush;
    // Might also have reported and cleared the latency log.
  }

  // Update throughput logging.
  _sumtimer->add(timer);
  _sumtimerend = now;
  // ASSERT timer.getCount() == 1 otherwise we don't know the first start.
  // class Timer would need to keep track of first_start.
  _sumtimerbeg = std::min(_sumtimerbeg, _sumtimerend - timer.getTotal());
  _sumbytes += blocksize;

  // Update latency a.k.a. round trip time logging.
  // Note _hist{min,max} is the center of first and last bin.
  // Example: min=0, max=1000, nbuckets = 101
  // gives binwidth=10, bin(3) = 0, bin(7)=1
  // which is correct because bin 0 covers 0 +/- 5.
  const double value = timer.getTotal()*1000;
  const int    nbuckets = static_cast<int>(_histbins.size());
  const double binwidth = (_histmax - _histmin) / (nbuckets - 1);
  const double fbin = (value - _histmin) / binwidth;
  const int    ibin = (fbin < 0 ? 0 :
                       fbin > (nbuckets-1) ? nbuckets-1 :
                       static_cast<int>(std::floor(fbin+0.5)));
  _histbins[ibin] += 1;
  _nsamples += 1;
  _statsum += value;
  _statssq += value*value;
  if (_statmin > value) _statmin = value;
  if (_statmax < value) _statmax = value;
}

/**
 * Format the results as CSV that is simple to paste into a spreadsheet.
 * Use a distinct start/end marker so the csv data can be easily
 * extracted from a log file. Even if the previous line did not end in
 * a newline.
 *
 * Quick instructions for the spreadsheet:
 * Select the CSV3,CSV4,CSV5 cells except for the first and last columns.
 * Make a line graph: data series in rows, first row and first column as label.
 * This shows distibution of latency as well as cumulative counts; the latter
 * can be used to get the percentile of samples below a certain latency.
 */
std::string
PerformanceLogger::dumpLatency(bool clear)
{
  std::lock_guard<std::mutex> lk(_mutex);
  std::stringstream ss;
  if (_nsamples != 0) {
    const int    nbins = static_cast<int>(_histbins.size());
    const double binwidth = (_histmax - _histmin) / (nbins - 1);

    ss << "CSV1,ID,Samplecount,Histogram min,Histogram max"
       << ",Statistic min,Statistic max,Statistic average,Filename,END\n"
       << "CSV2," << _id << "," << _nsamples
       << "," << _histmin << "," << _histmax
       << "," << _statmin << "," << _statmax
       << "," << _statsum / _nsamples
       << ",\"" << _srcname << "\""
       << ",END\n";

    // Dump center value of each histogram bin.
    // Note that CSV3, CSV4, CSV5 can all be computed in a
    // spreadsheet using a simple formula. But I'd like to
    // have an (almost) single-click way of making the graph.
    ss << "CSV3," << _id << ",Latency";
    for (int ii=0; ii<nbins; ++ii)
      ss << "," << _histmin + binwidth * ii;
    ss << ",END\n";

    // Dump the bins as percentages
    ss << "CSV4," << _id << ",Frequency%";
    for (const std::int64_t count : _histbins)
      ss << "," << count / (double)_nsamples;
    ss << ",END\n";

    // Dump the cumulative counts.
    ss << "CSV5," << _id << ",Cumulative%";
    std::int64_t summed = 0;
    for (const std::int64_t count : _histbins) {
      summed += count;
      ss << "," << summed / (double)_nsamples;
    }
    ss << ",END\n";

    // Dump the bins as raw counts
    ss << "CSV6," << _id << ",Count";
    for (const std::int64_t count : _histbins)
      ss << "," << count;
    ss << ",END\n";

    if (clear) {
      _nsamples = 0;
      _statmin = std::numeric_limits<double>::max();
      _statmax = std::numeric_limits<double>::lowest();
      _statsum= 0;
      _statssq = 0;
      _histbins = std::vector<std::int64_t>(_histbins.size(), 0);
    }
  }
  return ss.str();
}

std::string
PerformanceLogger::dumpThroughput(bool clear)
{
  std::lock_guard<std::mutex> lk(_mutex);
  std::stringstream ss;
  if (_first) {
    ss << "CSV0,Timestamp,ID,Data(MB),Time(sec),Speed(MB/s),Readcount,Mean latency(ms),Filename,END\n";
    _first = false;
  }
  if (_sumtimerbeg < _sumtimerend && _sumbytes > 0) {
    const double bytecount = static_cast<double>(_sumbytes);
    const double elapsed = _sumtimerend - _sumtimerbeg;
    // Note, if you want to emulate an interval timer running
    // regardless of whether there was any traffic then round
    // begin and end to a multiple of interval. And add code
    // to output a lot of zeros for intervals with no traffic.
    ss << "CSV8"
       << std::setprecision(3) << std::fixed
       << "," << _sumtimerbeg
       << "," << _id
       << "," << bytecount / (1024.0*1024.0)
       << "," << elapsed
       << "," << (bytecount / (1024.0*1024.0)) / elapsed
       << "," << _sumtimer->getCount()
       << "," << std::setprecision(0)
       << (_sumtimer->getCount() == 0 ? 0 :
           1000.0 * _sumtimer->getTotal() / _sumtimer->getCount())
       << ",\"" << _srcname << "\""
       << ",END\n";
  }
  if (clear) {
    _sumtimer->reset();
    _sumtimerbeg = std::numeric_limits<double>::max();
    _sumtimerend = std::numeric_limits<double>::lowest();
    _sumbytes = 0;
  }
  return ss.str();
}

void
PerformanceLogger::dumpToFile(const std::string& comment)
{
  std::string str1 = dumpThroughput(true);
  std::string str2 = dumpLatency(true);
  std::string str3 = "CSV9," + std::to_string(_id) + "," + comment + ",END\n";
  if (_outfile && (!str1.empty() || !str2.empty()))
    *_outfile << str1 + str2 + str3 << std::flush;
}

} // namespace
