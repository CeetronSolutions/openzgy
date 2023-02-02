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

#pragma once

#include "declspec.h"

namespace InternalZGY {
#if 0
}
#endif

/**
 * \file timer.h
 * \brief Easy to use performance measurement.
 *
 * The Timer class is implemented using basic, old style techniques.
 * The goal is to make it as safe as possible, usable also if the stop
 * & print happens in a static destructor where stdio and (if using
 * .NET) the CLR might be shut down.
 *
 * It is also a goal to avoid as many includes as possible in the
 * header file, so as to not introduce additional dependencies. There
 * are a few enhancements that have been rejected because of this. I
 * might need to revist those, though.
 *
 * The SummaryTimer and derived classes loosen these requirements
 * somewhat. A pimpl is used (since those classes are less likely
 * to be performance critical) and atomics are used to make add()
 * threadsafe.
 *
 * Thread safety: Not thread safe. Not a problem if used as recommended.
 * See comments in SimpleTimer.
 *
 * \li Allowing \#include\<functional\> in the header would make PrintingTimer more general.
 *     Or maybe template PrintingTimer on the functor that does the printing
 *     and the data type of the counters (to hide atomic).
 *     typedef PrintingTimerT<std::atomic<long long>, std::functor<const std::string&>> PrintingTimer
 * \li Allowing \#include\<chrono\> would make the code more portable but requiring c++11.
 *     Note that the actual timers should probably be long long microseconds
       since epoch so the chrono include isn't visible in the header.
 */

class OPENZGY_TEST_API Timer {

private: // readonly, i.e. set only in constructor
  bool enabled_;      // If false, the class is just a stub.
  int  skip_;         // Do not accumulate the N first results
  long long frequency_; // Frequency of internal timer, in Hertz
  long long overhead_;  // Estimated overhead of reading timer, in usec
  char name_[256];    // Optional name of this timer
  char buff_[256];    // Returned from getValue().

private: // mutable
  int  laps_;         // Number of accumulated results
  long long last_;      // Time for last lap
  long long total_;     // Accumulated total time
  long long adjusted_;  // Accumulated total time adjusted for thread count
  long long begin_;     // Time of last call to Start();
  long long end_;       // Time of last call to Stop();
  bool running_;      // Guard against two Stop() in a row.
  int  verbose_;      // Used by higher levels only.

private: // non-inlined versions of operations
  void doStart();
  void doStop();
  void doReset();
  long long getNativeTime();
  long long getNativeFrequency();

public:
  Timer(bool enable = true, const char* name = 0, int skip = 0, bool startrunning = true);

  // ACCESSORS
  bool   getEnabled()    const { return enabled_; }
  double getFrequency()  const { return static_cast<double>(frequency_); }
  double getLast()       const { return static_cast<double>(last_)     / static_cast<double>(frequency_); }
  double getTotal()      const { return static_cast<double>(total_)    / static_cast<double>(frequency_); }
  double getAdjusted()   const { return static_cast<double>(adjusted_) / static_cast<double>(frequency_); }
  double getOverhead()   const { return static_cast<double>(overhead_) / static_cast<double>(frequency_); }
  int    getCount()      const { return laps_ < skip_ ? 0 : laps_ - skip_; }
  const char* getName()  const { return name_; }
  int    getSkip()       const { return skip_ < laps_ ? skip_ : laps_; }
  bool   getRunning()    const { return running_; }
  int    getVerbose()    const { return verbose_; }
  static void   getValue_s(char *buf, int len, const char *name, int count, double total, double adjusted, bool running, bool details, bool msonly);
  const char* getValue(bool details = false, bool msonly = false);
  // Shouldn't be public but might come in handy if extending Timer.
  double getLastStop()   const { return static_cast<double>(end_); }

  // OPERATIONS
  void setVerbose(int v) { verbose_ = v; }
  void start()           { if (enabled_) doStart(); }
  void stop()            { if (enabled_) doStop();  }
  void reset()           { if (enabled_) doReset(); }
};

/**
 * \brief Timer that prints its result when going out of scope.
 *
 * Convenience class that stops (if needed) the timer and logs the result
 * when the instance goes out of scope.
 *
 * Hopefully the logging is safe even when run at exit (e.g. in a static
 * destructor) when the higher level logging mechanism or even iostreams
 * might be shut down. Use one of the fancier timers if you need to
 * control where the output ends up.
 *
 * Thread safety:
 * Modification may lead to a data race. This should not be an issue,
 * because instances are meant to be short lived. Typically used inside
 * a method and not acceesible outside. If used for something else then
 * take care. Consider using a long lived SummaryPrintingTimer with a
 * short lived Timer to handle the actual measurements.
 */
class OPENZGY_TEST_API PrintingTimer : public Timer
{
  int level_;
public:
  explicit PrintingTimer(const char *name, int level = 1, bool startrunning = true);
  ~PrintingTimer();
  void print();
};

/**
 * \brief Hold the timing results from zero or more Timer instances.
 * \copydetails SimpleTimer
 */
class OPENZGY_TEST_API SummaryTimer
{
  class Impl;
  Impl *pimpl_;
public:
  explicit SummaryTimer(const char* name);
  virtual ~SummaryTimer();
  SummaryTimer(const SummaryTimer&) = delete;
  SummaryTimer& operator=(const SummaryTimer&) = delete;

  // ACCESSORS
  double getFrequency() const;
  int    getCount()     const;
  double getTotal()     const;
  double getAdjusted()  const;
  double getLast()      const;
  const char* getName() const;
  const char* getValue(bool details, bool msonly) const;
  const char* getCSV() const;

  // OPERATIONS
  void reset();
  void add(int count, double total, double adjusted, double last);
  void add(const Timer& t);
};

/**
 * \brief SummaryTimer that prints its result when going out of scope.
 * \copydetails SimpleTimer
 */
class OPENZGY_TEST_API SummaryPrintingTimer : public SummaryTimer
{
private:
  bool csv_;
public:
  explicit SummaryPrintingTimer(const char *name, bool csv = false);
  virtual ~SummaryPrintingTimer();
  virtual void print();
};

/**
 * \brief Timer that knows where to store the result.
 *
 * SimpleTimer is normally short lived. The SimpleTimer constructor
 * takes a SummaryTimer as an argument. When the SimpleTimer goes out
 * of scope the result is added to the accumulate timer. The caller is
 * responsible for not allowing the SummaryTimer to go out of scope
 * before the SimpleTimer does.
 *
 * There is also a SummaryPrintingTimer class that extends SummaryTimer
 * to print a single line report when it goes out of scope.
 *
 * Thread safety:
 * Note that Timer and SimpleTimer are not designed to be threadsafe.
 * This makes no sense when measuring a single piece of code execution.
 * SummaryTimer::add() is threadsafe. This means the following example
 * with the SummaryTimer being shared among threads will work.
 * The example  measures every execution of the function and prints a
 * report on application exit.
 *
 * \code
 * void SomeFunction() {
 *    static SummaryPrintingTimer pt("MyTimer");
 *    SimpleTimer tt(pt);
 *    // timing the following code...
 *}
 * \endcode
 *
 * A more obscure feature: If you want to output the results of a
 * SummaryPrintngTimer at a specific point then you can call done()
 * on any SimpleTimer that is linked to this timer and is still
 * in scope, and then call print() on the PrintingTimer to pretend
 * that it went out of scope.
 */
class OPENZGY_TEST_API SimpleTimer : public Timer
{
  SummaryTimer& owner_;
public:
  explicit SimpleTimer(SummaryTimer& owner, bool enabled = true)
    : Timer(enabled)
    , owner_(owner)
  {
  }
  ~SimpleTimer()
  {
    done();
  }
  void done() {
    if (getEnabled()) {
      stop();
      owner_.add(*this);
      reset();
    }
  }
};

} // end namespace
