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

#include "timer.h"

#ifdef _WIN32
#include <Windows.h>
#endif

#include <sstream>
#include <stdarg.h>
#include <errno.h>
#ifndef _WIN32
#include <sys/time.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <atomic>
#include <omp.h>

namespace {
  int sprintf_s(char *buffer, int size, const char *format, ...)
  {
    if (size <= 0) {
      errno = EINVAL;
      return -1;
    }
    va_list ap;
    va_start(ap, format);
    int n = vsnprintf(buffer, size, format, ap);
    va_end(ap);
    if (n < 0 || n >= size) {
      buffer[0] = '\0'; // error or overflow. Don't write anything at all.
      errno = EINVAL;
      return -1;
    }
    return n; // Number of characters written, excluding terminatung null.
  }
}

namespace InternalZGY {
#if 0
}
#endif

/**
 * High resolution timer for performance measurement.
 *
 * The timer starts running as soon as it is constructed,
 * but for added accuracy it is recommended to invoke start() explicity.
 * The previous start time is then ignored.
 *
 * After calling stop(), the formatted elapsed time is available in getValue().
 * start() and stop() may be called multiple times to accumulate statistics.
 * In this case, it might be useful to pass a "skip" argument to the constructor
 * telling it that the first N laps are not representative and should be ignored.
 *
 * The resolution and accuracy for timing short duration calls is approx.
 * +/- 1 microsecond. Note that the class tries to adjust for the overhead
 * of calling the Windows performance timer.  This is done heuristically,
 * so if you time something that takes just a few nanoseconds to execute then
 * the result could end up negative.
 *
 * Applications may use Timer directly, but the higher level NamedTimer
 * will often be more convenient. The following examples show the use:
 *
 * \code
 * // One-shot timer:
 *    static bool enable = Timer.getNumericEnv("SALMON_TIMERS")>0;
 *    Timer t(enable);
 *    DoTheWork();
 *    t.Stop();
 *    printf("It took %s\n", t.getValue());
 *
 * //Statistics timer:
 *    static Timer t(enable);
 *    DoTheWork();
 *    t.Stop();
 *    if (last_call)
 *      printf("It took %s\n", t.getValue());
 * \endcode
 *
 * The first example is often too noisy, and there is no way for the
 * user to configure the timer to only print statistics.
 * The second approach may be tricky because it is not always obvious
 * when to print the gathered statistics. The higher level NamedTimer
 * and TimerPool solves these problems.
 */


/**
 * Create a new Timer instance, optionally giving it a name.
 * Optionally pass enabled=false to create a low overhead stub.
 * Optionally pass an initial count not to be included in statistics.
 * Once the instance is created the enabled state cannot be changed.
 */
Timer::Timer(bool enabled, const char* name, int skip, bool startrunning)
  : enabled_(enabled), skip_(skip), frequency_(1), overhead_(0),
    laps_(0), last_(0), total_(0), adjusted_(0), begin_(0), end_(0), running_(false), verbose_(enabled ? 1 : 0)
{
  name_[0] = '\0';
  buff_[0] = '\0';
  if (enabled_) {
    if (name != 0) {
#ifdef _WIN32
      strncpy_s(name_, name, sizeof(name_) - 1);
#else
      strncpy(name_, name, sizeof(name_)-1);
#endif
      name_[sizeof(name_)-1] = '\0';
    }
    frequency_ = getNativeFrequency();
    frequency_ = getNativeFrequency();
    if (frequency_ == 0) frequency_ = 1; // prevent divide-by-zero

    // Start and stop a couple of times to remove artifacts from JITting and cache
    start();
    stop();
    start();
    stop();
    start();
    stop();

    // The time for the last lap should be 0, but will probably be a small
    // number of microseconds representing the overhead of these calls.
    // Fudge subsequent numbers by 90% of this amount.
    overhead_ += (last_ * 9) / 10;

    // Reset the counters again, and then start them for real.
    // If the application does a start also, this first one is ignored.
    doReset();
    if (startrunning)
      doStart();
  }
}


long long
Timer::getNativeTime() {
#if 0
  // Not synchronized to UTC (bad), not affected by clock adjustments (good).
  LARGE_INTEGER now;
  QueryPerformanceCounter(&now);
  return now.QuadPart;
#elif defined(_WIN32)
  FILETIME ft{};
  GetSystemTimePreciseAsFileTime(&ft);
  return ((static_cast<long long>(ft.dwHighDateTime) << 32)
          + ft.dwLowDateTime
          - 116444736000000000LL);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (long long)tv.tv_sec * 1000000 + tv.tv_usec;
#endif
}


long long
Timer::getNativeFrequency() {
#if 0
  LARGE_INTEGER freq;
  QueryPerformanceFrequency(&freq);
  return freq.QuadPart;
#elif defined(_WIN32)
  return 10000000;
#else
  return 1000000;
#endif
}


/**
 * Start the timer. If already running, discard the previous start time.
 */
void
Timer::doStart()
{
  if (enabled_) {
                running_ = true;
                begin_ = getNativeTime();
  }
}


/**
 * Stop the timer. Add the elapsed time to the accumulated result.
 */
void
Timer::doStop()
{
  if (enabled_) {
    end_ = getNativeTime();
    if (running_) {
      last_ = end_ - begin_ - overhead_;
      if (laps_ >= skip_) {
        total_ += last_;
        adjusted_ += last_ / omp_get_num_threads();
      }
      ++laps_;
      running_ = false;
    }
  }
}


/**
 * Reset all counters, as if the timer instance was deleted and re-created.
 * The timer will be in a stopped state after this.
 */
void
Timer::doReset()
{
  last_    = 0;
  total_   = 0;
  adjusted_= 0;
  laps_    = 0;
  running_ = false;
}


/**
 * Show the elapsed time and lap count as a human readable string. The
 * function only uses public information, so callers need not use this
 * method if they don't like the particular formatting being done.
 *
 * The result of getValue() points to a class member that is
 * overwritten at each call. So this method is definitely not thread
 * safe. But safe to call in two different instances at the same time.
 * The result from getValue_s() is thread safe because caller supplies
 * the buffer,
 *
 * By default the result is shown in a "reasonable" unit. us, ms, s,
 * etc. depending on the elapsed time. This can be confusing if
 * printing a list of timers wanting to quickly see where the time is
 * being spent. Pass msonly=true for those cases.
 *
 * \param details If true, include the timer name.
 * \param msonly If true, always show result in milliseconds.
 */
const char*
Timer::getValue(bool details, bool msonly)
{
  if (!enabled_)
    return "";
  getValue_s(buff_, sizeof(buff_),
             getName(), getCount(), getTotal(), getAdjusted(), getRunning(),
             details, msonly);
  return buff_;
}


/**
 * \copydoc Timer::getValue
 */
void
Timer::getValue_s(char *result, int buffer_size, const char *name, int count, double total, double adjusted, bool running, bool details, bool msonly)
{
  result[0] = '\0';

  if (details) {
    int n = sprintf_s(result, buffer_size, "Time for %-16s: ", name[0] ? name : "(anonymous)");
    if (n > 0 && n < buffer_size) {
      buffer_size -= n;
      result += n;
    }
  }

  if (count == 0)
    sprintf_s(result, buffer_size, "No calls logged");
  else if (msonly)
    sprintf_s(result, buffer_size, "%9.1lf ms", total * 1000.0);
  else if (total < 0.00999)
    sprintf_s(result, buffer_size, "%6.1lf us", total * 1000000.0);
  else if (total < 9.99)
    sprintf_s(result, buffer_size, "%6.1lf ms", total * 1000.0);
  else
    sprintf_s(result, buffer_size, "%6.1lf s ", total * 1.0);

  if (count > 1)
    sprintf_s(result+strlen(result), buffer_size-strlen(result), " in %d calls", count);

  if (adjusted != total && adjusted != 0) {
    if (msonly)
      sprintf_s(result+strlen(result), buffer_size-strlen(result), " adjusted %9.1lf ms", adjusted * 1000.0);
    else if (total < 0.00999)
      sprintf_s(result+strlen(result), buffer_size-strlen(result), " adjusted %6.1lf us", adjusted * 1000000.0);
    else if (total < 9.99)
      sprintf_s(result+strlen(result), buffer_size-strlen(result), " adjusted %6.1lf ms", adjusted * 1000.0);
    else
      sprintf_s(result+strlen(result), buffer_size-strlen(result), " adjusted %6.1lf s ", adjusted * 1.0);
  }

  //if (getSkip() != 0)
  //  sprintf_s(result+strlen(result), buffer_size-strlen(result), " (first %d not counted)", getSkip());

  if (running)
    sprintf_s(result+strlen(result), buffer_size-strlen(result), " STILL RUNNING (forgot to Stop?)");

  if (details)
    sprintf_s(result+strlen(result), buffer_size-strlen(result), "\n");
}

/*=========================================================================*/
/*===   PrintingTimer   ===================================================*/
/*=========================================================================*/

PrintingTimer::PrintingTimer(const char *name, int level, bool startrunning)
  : Timer(true, name, 0, startrunning), level_(level)
{
}

PrintingTimer::~PrintingTimer()
{
  print();
}

void
PrintingTimer::print()
{
  if (getRunning())
    stop();
  // Most code will use one of the fancier timers that can redirect
  // output to a file or possibly interface with a proper logger.
  // What you see here is a bulletproof variant for ad-hoc debugging,
  // unit tests, etc. made to work also in a static destructor.
  if (getCount() != 0) {
    const char *msg = getValue(true, true);
    fwrite(msg, 1, strlen(msg), stderr);
  }
  reset();
}

/*=========================================================================*/
/*   SimpleTimer   ========================================================*/
/*=========================================================================*/

// Inline code only

/*=========================================================================*/
/*   SummaryTimer   =======================================================*/
/*=========================================================================*/

/*
 * \details Thread safety:
 * Modification other than the statistics may lead to a data race.
 * Planned usage is that most information will be only set in create.
 * Updating statistics, as is done in SummaryTimer::add() is thread
 * safe. Note in particular the buff_ member which is not safe and
 * which causes SummaryTimer::getValue() to not be threads safe.
 */
class SummaryTimer::Impl
{
public:
  long long frequency_;
  std::atomic<int> count_;
  std::atomic<long long> total_;
  std::atomic<long long> adjusted_;
  std::atomic<long long> last_;
  char name_[256];    // Optional name of this timer
  char buff_[256];    // Returned from getValue().
  explicit Impl(const char *name)
    : frequency_(1000*1000), count_(0), total_(0), adjusted_(0), last_(0)
  {
    name_[0] = '\0';
    buff_[0] = '\0';
    if (name != 0) {
#ifdef _WIN32
      strncpy_s(name_, name, sizeof(name_) - 1);
#else
      strncpy(name_, name, sizeof(name_)-1);
#endif
      name_[sizeof(name_)-1] = '\0';
    }
  }
};

SummaryTimer::SummaryTimer(const char *name)
  : pimpl_(new Impl(name))
{
}

SummaryTimer::~SummaryTimer()
{
  delete pimpl_;
  pimpl_ =  nullptr;
}

double
SummaryTimer::getFrequency() const
{
  // Reason the cast is safe: Even if there is a loss of precision,
  // this would only cause a slightly less precise measurement.
  return static_cast<double>(pimpl_->frequency_);
}

int
SummaryTimer::getCount() const
{
  return pimpl_->count_.load();
}

double
SummaryTimer::getTotal() const
{
  return static_cast<double>(pimpl_->total_.load()) / pimpl_->frequency_;
}

double
SummaryTimer::getAdjusted() const
{
  return static_cast<double>(pimpl_->adjusted_.load()) / pimpl_->frequency_;
}

double
SummaryTimer::getLast() const
{
  return static_cast<double>(pimpl_->last_.load())  / pimpl_->frequency_;
}

const char*
SummaryTimer::getName() const
{
  return pimpl_->name_;
}

const char*
SummaryTimer::getValue(bool details, bool msonly) const
{
  Timer::getValue_s(pimpl_->buff_, sizeof(pimpl_->buff_),
                    getName(), getCount(), getTotal(), getAdjusted(), false,
                    details, msonly);
  return pimpl_->buff_;
}

const char*
SummaryTimer::getCSV() const
{
  // Avoid odd characters in name to make it easier for spreadsheets.
  char *safename = new char[strlen(getName())+1];
  char *dst = safename;
  const char *src = getName(); // Instance member, remains valid long enough.
  for (; *src; ++src, ++dst)
    *dst = isalnum(*src) ? *src : '_';
  *dst = '\0';
  sprintf_s(pimpl_->buff_, sizeof(pimpl_->buff_),
            "TIMER,\"%s\",%d,%.3lf,%.3lf\n",
            safename, getCount(), getTotal(), getAdjusted());
  delete[] safename;
  return pimpl_->buff_;
}

/**
 * \brief Clear all counters.
 * \details The timer name is unchanged.
 */
void
SummaryTimer::reset()
{
  pimpl_->count_ = 0;
  pimpl_->total_ = 0;
  pimpl_->adjusted_  = 0;
  pimpl_->last_  = 0;
}

/**
 * \brief Add the contents of the specified Timer to our summary.
 * \details The information to add is passed as discrete arguments,
 * so it is possible to e.g. add an managed timer to an umanaged one.
 */
void
SummaryTimer::add(int count, double total, double adjusted, double last)
{
  pimpl_->count_.fetch_add(count);
  pimpl_->total_.fetch_add((long long)(total * pimpl_->frequency_));
  pimpl_->adjusted_.fetch_add((long long)(adjusted * pimpl_->frequency_));
  pimpl_->last_.store((long long)(last * pimpl_->frequency_));
}

/**
 * \brief Add the contents of the specified Timer to our summary.
 * \details The caller is responsible for stopping the timer first.
 */
void
SummaryTimer::add(const Timer& t)
{
  if (t.getEnabled()) {
    add(t.getCount(), t.getTotal(), t.getAdjusted(), t.getLast());
  }
}

/*=========================================================================*/
/*   SummaryPrintingTimer   ===============================================*/
/*=========================================================================*/

SummaryPrintingTimer::SummaryPrintingTimer(const char *name, bool csv)
  : SummaryTimer(name)
  , csv_(csv)
{
}

SummaryPrintingTimer::~SummaryPrintingTimer()
{
  print();
}

void
SummaryPrintingTimer::print()
{
  if (getCount() != 0) {
    const char *msg = csv_ ? getCSV() : getValue(true, true);
    fwrite(msg, 1, strlen(msg), stderr);
  }
  reset();
}

} // end namespace
