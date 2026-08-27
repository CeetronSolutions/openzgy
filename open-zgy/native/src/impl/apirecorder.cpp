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

#include "apirecorder.h"
#include "environment.h"

#include <fstream>
#include <chrono>

// Ugh! All this to get the thread id.
#ifdef _WIN32
#define WINDOWS_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#endif

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY { namespace Telemetry {
#if 0
}}
#endif

namespace {
  // gcc 9 and above have gettid() in the system headers, but
  // devtoolset-10 doesn't. Use a different name to avoid conflicts.
#if defined(_WIN32)
  static int my_gettid() { return ::GetCurrentThreadId(); }
  static int my_getpid() { return ::GetCurrentProcessId(); }
#elif defined(__GNUC__) && defined(SYS_gettid)
  static int my_gettid() { return syscall(SYS_gettid); }
  static int my_getpid() { return getpid(); }
#else
  // If the toolset has both SYS_gettid and gettid(), the former
  // will be used. No big deal, I think.
  static int my_gettid() { return gettid(); }
  static int my_getpid() { return getpid(); }
#endif
}

/////////////////////////////////////////////////////////////////////////////
///   ApiRecorder - do the actual logging                                 ///
/////////////////////////////////////////////////////////////////////////////

std::atomic<int> ApiRecorder::next_id_(1);
std::mutex ApiRecorder::global_mutex_;
std::weak_ptr<std::ostream> ApiRecorder::global_os_;

ApiRecorder::ApiRecorder(const std::string& name)
  : id_(next_id_.fetch_add(1))
  , name_(name)
  , os_(getVerboseFileFromEnv(id_))
{
  logger("ApiRecorder\topen\t\"" + name_ + "\"");
}

ApiRecorder::~ApiRecorder()
{
  logger("ApiRecorder\tdestruct");
  if (os_.use_count() == 1)
    logger("ApiRecorder\tcloselog");
}

void
ApiRecorder::logger(const std::string& str) const
{
  std::stringstream ss;
  ss.precision(3);
  ss << id_
     << "\t" << my_getpid()
     << "\t" << my_gettid()
     << "\t" << std::fixed << timestamp()
     << "\t" << str
     << "\n";
  std::lock_guard<std::mutex> lk(global_mutex_);
  if (os_ && os_->good())
    *os_ << ss.str();
}

void
ApiRecorder::logger(const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  logger(sstream ? sstream->str() : std::string());
}

void
ApiRecorder::logio(
     const std::string prefix,
     const size3i_t& start,
     const size3i_t& size,
     int lod) const
{
  logger(std::stringstream().flush()
         << prefix
         << "\t" << start[0] << "\t" << start[1] << "\t" << start[2]
         << "\t" << size[0]  << "\t" << size[1]  << "\t" << size[2]
         << "\t" << lod);
}

bool
ApiRecorder::enabled()
{
  static std::string env = InternalZGY::Environment::getStringEnv("OPENZGY_RECORD_LOGFILE");
  return !env.empty();
}

std::shared_ptr<std::ostream>
ApiRecorder::getVerboseFileFromEnv(int id)
{
  std::lock_guard<std::mutex> lk(global_mutex_);
  std::string name = InternalZGY::Environment::getStringEnv("OPENZGY_RECORD_LOGFILE");
  if (name.empty() || name == "cerr" || name == "stderr") {
    return std::shared_ptr<std::ostream>(&std::cerr, [](std::ostream*){});
    // Never closed. The smart pointer has a dummy destructor.
  }
  else if (name == "cout" || name == "stdout" || name == "con:") {
    return std::shared_ptr<std::ostream>(&std::cout, [](std::ostream*){});
  }
  else {
    std::size_t pos = name.find("{pid}");
    if (pos != std::string::npos) {
      name = name.substr(0, pos) + std::to_string(my_getpid()) + name.substr(pos+5);
    }
    pos = name.find("{id}");
    if (pos != std::string::npos) {
      name = name.substr(0, pos) + std::to_string(id) + name.substr(pos+4);
      return std::make_shared<std::ofstream>(name, std::ofstream::app);
      // One file per ApiRecorder instance, but not per thread.
    }
    else {
      std::shared_ptr<std::ostream> os = global_os_.lock();
      if (!os) {
        os = std::make_shared<std::ofstream>(name, std::ofstream::app);
        global_os_ = os;
        // Weak pointer. The file is closed when the last ApiRecorder exits.
      }
      return os;
    }
  }
}

double
ApiRecorder::timestamp()
{
  // Remember the good old days with ::time(nullptr) ?
  typedef std::chrono::high_resolution_clock clock;
  static constexpr double factor =
    static_cast<double>(clock::duration::period::num) /
    static_cast<double>(clock::duration::period::den);
  return clock::now().time_since_epoch().count() * factor;
}

/////////////////////////////////////////////////////////////////////////////
///   ApiReaderRecorder - intercept usage.                                ///
/////////////////////////////////////////////////////////////////////////////

ApiReaderRecorder::ApiReaderRecorder(
     const std::shared_ptr<IZgyReader>& relay,
     const std::string& name)
  : ZgyReaderRelay(relay)
  , recorder_(std::make_shared<ApiRecorder>(name))
{
}

ApiReaderRecorder::~ApiReaderRecorder()
{
}

std::shared_ptr<IZgyReader>
ApiReaderRecorder::inject(
     std::shared_ptr<IZgyReader> in,
     const std::string& name)
{
  if (ApiRecorder::enabled())
    in.reset(new ApiReaderRecorder(in, name));
  return in;
}

void ApiReaderRecorder::read(const size3i_t& start, const size3i_t& size, float* data, int lod) const
{
  logio("ZgyReader\tread\tfloat", start, size, lod);
  relay().read(start, size, data, lod);
}

void ApiReaderRecorder::read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod) const
{
  logio("ZgyReader\tread\tint16", start, size, lod);
  relay().read(start, size, data, lod);
}

void ApiReaderRecorder::read(const size3i_t& start, const size3i_t& size, std::int8_t* data, int lod) const
{
  logio("ZgyReader\tread\tint8", start, size, lod);
  relay().read(start, size, data, lod);
}

std::pair<bool,double> ApiReaderRecorder::readconst(const size3i_t& start, const size3i_t& size, int lod, bool as_float = true) const
{
  logio("ZgyReader\treadconst\t*", start, size, lod);
  return relay().readconst(start, size, lod, as_float);
}

void ApiReaderRecorder::close()
{
  relay().close();
}

/////////////////////////////////////////////////////////////////////////////
///   ApiWriterRecorder - intercept usage.                                ///
/////////////////////////////////////////////////////////////////////////////

ApiWriterRecorder::ApiWriterRecorder(
     const std::shared_ptr<IZgyWriter>& relay,
     const std::string& name)
  : ZgyWriterRelay(relay)
  , recorder_(std::make_shared<ApiRecorder>(name))
{
}

ApiWriterRecorder::~ApiWriterRecorder()
{
}

std::shared_ptr<IZgyWriter>
ApiWriterRecorder::inject(
     std::shared_ptr<IZgyWriter> in,
     const std::string& name)
{
  if (ApiRecorder::enabled())
    in.reset(new ApiWriterRecorder(in, name));
  return in;
}

void ApiWriterRecorder::read(const size3i_t& start, const size3i_t& size, float* data) const
{
  logio("ZgyWriter\tread\tfloat", start, size, 0);
  relay().read(start, size, data);
}

void ApiWriterRecorder::read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const
{
  logio("ZgyWriter\tread\tint16", start, size, 0);
  relay().read(start, size, data);
}

void ApiWriterRecorder::read(const size3i_t& start, const size3i_t& size, std::int8_t* data) const
{
  logio("ZgyWriter\tread\tint8", start, size, 0);
  relay().read(start, size, data);
}

std::pair<bool,double> ApiWriterRecorder::readconst(const size3i_t& start, const size3i_t& size, bool as_float) const
{
  logio("ZgyWriter\treadconst\t*", start, size, 0);
  return relay().readconst(start, size, as_float);
}

void ApiWriterRecorder::write(const size3i_t& start, const size3i_t& size, const float* data)
{
  logio("ZgyWriter\twrite\tfloat", start, size, 0);
  relay().write(start, size, data);
}

void ApiWriterRecorder::write(const size3i_t& start, const size3i_t& size, const std::int16_t *data)
{
  logio("ZgyWriter\twrite\tint16", start, size, 0);
  relay().write(start, size, data);
}

void ApiWriterRecorder::write(const size3i_t& start, const size3i_t& size, const std::int8_t* data)
{
  logio("ZgyWriter\twrite\tint8", start, size, 0);
  relay().write(start, size, data);
}

void ApiWriterRecorder::writeconst(const size3i_t& start, const size3i_t& size, const float* data)
{
  logio("ZgyWriter\twriteconst\tfloat", start, size, 0);
  relay().writeconst(start, size, data);
}

void ApiWriterRecorder::writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data)
{
  logio("ZgyWriter\twriteconst\tint16", start, size, 0);
  relay().writeconst(start, size, data);
}

void ApiWriterRecorder::writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data)
{
  logio("ZgyWriter\twriteconst\tint8", start, size, 0);
  relay().writeconst(start, size, data);
}

void ApiWriterRecorder::finalize(
     const std::vector<DecimationType>& decimation,
     const std::function<bool(std::int64_t,std::int64_t)>& progress,
     FinalizeAction action,
     bool force)
{
  logger("ZgyWriter\tfinalize");
  relay().finalize(decimation, progress, action, force);
}

void ApiWriterRecorder::close_incomplete()
{
  logger("ZgyWriter\tclose_incomplete");
  relay().close_incomplete();
}

void ApiWriterRecorder::close()
{
  logger("ZgyWriter\tclose");
  relay().close();
}

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
