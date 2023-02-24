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

/**
 * \file zgycopyc.cpp
 * \brief Copy one ZGY file to another using the OpenZGY/C++ API.
 *
 * The initial version of this program was based on "example.h"
 * I will need to add a lot more bells and whistles because I need
 * the program to measure performance. This is also why I cannot
 * simply use the Python API.
 * See src/Salmon/Shared/Tools/ZgyCopy/ZgyCopy.cpp in The Salmon baseline.
 */

#include "../api.h"
#include "../impl/fancy_timers.h"
#include "../test/mock.h"
#include "../impl/environment.h"
#include "../impl/mtguard.h"
#include "../iocontext.h"
#include "readwritemirror.h"
#include "readlodcrop.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <stdlib.h>
#include <chrono>
#include <random>
#include <malloc.h>
#include <omp.h>
#include <atomic>
#include <cstring>
#include <mutex>
#include <thread>
#include <chrono>
#ifndef _WIN32
#include <signal.h>
#endif

#ifndef _WIN32
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#define HAVE_GETOPT
#else
#include <io.h>
#undef HAVE_GETOPT
#endif

// Enable the --dropcache option for testing on-prem I/O.
// It is not normally present because execuring an external
// program will be flagged as a potential vulnerability.
// When enabling this, also change SYSTEM() to system().
// The vulnerability scanner apparently doesn't parse #if properly
// so it still sees the call to system() even when ifdeffed out.
//#define HAVE_DROPCACHE 1

/*
 * OpenMP note:
 * Most of the Linux configs I am currently building have OpenMP 2015.11.
 * Ubuntu Xenial was stuck at 2013.07, as was devtoolset-4 that I used when
 * building SDAPI for CentOS 7. The latest Windows compilers use 2002.03.
 * Linux check: g++ -fopenmp -dM -E -x c++ - < /dev/null | grep _OPENMP
 */

using InternalZGY::SummaryPrintingTimerEx;
using InternalZGY::SimpleTimerEx;
using InternalZGY::Environment;
using InternalZGY::MTGuardWithProgress;
using OpenZGY::SeismicStoreIOContext;
using OpenZGY::Tools::ZgyReaderMirror;
using OpenZGY::Tools::ZgyWriterMirror;
using OpenZGY::Tools::ZgyReadLodCrop;

/*=========================================================================*/
/*   OPTION PROCESSING   ==================================================*/
/*=========================================================================*/

#ifndef HAVE_GETOPT
  // getopt_long() missing from windows, so fake it.
  struct option {
    const char *name;
    int         has_arg;
    int        *flag;
    int         val;
  };
  enum { no_argument = 0, required_argument = 1, optional_argument = 2 };
#endif

/**
 * This is for token stored in a file, which is a kludge for testing.
 *
 * Assume the token that was written to the file had at least 40
 * minutes till expiry. If the file is more than 40 minutes old then
 * assume the token has expired. Actually checking the token would be
 * a serious pain. Wait for the file to be updated with a new token.
 * recheck every 30 seconds. Print a reminder every 5 minutes.
 *
 * The caller must *either* ensure that it holds a lock causing all
 * threads in this file descriptor (or alternatively in all files)
 * to block when this methods sleeps. *Or* uncomment the mutex in the
 * code below. Using both will increase the risk of deadlocks.
 * Especially if caller has a per-file mutes while this method has
 * a global one.
 *
 * Note that specifying a token file in a local to local copy is a
 * bad idea, because the code here will still demand that the file
 * should be updated regularly. Unless waitForFreshToken() is called
 * from inside the SDAPI callback.
 *
 * Caveats if calling from inside the SDAPI callback:
 *  - Higher risk of deadlocks because we cannot know what SDAPI is doing.
 *  - Risk of timeouts if SDAPI expects token callbacks to be quick.
 *  - Multiple scenarios are possible, depending on timing. All need testing.
 *     - Both the reader and the writer might be blocked waiting on the file.
 *     - The reader might block here, with the writer starving for missing data,
 *     - Probably other scenarios.
 *
 * The code assumes that a call to ::time() is a lot cheaper than
 * ::stat(), and will not check more often that every minute.
 * Except when waiting for the file to be updated.
 */
static void
waitForFreshToken(const std::string& token, std::int64_t *last_check)
{
#ifndef _WIN32
  static const auto isfile = [](const std::string& s) {
    return !s.empty() && (s.front() == '/' || s.front() == '\\');
  };
  if (isfile(token)) {
    const std::string filename = token;

    std::int64_t now = (std::int64_t)::time(nullptr);
    if (now - *last_check < 60)
      return;

    struct stat st{};
    for (int loops = 0;; ++loops) {
      if (::stat(filename.c_str(), &st) >= 0) {
        now = (std::int64_t)::time(nullptr);
        std::int64_t age = now - (std::int64_t)st.st_mtime;
        if (age < 40*60) {
          if (loops != 0)
          std::cerr << "Token in \"" << filename << "\""
                    << " Age is " << age / 60 << " minutes." << std::endl;
          *last_check = now;
          return;
        }
        if (loops == 0)
          std::cerr << "\n"; // There is probably a progress bar.
        if ((loops % 10) == 0)
          std::cerr << "Waiting for new token in \"" << filename << "\"."
                    << " Age is " << age / 60 << " minutes." << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(30));
      }
      else {
        // File not found. Treat this the same way as a too old file.
        if ((loops % 10) == 0)
          std::cerr << "Waiting for \"" << filename << "\""
                    << " to be cratated." << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(30));
      }
    }
  }
  #endif
}

static void
waitForTokenInTokenCallback(
     const std::string& token, std::int64_t *last_check)
{
  // NOTE: Only one of waitForTokenInTokenCallback and
  // waitForTokenInCopyLoop should be enabled.
  // NOTE, a mutex is held in SDTokenUpdater::operator(). None here.
  waitForFreshToken(token, last_check);
}

static void
waitForTokenInCopyLoop(
     const std::string& token, std::int64_t *last_check)
{
  //static std::mutex mutex;
  //std::lock_guard<std::mutex> lk(mutex);
  //waitForFreshToken(token, last_check);
}

/**
 * Static methods needed for option processing.
 * The entire class cam probably be re-used in other command line
 * applications that copies from some format to ZGY.
 */
class OptionsTools
{
public:
  typedef std::pair<OpenZGY::FinalizeAction, bool> finalize_t;

protected:
  OptionsTools();
  static int geti(const char *str);
  static double getd(const char *str);
  static std::vector<int> getlist(const char *str);
  static std::vector<double> getdoublelist(const char *str);
  static std::vector<std::pair<OpenZGY::DecimationType, const char*>> decimationTypes();
  static OpenZGY::DecimationType getDecimationType(const char *name);
  static std::vector<OpenZGY::DecimationType> getDecimationTypes(const char *names);
  static std::vector<std::pair<finalize_t, const char*>> finalizeTypes();
  static finalize_t getFinalizeAction(const char *name);
  static const char* showFinalizeAction(const finalize_t& finalize);
};

OptionsTools::OptionsTools()
{
}

int
OptionsTools::geti(const char *str)
{
  char *end;
  int result = static_cast<int>(strtol(str, &end, 10));
  if (end == str || *end != '\0')
    throw std::runtime_error("command line: expected a number, found \"" + std::string(str) + "\"");
  return result;
}

/**
 * Parse a list of integers separated by comma or x
 */
std::vector<int>
OptionsTools::getlist(const char *str)
{
  std::vector<int> result;
  if (str && *str) {
    for (;;) {
      char *end;
      result.push_back(static_cast<int>(strtol(str, &end, 10)));
      if (end == str)
        throw std::runtime_error("command line: expected a number, found \"" + std::string(str) + "\"");
      if (*end != '\0' && *end != 'x' && *end != ',')
        throw std::runtime_error("command line: numbers should be separated by comma or x, not '" + std::string(end).substr(0,1) + "'");
      if (*end == '\0')
        break;
      str = end + 1;
    }
  }
  return result;
}

double
OptionsTools::getd(const char *str)
{
  char *end;
  double result = (strtod(str, &end));
  if (end == str || *end != '\0')
    throw std::runtime_error("command line: expected a decimal number, found \"" + std::string(str) + "\"");
  return result;
}

/**
 * Parse a list of decimal numbers separated by comma or x
 */
std::vector<double>
OptionsTools::getdoublelist(const char *str)
{
  std::vector<double> result;
  if (str && *str) {
    for (;;) {
      char *end;
      result.push_back(strtod(str, &end));
      if (end == str)
        throw std::runtime_error("command line: expected a decimal number, found \"" + std::string(str) + "\"");
      if (*end != '\0' && *end != 'x' && *end != ',')
        throw std::runtime_error("command line: numbers should be separated by comma or x, not '" + std::string(end).substr(0,1) + "'");
      if (*end == '\0')
        break;
      str = end + 1;
    }
  }
  return result;
}

std::vector<std::pair<OpenZGY::DecimationType, const char*>>
OptionsTools::decimationTypes()
{
  using OpenZGY::DecimationType;
  static std::vector<std::pair<DecimationType, const char*>> result
    {
     {DecimationType::LowPass,"LowPass"},
     {DecimationType::LowPassNew,"LowPassNew"},
     {DecimationType::WeightedAverage,"WeightedAverage"},
     {DecimationType::Average,"Average"},
     {DecimationType::Median,"Median"},
     {DecimationType::Minimum,"Minimum"},
     {DecimationType::Maximum,"Maximum"},
     //{DecimationType::MinMax,"MinMax"},
     {DecimationType::Decimate,"Decimate"},
     {DecimationType::DecimateSkipNaN,"DecimateSkipNaN"},
     //{DecimationType::DecimateRandom,"DecimateRandom"},
     {DecimationType::AllZero,"AllZero"},
     //{DecimationType::WhiteNoise,"WhiteNoise"},
     {DecimationType::MostFrequent,"MostFrequent"},
     {DecimationType::MostFrequentNon0,"MostFrequentNon0"},
     {DecimationType::AverageNon0,"AverageNon0"},
    };
  return result;
}

OpenZGY::DecimationType
OptionsTools::getDecimationType(const char *name)
{
  for (const auto& it : decimationTypes())
    if (0==strcmp(name, it.second))
      return it.first;
  throw std::runtime_error("Invalid decimation type \"" +
                           std::string(name) + "\".");
}

std::vector<OpenZGY::DecimationType>
OptionsTools::getDecimationTypes(const char *names)
{
  std::vector<OpenZGY::DecimationType> result;
  std::string s(names);
  for (;;) {
    auto delim = s.find_first_of(",/ ");
    result.push_back(getDecimationType(s.substr(0, delim).c_str()));
    if (delim == std::string::npos)
      break;
    s = s.substr(delim+1);
  }
  return result;
}

std::vector<std::pair<OptionsTools::finalize_t, const char*>>
OptionsTools::finalizeTypes()
{
  using OpenZGY::FinalizeAction;
  static std::vector<std::pair<finalize_t, const char*>> result
    {
     {{FinalizeAction::BuildDefault, false},     "default"},
     {{FinalizeAction::Delete,       false},     "delete"},
     {{FinalizeAction::Delete,       true},      "force-delete"},
     {{FinalizeAction::Keep,         false},     "keep"},
     {{FinalizeAction::Keep,         true},      "force-keep"},
     {{FinalizeAction::BuildIncremental, false}, "incremental"},
     {{FinalizeAction::BuildIncremental, true},  "force-incremental"},
     {{FinalizeAction::BuildFull,        false}, "full"},
     {{FinalizeAction::BuildFull,        true},  "force-full"},
     {{FinalizeAction::BuildNoHistogram, false}, "nohistogram"},
     {{FinalizeAction::BuildNoHistogram, true},  "force-nohistogram"},
    };
  return result;
}

OptionsTools::finalize_t
OptionsTools::getFinalizeAction(const char *name)
{
#ifdef _WIN32
  for (const auto& it : finalizeTypes())
    if (0 == _stricmp(name, it.second))
      return it.first;
#else
  for (const auto& it : finalizeTypes())
    if (0==strcasecmp(name, it.second))
      return it.first;
#endif
  throw std::runtime_error("Invalid finalize type \"" +
                           std::string(name) + "\".");
}

const char*
OptionsTools::showFinalizeAction(const finalize_t& finalize)
{
  for (const auto& it : finalizeTypes())
    if (it.first == finalize)
      return it.second;
  return "(invalid)";
}

/**
 * TODO-Low: Move content here that might be re-used in similar apps.
 *
 * Specific to vds2zgy:
 *   --connection
 * Specific to zgycopyc, as the feature might be too expensive for vds2zgy.
 *   --native, --float
 * Specific to zgycopyc, used for testing OpenZGY itself.
 *   --alpha, --lod, --dumpsqnr (actually not implemented yet)
 *   --update, --size, --noise
 *   --lod, --cropstart, --cropsize, --mirror, (--upsample)
 * Specific to zgycopyc because it uses OpenMP.
 *   --omp-nest
 * Features specific to zgycopyc, no explicit option.
 *   (Allow updating token during a long copy)
 *   (Correctly handle an input file with all constant data)
 *   (Debug the option parsing by re-constructing the command line)
 */
class OptionsCommon: protected OptionsTools
{
public:
};

/**
 * Option processing based on the old ZgyCopy.
 */
class Options: public OptionsCommon
{
public:
  std::string myname;
  int verbose;
  bool nolod;
  bool update;
  bool randomize;
  bool alpha;                   // Still unused.
  bool dumpsqnr;                // Still unused.
  bool native;
#if HAVE_DROPCACHE
  bool dropcache;
#endif
  bool ordered_write;
  std::string sigpipe;
  std::string input;
  std::string output;
  std::string src_sdurl;
  std::string src_sdkey;
  std::string src_token;
  std::string dst_sdurl;
  std::string dst_sdkey;
  std::string dst_token;
  std::array<std::int64_t,4> fakesize;
  std::array<std::int64_t,3> chunksize;
  std::array<std::int64_t,3> obricksize;
  OpenZGY::SampleDataType osamplesize;
  std::array<int,3>       mirror;
  std::array<float,3>     upsample;
  std::array<std::int64_t,3> cropstart;
  std::array<std::int64_t,3> cropsize;
  std::vector<OpenZGY::DecimationType> algorithm;
  finalize_t finalize;
  int lod;                      // Still unused. Maybe not useful.
  int threads;
  int brickcount;
  int sqnr;
  int omp_nest;
  int noisefactor;

public:
  Options(int argc, char **argv);

private:
  static void help(const std::string& myname);
  void show(std::ostream& os) const;
  static const char* short_options();
  static const struct option *long_options();
  bool setoptCommon(int ch, const char *optarg);
  void setopt(int ch, const char *optarg);
  void parse(int argc, char** argv);
  void check();
};

Options::Options(int argc, char **argv)
  : myname(argc >= 1 ? argv[0] : "zgycopyc")
  , verbose(2)
  , nolod(false)
  , update(false)
  , randomize(false)
  , alpha(false)
  , dumpsqnr(false)
  , native(true)
#if HAVE_DROPCACHE
  , dropcache(false)
#endif
  , ordered_write(false)
  , sigpipe()
  , input()
  , output()
  , src_sdurl()
  , src_sdkey()
  , src_token()
  , dst_sdurl()
  , dst_sdkey()
  , dst_token()
  , fakesize(std::array<std::int64_t,4>{0,0,0,0})
  , chunksize(std::array<std::int64_t,3>{64,256,0})
  , obricksize(std::array<std::int64_t,3>{64,64,64})
  , osamplesize(OpenZGY::SampleDataType::unknown)
  , mirror(std::array<int,3>{1,1,1})
  , upsample(std::array<float,3>{1,1,1})
  , cropstart(std::array<std::int64_t,3>{-1,-1,-1})
  , cropsize(std::array<std::int64_t,3>{-1,-1,-1})
  , algorithm()
  , finalize(std::make_pair(OpenZGY::FinalizeAction::BuildDefault, false))
  , lod()
  , threads(1)
  , brickcount(-1)
  , sqnr(0)
  , omp_nest(-1)
  , noisefactor(0)
{
  if (myname.find_last_of("/\\") != std::string::npos)
    myname = myname.substr(myname.find_last_of("/\\")+1);
  if (argc > 1 && argv != nullptr) {
    parse(argc, argv);
    check();
    if (verbose >= 3)
      show(std::cerr);
  }
  else {
    help(myname);
    exit(1);
  }
}

void
Options::help(const std::string& myname)
{
  static std::vector<std::string> help_options {
     "-h, --help:                Output this help text.",
     "-v, --verbose:             Verbose output. May be repeated",
     "-q, --quiet:               No output except errors.",
     "-G, --nolod:               Completely disable lod generation on write.",
     "-u, --update:             *Update (don't truncate) an existing file.",
     "-r, --randomize:           Read and write bricks in random order.",
     //"-a, --alpha:              *Also copy the alpha plane.",
     //"-D, --dumpsqnr:           *Dump table of sqnr vs. compression ratio.",
     "-N, --native:              Read/write as file's type (default).",
     "-F, --float:               Read/write as float.",
#if HAVE_DROPCACHE
     "-U, --dropcache:           Invoke \"dropcache\" before finalize.",
#endif
     "-p, --sigpipe              SIGPIPE disposition.",
     "-i, --input      FILE.zgy: Input file name. If missing, use random data.",
     "-o, --output     FILE.zgy: Output file name. If missing, discard data.",
     "-s, --size       I,J,K,S:  Size, if no input file. E.g. 16x16x16x2.",
     "-b, --bricksize  I,J,K:    Chunk size when copying. E.g. 64x64x64 samples.",
     "-B, --obricksize I,J,K     Brick size in output. E.g. 64x64x64 samples.",
     "-O, --osamplesize type     Output float int16, or int8.",
     "-M, --mirror     I,J,K     Fake a larger survey by mirroring.",
     "-S, --upsample   I,J,K     Fake a larger survey by upsampling with sinc.",
     "-c, --cropstart  I,J,K     Crop from this position. Default center.",
     "-C, --cropsize   I,J,K     Crop to this size.",
     "-g, --algorithm  1,2,N:    LOD algorithms as 3 int: lod1,lod2,lodN.",
     "-f, --finalize   type      full, incremental, keep, etc.",
     "-l, --lod        N:       *Level of detail, 0 = full resolution.",
     "-t, --threads    N:        Number of threads to use for reading.",
     "-T, --uthreads   N:        As -t but writes may be unordered",
     "-n, --brickcount N:        Only copy the first N bricks.",
     "-Q, --sqnr       QUALITY:  Compression quality. Uncompressed if absent.",
     "-z, --omp-nest   N:        Levels of OpenMP nesting.",
     "    --noise      factor    Add noise when copying.",
     "    --src-sdurl  url:      Seismic Store endpoint to copy from.",
     "    --src-sdkey  key:      Seismic Store API key.",
     "    --src-token  token:    Seismic Store authorization.",
     "    --dst-sdurl  url:      Seismic Store endpoint to copy from.",
     "    --dst-sdkey  key:      Seismic Store API key.",
     "    --dst-token  token:    Seismic Store authorization.",
    };
  std::cerr << "Usage: " << myname << " options...\n";
  for (const std::string& s : help_options)
    std::cerr << "    " << s << "\n";
}

void
Options::show(std::ostream& os) const
{
  os << myname << " "
     << (verbose < 2 ? "-" + std::string(2-verbose, 'q') + " ":
         verbose > 2 ? "-" + std::string(verbose-2, 'v') + " ":
         std::string())
     << "--input \"" << input << "\" "
     << "--output \"" << output << "\" "
     << (nolod ? "--nolod " : "")
     << (update ? "--update " : "")
     << (randomize ? "--randomize " : "")
     << (alpha ? "--alpha " : "")
     << (dumpsqnr ? "--dumpsqnr " : "")
     << (native ? "--native " : "--float ")
#if HAVE_DROPCACHE
     << (dropcache ? "--dropcache " : "")
#endif
     << (sigpipe.empty() ? "" : ("--sigpipe " + sigpipe + " "))
     << (noisefactor ? "--noise " + std::to_string(noisefactor) : std::string())
     << (src_sdurl.empty() ? "" : "--src-sdurl " + src_sdurl)
     << (src_sdkey.empty() ? "" : "--src-sdkey " + src_sdkey)
     << (src_token.empty() ? "" : "--src-token " + src_token)
     << (dst_sdurl.empty() ? "" : "--dst-sdurl " + dst_sdurl)
     << (dst_sdkey.empty() ? "" : "--dst-sdkey " + dst_sdkey)
     << (dst_token.empty() ? "" : "--dst-token " + dst_token);

  if (fakesize[0]||fakesize[1]||fakesize[2]||fakesize[3])
    os << "--size "
       << fakesize[0] << "x"
       << fakesize[1] << "x"
       << fakesize[2] << "x"
       << fakesize[3] << " ";

  os << "--bricksize "
     << chunksize[0] << "x"
     << chunksize[1] << "x"
     << chunksize[2] << " ";

  if (obricksize[0] != 64 || obricksize[1] != 64 || obricksize[2] != 64)
    os << "--obricksize "
       << obricksize[0] << "x"
       << obricksize[1] << "x"
       << obricksize[2] << " ";

  if (mirror[0] != 1 || mirror[1] != 1 || mirror[2] != 1)
    os << "--mirror "
       << mirror[0] << "x"
       << mirror[1] << "x"
       << mirror[2] << " ";

  if (upsample[0] != 1 || upsample[1] != 1 || upsample[2] != 1)
    os << "--upsample "
       << upsample[0] << "x"
       << upsample[1] << "x"
       << upsample[2] << " ";

  if (cropstart[0] >= 0 || cropstart[1] >= 0 || cropstart[2] >= 0)
    os << "--cropstart "
       << cropstart[0] << "x"
       << cropstart[1] << "x"
       << cropstart[2] << " ";

  if (cropsize[0] >= 0 || cropsize[1] >= 0 || cropsize[2] >= 0)
    os << "--cropsize "
       << cropsize[0] << "x"
       << cropsize[1] << "x"
       << cropsize[2] << " ";

  switch (osamplesize) {
  default:
  case OpenZGY::SampleDataType::unknown:
    break;
  case OpenZGY::SampleDataType::int8:
    os << "--osamplesize int8 ";
    break;
  case OpenZGY::SampleDataType::int16:
    os << "--osamplesize int16 ";
    break;
  case OpenZGY::SampleDataType::float32:
    os << "--osamplesize float32 ";
    break;
  }

  if (!algorithm.empty()) {
    for (size_t ii=0; ii<algorithm.size(); ++ii)
      if (ii == 0)
        os << "--algorithm " << (int)algorithm[ii];
      else
        os << "," << (int)algorithm[ii];
    os << " ";
  }

  if (finalize != finalize_t(OpenZGY::FinalizeAction::BuildDefault, false))
    os << "--finalize=" << showFinalizeAction(finalize) << " ";

  if (lod != 0)
    os << "--lod=" << lod << " ";
  if (threads != 1)
    os << (ordered_write ? "--threads=" : "--uthreads=") << threads << " ";
  if (brickcount >= 0)
    os << "--brickcount=" << brickcount << " ";
  if (sqnr > 0)
    os << "--sqnr=" << sqnr << " ";
  if (omp_nest != -1)
    os << "--omp-nest=" << omp_nest << " ";
  os  << std::endl;
}

const char*
Options::short_options()
{
  return "hvqGuraDNFUp:i:o:s:l:b:B:O:g:t:T:n:Q:M:S:c:C:";
}

const struct option *
Options::long_options()
{
  static const struct option result[] = {
     {"help",       no_argument,       0,  'h' },
     {"verbose",    no_argument,       0,  'v' },
     {"quiet",      no_argument,       0,  'q' },
     {"nolod",      no_argument,       0,  'G' },
     {"update",     no_argument,       0,  'u' },
     {"randomize",  no_argument,       0,  'r' },
     {"alpha",      no_argument,       0,  'a' },
     {"dumpsqnr",   no_argument,       0,  'D' },
     {"native",     no_argument,       0,  'N' },
     {"float",      no_argument,       0,  'F' },
#if HAVE_DROPCACHE
     {"dropcache",  no_argument,       0,  'U' },
#endif
     {"sigpipe",    required_argument, 0,  'p' },
     {"input",      required_argument, 0,  'i' },
     {"output",     required_argument, 0,  'o' },
     {"size",       required_argument, 0,  's' },
     {"bricksize",  required_argument, 0,  'b' },
     {"obricksize", required_argument, 0,  'B' },
     {"osamplesize",required_argument, 0,  'O' },
     {"mirror",     required_argument, 0,  'M' },
     {"upsample",   required_argument, 0,  'S' },
     {"cropstart",  required_argument, 0,  'c' },
     {"cropsize",   required_argument, 0,  'C' },
     {"algorithm",  required_argument, 0,  'g' },
     {"finalize" ,  required_argument, 0,  'f' },
     {"threads",    required_argument, 0,  't' },
     {"uthreads",   required_argument, 0,  'T' },
     {"brickcount", required_argument, 0,  'n' },
     {"lod",        required_argument, 0,  'l' },
     {"sqnr",       required_argument, 0,  'Q' },
     {"snr",        required_argument, 0,  'Q' },
     {"omp-nest",   required_argument, 0,  'z' },
     {"noise",      required_argument, 0,  '\001' },
     {"src-sdurl",  required_argument, 0,  '\002' },
     {"src-sdkey",  required_argument, 0,  '\003' },
     {"src-token",  required_argument, 0,  '\004' },
     {"dst-sdurl",  required_argument, 0,  '\005' },
     {"dst-sdkey",  required_argument, 0,  '\006' },
     {"dst-token",  required_argument, 0,  '\007' },
     {0,            0,                 0,  0 }
    };
  return result;
}

bool
Options::setoptCommon(int ch, const char *optarg)
{
  switch (ch) {
  case 'v': ++verbose; break;
  case 'q': --verbose; break;
  case 'G': nolod     = true; break;
  case 'u': update    = true; break;
  case 'r': randomize = true; break;
  case 'a': throw std::runtime_error("--alpha not supported"); //alpha     = true; break;
  case 'D': throw std::runtime_error("--dumpsqnr not supported"); //dumpsqnr  = true; break;
  case 'N': native    = true; break;
  case 'F': native    = false; break;
#if HAVE_DROPCACHE
  case 'U': dropcache = true; break;
#endif
  case 'p': sigpipe   = optarg; break;

  case 'i':
    if (fakesize[0]||fakesize[1]||fakesize[2]||fakesize[3])
      throw std::runtime_error("Do not specify both --size and --input");
    input = std::string(optarg ? optarg : "");
    break;

  case 'o':
    output = std::string(optarg ? optarg : "");
    break;

  case 's':
    {
      if (!input.empty())
        throw std::runtime_error("Do not specify both --size and --input");
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 4)
        throw std::runtime_error("command line: size option needs 4 ints");
      fakesize[0] = tmp.at(0);
      fakesize[1] = tmp.at(1);
      fakesize[2] = tmp.at(2);
      fakesize[3] = tmp.at(3);
    }
    break;

  case 'b':
    {
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -b option needs 3 ints");
      chunksize[0] = tmp.at(0);
      chunksize[1] = tmp.at(1);
      chunksize[2] = tmp.at(2);
    }
    break;

  case 'B':
    {
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -B option needs 3 ints");
      obricksize[0] = tmp.at(0);
      obricksize[1] = tmp.at(1);
      obricksize[2] = tmp.at(2);
    }
    break;

  case 'M':
    {
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -M option needs 3 ints");
      if (tmp.at(0) < 1 || tmp.at(1) < 1 || tmp.at(2) < 1)
        throw std::runtime_error("command line: Mirrors must be at least 1");
      mirror[0] = tmp.at(0);
      mirror[1] = tmp.at(1);
      mirror[2] = tmp.at(2);
    }
    break;

  case 'S':
    {
      std::vector<double> tmp = getdoublelist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -S option needs 2 or 3 numbers");
      if (tmp.at(0) < 1 || tmp.at(1) < 1 || tmp.at(2) < 1)
        throw std::runtime_error("command line: Upsampling must be at least 1");
      upsample[0] = (float)tmp.at(0);
      upsample[1] = (float)tmp.at(1);
      upsample[2] = (float)tmp.at(2);
    }
    break;

  case 'c':
    {
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -c option needs 3 ints");
      cropstart[0] = tmp.at(0);
      cropstart[1] = tmp.at(1);
      cropstart[2] = tmp.at(2);
    }
    break;

  case 'C':
    {
      std::vector<int> tmp = getlist(optarg);
      if (tmp.size() != 3)
        throw std::runtime_error("command line: -C option needs 3 ints");
      cropsize[0] = tmp.at(0);
      cropsize[1] = tmp.at(1);
      cropsize[2] = tmp.at(2);
    }
    break;

  case 'O':
    if (0==strcmp(optarg, "float32") || 0==strcmp(optarg, "float"))
      osamplesize = OpenZGY::SampleDataType::float32;
    else if (0==strcmp(optarg, "int16"))
      osamplesize = OpenZGY::SampleDataType::int16;
    else if (0==strcmp(optarg, "int8"))
      osamplesize = OpenZGY::SampleDataType::int8;
    else
      throw std::runtime_error("command line: sample type must be float32, int16, or int8");
    // --osamplesize implies --float, so in fact the --float option is
    // now almost redundant. If you want "don't change value type but
    // do read/write as float" then --float means you don't need to
    // state the existing data type.
    native = false;
    break;

  case 'g': algorithm  = getDecimationTypes(optarg); break;
  case 'f': finalize   = getFinalizeAction(optarg); break;
  case 'l': lod        = geti(optarg); break;
  case 't': threads    = geti(optarg); ordered_write = true; break;
  case 'T': threads    = geti(optarg); ordered_write = false; break;
  case 'n': brickcount = geti(optarg); break;
  case 'Q': sqnr       = geti(optarg); break;
  case 'z': omp_nest   = geti(optarg); break;

  case 1: noisefactor  = geti(optarg); break;
  case 2: src_sdurl    = optarg; break;
  case 3: src_sdkey    = optarg; break;
  case 4: src_token    = optarg; break;
  case 5: dst_sdurl    = optarg; break;
  case 6: dst_sdkey    = optarg; break;
  case 7: dst_token    = optarg; break;
  default: return false;
  }
  return true;
}

void
Options::setopt(int ch, const char *optarg)
{
  if (!setoptCommon(ch, optarg)) {
    switch (ch) {
    case 'h':
      help(myname);
      exit(1);
      break;
    default:
      throw std::runtime_error("command line: unknown option");
      help(myname);
      exit(1);
    }
  }
  if (sqnr > 0 &&
      osamplesize != OpenZGY::SampleDataType::float32 &&
      osamplesize != OpenZGY::SampleDataType::unknown) {
    throw std::runtime_error("command line: Don't use --osamplesize with compressed output files");
  }
  if (update && osamplesize != OpenZGY::SampleDataType::unknown)
    throw std::runtime_error("command line: Don't use --osamplesize with --update");
  if (update && noisefactor)
    throw std::runtime_error("command line: Don't use --noise with --update");
}

#ifdef HAVE_GETOPT

void
Options::parse(int argc, char** argv)
{
  int ch;
  while ((ch = getopt_long(argc, argv, short_options(), long_options(), nullptr)) >= 0) {
    setopt(ch, optarg);
  }
  // Deliberately not allowing "zgycopy infile outfile". It is safer to
  // have the user be explicit in designating which is the output.
  // I.e. use: "zgycopy -i infile -o outfile"
  if (optind != argc)
    throw std::runtime_error("command line: no non-option arguments expected.");
}

#else

void
Options::parse(int argc, char** argv)
{
  // Use a crude replacement for getopt_long(). No combining of options, so
  // e.g. -rq is not allowed, use -r -q. Long options cannot be abbreviated.
  // Yes there are a gazillion proper implementations on the net, but
  // figuring out the license stuff is not worth the trouble.
  for (int ii=1; ii<argc; ++ii) {
    const char *arg = argv[ii];
    bool shortopt = (strlen(arg) == 2 && arg[0] == '-' && arg[1] != '-');
    bool longopt  = (strlen(arg) >= 3 && arg[0] == '-' && arg[1] == '-' && arg[2] != '-');
    const option* opt;
    for (opt = long_options(); opt->name; ++opt) {
      if ((shortopt && arg[1] == opt->val) ||
          (longopt && 0 == strcmp(opt->name, &arg[2])))
        break;
    }
    if (!opt->name) {
      if (arg[0] != '\0' && arg[0] != '-')
        throw std::runtime_error("Not an option: " + std::string(arg));
      else
        throw std::runtime_error("Unknown option: " + std::string(arg));
    }
    switch (opt->has_arg) {
    case no_argument:
      setopt(opt->val, nullptr);
      break;
    case required_argument:
      if (ii+1 >= argc || argv[ii+1][0] == '-')
        throw std::runtime_error("Missing mandatory argument to " +
                                 std::string(arg));
      setopt(opt->val, argv[++ii]);
      break;
    case optional_argument:
    default:
      throw std::runtime_error("This argument type not supported.");
    }
  }
}

#endif

void
Options::check()
{
  if (input.empty() && fakesize[0]*fakesize[1]*fakesize[2]*fakesize[3] == 0)
    throw std::runtime_error("Must specify either --input or --size");
  //if (!output.empty() && threads != 1)
  //  throw std::runtime_error("Multi threading only allowed when discarding data.");
}

/*=========================================================================*/
/*   END OPTION PROCESSING   ==============================================*/
/*=========================================================================*/

/**
 * Read from file and cache the SAuth token to use.
 * Refresh the cache every 5 minutes.
 *
 * This allows the token to be updated while the copy is in progress.
 * Which may be needed if the copy takes more than 60 minutes and
 * the user is using a regular token.
 *
 * The way the user updates the token is a pretty huge kludge:
 *
 * repeat as needed (and at least 5 minutes before token expires)
 *
 *   sdutil auth login
 *   sdutil auth idtoken > /tmp/token$$.new
 *   mv /tmp/tojen$$.new /tmp/token$$
 *
 * zgycopyc --src-token /tmp/token$$ --dst-token /tmp/token$$
 */
class SDTokenUpdater
{
  mutable std::mutex mutex_;
  std::string filename_;
  std::string token_;
  time_t lastrefresh_;
  std::int64_t last_file_check_;
  int verbose_;
public:
  explicit SDTokenUpdater(const std::string& filename, int verbose)
    : mutex_()
    , filename_(filename)
    , token_()
    , lastrefresh_(0)
    , last_file_check_(0)
    , verbose_(verbose)
  {
  }
  SDTokenUpdater(const SDTokenUpdater&) = delete;
  SDTokenUpdater& operator=(const SDTokenUpdater&) = delete;
  ~SDTokenUpdater() {
    if (verbose_ >= 3)
      std::cerr << "Token updater for \"" << filename_
                << "\"is done." << std::endl;
  }

  /**
   * Return the cached token if it is less then one minute since we cached it.
   * Otherwise read the possibly refreshed token from file.
   */
  std::string operator()() {
    std::lock_guard<std::mutex> lk(mutex_);
    waitForTokenInTokenCallback(filename_, &last_file_check_);
    if (::time(nullptr) - lastrefresh_ >= 60) {
      std::ifstream file(filename_);
      if (!file.good())
        throw std::runtime_error("Cannot open token file \"" + filename_ + "\".");
      std::string old = token_;
      std::getline(file, token_);
      if (file.good() && token_.empty())
        std::getline(file, token_);
      if (file.fail() || file.bad())
        throw std::runtime_error("Cannot read token file \"" + filename_ + "\".");
      lastrefresh_ = time(nullptr);
      // Note: No longer skip the verbose logging if the token did not change.
      // Because, the callback now only gets invoked when the token will soon
      // expire. Which means there had better be a new one available.
      if (verbose_ >= 3 || (verbose_ >= 2 && !old.empty())) {
        std::cerr << "\nToken refresh: "
                  << (token_.empty() ? "<empty>" :
                      token_.size() < 20 ? token_ :
                      "..." + token_.substr(token_.size()-5))
                  << std::endl;
      }
    }
    return token_;
  }

  /**
   * Create an instance and wrap it in a lambda in a way that handles
   * reference counting properly. The instance should be destructed
   * when the last functor referencing it goes out of scope.
   * Note that IOContext is cloneable, hence the functor must be
   * copyable. While SDTokenUpdater is noncopyable.
   *
   * Yes this is rather pedantic. If I just leak the SDTokenUpdater
   * instances this would not have been needed. But just in case I
   * need to copy the code somewhere else I'll try to do it right.
   */
  static SeismicStoreIOContext::tokencb_t create(const std::string& filename, int verbose)
  {
    auto updater = std::make_shared<SDTokenUpdater>(filename, verbose);
    auto lambda = [updater]() -> std::string {return (*updater)();};
    return lambda;
  }
};

SeismicStoreIOContext getContext(const Options& opt, bool read, int verbose)
{
  static const auto isfile = [](const std::string& s) {
    return !s.empty() && (s.front() == '/' || s.front() == '\\');
  };
  SeismicStoreIOContext context;
  if (read) {
    context.sdurl   (!opt.src_sdurl.empty() ? opt.src_sdurl :
                     Environment::getStringEnv("OPENZGY_SDURL"));
    context.sdapikey(!opt.src_sdkey.empty() ? opt.src_sdkey :
                     Environment::getStringEnv("OPENZGY_SDAPIKEY"));
    if (isfile(opt.src_token)) {
      context.sdtokencb(SDTokenUpdater::create(opt.src_token, opt.verbose), "");
    }
    else if (!opt.src_token.empty()) {
      context.sdtoken(opt.src_token);
    }
    else {
      context.sdtoken(Environment::getStringEnv("OPENZGY_TOKEN"), "");
    }
  }
  else {
    context.sdurl   (!opt.dst_sdurl.empty() ? opt.dst_sdurl :
                     Environment::getStringEnv("OPENZGY_SDURL"));
    context.sdapikey(!opt.dst_sdkey.empty() ? opt.dst_sdkey :
                     Environment::getStringEnv("OPENZGY_SDAPIKEY"));
    if (isfile(opt.dst_token)) {
      context.sdtokencb(SDTokenUpdater::create(opt.dst_token, opt.verbose), "");
    }
    else if (!opt.dst_token.empty()) {
      context.sdtoken(opt.dst_token);
    }
    else {
      context.sdtoken(Environment::getStringEnv("OPENZGY_TOKEN"), "");
    }
  }
  if (verbose >= 3)
    std::cerr << (read ? "read" : "write") << " context:\n"
              << context.toString();
  return context;
}

static std::int64_t
roundup(std::int64_t in, std::int64_t step)
{
  return ((in + step - 1) / step) * step;
}

static std::array<std::int64_t,3>
roundup(const std::array<std::int64_t,3>& in, const std::array<std::int64_t,3>& step)
{
  return std::array<std::int64_t,3>
    {roundup(in[0], step[0]),
     roundup(in[1], step[1]),
     roundup(in[2], step[2])};
}

static std::array<std::int64_t,3>
minimum(const std::array<std::int64_t,3>& a, const std::array<std::int64_t,3>& b)
{
  return std::array<std::int64_t,3>
    {std::min(a[0],b[0]),
     std::min(a[1],b[1]),
     std::min(a[2],b[2])};
}

static std::array<std::int64_t,3>
sub(const std::array<std::int64_t,3>& a,
    const std::array<std::int64_t,3>& b)
{
  return std::array<std::int64_t,3>
    {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

void randomize(std::vector<std::array<std::int64_t,3>>& vec)
{
  static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(static_cast<unsigned>(seed));
  for (std::size_t ii = vec.size()-1; ii>0; --ii) {
    std::size_t victim = std::uniform_int_distribution<std::size_t>(0, ii)(generator);
    if (victim != ii)
      std::swap(vec[victim], vec[ii]);
  }
}

template<typename T>
void
addNoise(T *data, const std::array<std::int64_t,3>& size3d, int factor, float valuerange)
{
  if (factor <= 0)
    return;
  const float noiselevel = valuerange / factor;
  const float noisefactor = 1.0f / factor;
  static auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  static std::default_random_engine generator(static_cast<unsigned>(seed));
  std::uniform_real_distribution<float> add_distribution(0.0, noiselevel);
  std::uniform_real_distribution<float> mul_distribution(1.0f-noisefactor, 1.0f);
  const T* end = data + (size3d[0] * size3d[1] * size3d[2]);
  for (T* ptr = data; ptr < end; ++ptr) {
    float add = add_distribution(generator);
    float mul = mul_distribution(generator);
    *ptr = static_cast<T>(*ptr * mul + (*ptr<=0 ? add : -add));
  }
}

#ifndef _WIN32

void sigpipe_report(int)
{
  static const char* msg = "WARNING: caught a SIGPIPE\n";
  (void)write(2, msg, strlen(msg));
}

void sigpipe_fatal(int)
{
  static const char* msg = "FATAL: caught a SIGPIPE\n";
  (void)write(2, msg, strlen(msg));
  _exit(3);
}

void signals(const Options& opt)
{
  if (opt.sigpipe.empty()) {
    // SIG_DFL or SIG_IGN inherited from parent
  }
  else if (opt.sigpipe == "inherit") {
    // SIG_DFL or SIG_IGN inherited from parent
  }
  else if (opt.sigpipe == "ignore") {
    ::signal(SIGPIPE, SIG_IGN);
  }
  else if (opt.sigpipe == "default") {
    ::signal(SIGPIPE, SIG_DFL);
  }
  else if (opt.sigpipe == "report") {
    struct sigaction act{{0}};
    act.sa_handler = sigpipe_report;
    ::sigaction(SIGPIPE, &act, NULL);
  }
  else if (opt.sigpipe == "fatal") {
    struct sigaction act{{0}};
    act.sa_handler = sigpipe_fatal;
    ::sigaction(SIGPIPE, &act, NULL);
  }
  else {
    std::cerr << "Allowable values for --sigpipe:\n"
              << "   inherit | ignore | default | report | fatal\n"
              << std::flush;
    throw std::runtime_error("Invalid argument to --sigpipe");
  }
}

#else

void signals(const Options& opt)
{
  // Signals might actually work on Windows but I don't have time
  // to test that now. And yagni unless SIGPIPE problems start
  // popping up on Windows as well.
  if (!opt.sigpipe.empty()) {
    throw std::runtime_error("The --signal option is not supported on Windows");
  }
}

#endif

/**
 * Suggest a data range to specify when the entire file is constant value.
 * Using low==high is not allowed because it triggers several corner cases.
 * The result should be zero-centric and map storage-zero to the single value.
 * Apart from that, just suggest something reasonable.
 *
 * This function might be used when creating a constant-cube from scratch.
 * It is also useful when encountering an old ZGY file from Petrel where
 * low==high was permitted. Adjusting old files could also be done inside
 * the access library. As of this writing it isn't.
 */
void
suggestRange(float value, float dt_lo, float dt_hi, float *lo, float *hi)
{
  if (value == 0) {
    // Choose the range -1..+1, but slightly wider at the low end
    // to make float-0 map exactly to int-0, zero centric.
    *lo = dt_lo / dt_hi;
    *hi = 1;
  }
  else if (value > 0) {
    // Choose the range 0..2*value, but slightly less at the high end
    // to make float-value map exactly to int-0, zero centric.
    *lo = 0;
    *hi = value * (1 - dt_hi / dt_lo);
  }
  else { // value < 0
    // Choose the range 2*value..0, but slightly wider at the low end
    // to make float-value map exactly to int-0, zero centric.
    *lo = value * (1 - dt_lo / dt_hi);
    *hi = 0;
  }
}

/**
 * See the overloaded version for details.
 */
void
suggestRange(float value, OpenZGY::SampleDataType dt, float *lo, float *hi)
{
  switch (dt) {
  case OpenZGY::SampleDataType::int8:
    suggestRange(value, -128, +127, lo, hi);
    break;
  case OpenZGY::SampleDataType::int16:
    suggestRange(value, -32768, +32767, lo, hi);
    break;
  default:
    *lo = *hi = value;
    break;
  }
}

/**
 * \brief Configure OpenMP nested loops.
 *
 * \details
 * Deciding whether to enable nested OpenMP loops is tricky.
 * If using --threads > 1 this causes both read and write calls to the
 * OpenZGY API to be invoked from inside a parallel region. Thus disabling
 * all the OpenMP loops inside OpenZGY by default. Writes are serialized
 * and even ordered when they reach the API. This is handled by locks. But
 * as long as the calls come from inside a parallel region they still
 * disable the loops at the lower level.
 *
 * So: With --threads > 1 the code should always allow at least one level
 * of nested loops. It may or may not be a good idea for writes. But having
 * separate settings for reads and writes doesn't seem feasible.
 *
 * Actually enabling nesting is also quite fun since it depends on the
 * OpenMP version used. I hope I got it right. And even so, the application
 * could technically pick up a newer version at runtime. So testing the
 * _OPENMP symbol isn't actually correct.
 *
 * Oh, and if not set explicitly then OpenMP might pick up an environent
 * variable that changes the behavior. Here is the behavior of this program:
 *
 *    <li>--omp_level=-1 or unset:
 *        One level if --threads=1, two otherwise.
 *    <li>--omp-level=0:
          OpenMP default, possibly affected by envirinment variables.
 *    <li>--omp-level=N, for 0<N<99:
 *        Exactly that many levels unless OpenMP version is too old.
 *    <li>--omp-level=99:
 *        Set as high as possible.
 *
 * TODO-Low: finalize() might need a different setting than read/write.
 */
void
openmp_config(const Options& opt)
{
  int nesting = opt.omp_nest >= 0 ? opt.omp_nest : opt.threads > 1 ? 2 : 1;
  if (nesting != 0) {
#if _OPENMP >= 201611
    {
      // OpenMP 5.0 preview 1 or later.
      // TODO-Test, not tested.
      // omp_set_nested(true); // Deprecated and prints warning at runtime.
      if (nesting < 99)
        omp_set_max_active_levels(nesting);
    }
#elif _OPENMP >= 200805
    {
      // OpenMP 3.x and 4.x. Currently used in most Linux distros.
      omp_set_nested(true); // In this case both must be used.
      if (nesting < 99)
        omp_set_max_active_levels(nesting);
    }
#else
    {
      // OpenMP 2.x or older. Currently the case in Windows.
      omp_set_nested(true); // In this case it is just an on/off switch.
      nesting = 99; // tell user we set it to infinite.
    }
#endif
    if (opt.verbose >= 2) {
      std::cerr << "OpenMP nesting set to "
                << (nesting >= 99 ? std::string("infinite") :
                    std::to_string(nesting))
                << ".\n";
    }
  }
}

void readchunk(
    const std::shared_ptr<OpenZGY::IZgyReader>& r,
    const std::shared_ptr<OpenZGY::IZgyWriter>& w,
    const std::array<std::int64_t,3>& pos,
    const std::array<std::int64_t,3>& chunksize,
    const std::array<std::int64_t,3>& surveysize,
    void *buffer, OpenZGY::SampleDataType dt,
    SummaryPrintingTimerEx& rtimer,
    int noisefactor)
{
  const float range = std::max(0.0f, r->datarange()[1] - r->datarange()[0]);
  const std::array<std::int64_t,3> realsize =
    minimum(chunksize, sub(surveysize, pos));
  if (realsize[0] > 0 && realsize[1] > 0 && realsize[2] > 0) {
    SimpleTimerEx rt(rtimer);
    int bps = 0;
    switch (dt) {
    case OpenZGY::SampleDataType::int8:
      bps=1;
      r->read(pos, realsize, static_cast<std::int8_t*>(buffer), 0);
      addNoise(static_cast<std::int8_t*>(buffer), realsize, noisefactor, range);
      break;
    case OpenZGY::SampleDataType::int16:
      r->read(pos, realsize, static_cast<std::int16_t*>(buffer), 0);
      addNoise(static_cast<std::int16_t*>(buffer), realsize, noisefactor, range);
      bps=2;
      break;
    case OpenZGY::SampleDataType::float32:
    default:
      r->read(pos, realsize, static_cast<float*>(buffer), 0);
      addNoise(static_cast<float*>(buffer), realsize, noisefactor, range);
      bps=4;
      break;
    }
    rtimer.addBytesRead(realsize[0]*realsize[1]*realsize[2]*bps);
  }
}

void writechunk(
    const std::shared_ptr<OpenZGY::IZgyReader>& r,
    const std::shared_ptr<OpenZGY::IZgyWriter>& w,
    const std::array<std::int64_t,3>& pos,
    const std::array<std::int64_t,3>& chunksize,
    const std::array<std::int64_t,3>& surveysize,
    const void *buffer, OpenZGY::SampleDataType dt,
    SummaryPrintingTimerEx& wtimer)
{
  const std::array<std::int64_t,3> realsize =
    minimum(chunksize, sub(surveysize, pos));
  if (realsize[0] > 0 && realsize[1] > 0 && realsize[2] > 0) {
    // OpenZGY will also use an exclusive lock for writes.
    // But by setting our own the "wtimer" will exclude the
    // wait time. Which would have made "wtimer" not very
    // useful in the multithreaded case.
    static std::mutex mutex;
    std::lock_guard<std::mutex> lk(mutex);
    SimpleTimerEx wt(wtimer);
    int bps = 0;
    switch (dt) {
    case OpenZGY::SampleDataType::int8:
      bps=1;
      w->write(pos, realsize, static_cast<const std::int8_t*>(buffer));
      break;
    case OpenZGY::SampleDataType::int16:
      bps=2;
      w->write(pos, realsize, static_cast<const std::int16_t*>(buffer));
      break;
    case OpenZGY::SampleDataType::float32:
    default:
      bps=4;
      w->write(pos, realsize, static_cast<const float*>(buffer));
      break;
    }
    wtimer.addBytesWritten(realsize[0]*realsize[1]*realsize[2]*bps);
  }
}

void
copy(const Options& opt, SummaryPrintingTimerEx& rtimer, SummaryPrintingTimerEx& wtimer, SummaryPrintingTimerEx& rwtimer, SummaryPrintingTimerEx& ftimer, SummaryPrintingTimerEx& stimer)
{
  using namespace OpenZGY;
  FancyProgressWithDots p1(opt.verbose >= 1 ? 51 : 0);
  FancyProgressWithDots p2(opt.verbose >= 1 ? 51 : 0);
  ZgyWriterArgs args;

  SeismicStoreIOContext rcontext(getContext(opt, true, opt.verbose));
  SeismicStoreIOContext wcontext(getContext(opt, false, opt.verbose));
  std::shared_ptr<IZgyReader> r = !opt.input.empty() ?
    IZgyReader::open(opt.input, &rcontext):
    Test::ZgyReaderMock::mock(opt.fakesize);
  if (opt.lod != 0 ||
      opt.cropstart[0] >= 0 || opt.cropstart[1] >= 0 || opt.cropstart[2] >= 0 ||
      opt.cropsize[0]  >= 0 || opt.cropsize[1]  >= 0 || opt.cropsize[2] >= 0)
  {
    r = std::make_shared<ZgyReadLodCrop>(r, opt.lod, opt.cropstart, opt.cropsize);
  }
  if (opt.upsample[0] != 1 || opt.upsample[1] != 1 || opt.upsample[2] != 1)
    //r = std::make_shared<ZgyReaderUpsample>(r, opt.upsample);
    throw std::runtime_error("Command line: --upsample not yet implemented");
  // This is for mirroring; even if there is a crop/lod modifier,
  // that reader is real enough for our purposes.
  std::shared_ptr<IZgyReader> real_r = r;
  if (opt.mirror[0] != 1 || opt.mirror[1] != 1 || opt.mirror[2] != 1)
    r = std::make_shared<ZgyReaderMirror>(r, opt.mirror);
  // TODO-MIRROR: See notes below about the zero chunk size.
  // TODO-MIRROR: Adding upsample might need to be done above
  // the mirror wrapper, to get the "aligned to chunk size" rule right.
  if (!opt.update) {
    args.metafrom(r);
    if (opt.osamplesize != OpenZGY::SampleDataType::unknown)
      args.datatype(opt.osamplesize);
    if (opt.obricksize[0]>0 && opt.obricksize[1]>0 && opt.obricksize[2]>0)
      args.bricksize(opt.obricksize[0], opt.obricksize[1], opt.obricksize[2]);
  }
  args.filename(opt.output);
  if (opt.sqnr > 0)
    args.zfp_compressor(static_cast<float>(opt.sqnr))
        .zfp_lodcompressor(static_cast<float>(opt.sqnr))
        .datatype(OpenZGY::SampleDataType::float32);
  args.iocontext(&wcontext);

  // Existing files with integral storage and all samples set to
  // the same value need special handling to appear consistent.
  // See explanation in the documentation of the format.
  // Note that w->datatype() or args.datatype isn't available.
  const OpenZGY::SampleDataType w_datatype =
    opt.update ? OpenZGY::SampleDataType::unknown :
    opt.osamplesize != OpenZGY::SampleDataType::unknown ? opt.osamplesize :
    r->datatype();
  if ((w_datatype == OpenZGY::SampleDataType::int8 ||
       w_datatype == OpenZGY::SampleDataType::int16) &&
      r->raw_datarange()[0] == r->raw_datarange()[1] &&
      !opt.output.empty()) {
    if (opt.update)
      throw std::runtime_error
        ("Update from an all-constant file is pointless and not implemented");
    const float value = r->raw_datarange()[0];
    float lo=0, hi=0;
    suggestRange(value, w_datatype, &lo, &hi);
    args.datarange(lo, hi);
    if (opt.verbose >= 2) {
      std::stringstream outstring;
      outstring << "Writing a constant-value file"
                << ", range " << lo << " .. " << hi
                << ", value " << value << "\n"
                << "With datarange() it would be "
                << r->datarange()[0] << " .. " << r->datarange()[1] << "\n";
      std::cerr << outstring.str() << std::flush;
    }
    std::shared_ptr<IZgyWriter> w = IZgyWriter::open(args);
    w->writeconst(std::array<int64_t,3>{0,0,0}, w->size(), &value);
    return;
  }

  std::shared_ptr<IZgyWriter> w =
    opt.output.empty() ? Test::ZgyWriterMock::mock(args) :
    opt.update ? IZgyWriter::reopen(args) :
    IZgyWriter::open(args);
  std::shared_ptr<IZgyWriter> real_w = w;

  if (true) {
    // Switch to reading from the real reader and writing to a
    // virtual one that will replicate data to be mirrored.
    // Both r and w will now appear to have the original size.
    // CAVEAT: ZgyWriterMirror is incomplete. It is supposed to
    // wrap an inflated file such that it appears to be the
    // small file we started with. Most of the metadata has not
    // been adjusted. The exception being size(). The bricksize()
    // and datatype() are also safe to call, since they don't change.
    //
    // finalize() should still be called on the underlying file.
    // Technically most of the lowres handling could also have
    // been optimized to only be computed once. But that would
    // become insanely complicated.
    //
    // If this code is not enabled, both r and w will have the
    // inflated size which means the copy loop will process
    // more data. And the input data will be read more than once.
    if (r != real_r) {
      w = std::make_shared<ZgyWriterMirror>(w, opt.mirror);
      r = real_r;
    }
  }

  // I normally want to read full bricks as stored on file. So, while
  // I need to ensure that I don't read too far past the survey's edge
  // I can allow reading up to the next boundary. Reader and Writer
  // might in the future disagree on size and bricksize.
  const std::array<std::int64_t,3> surveysize =
    minimum(roundup(r->size(), r->bricksize()),
            roundup(w->size(), w->bricksize()));

  // chunksize can be passed as zero, meaning as much as possible.
  // Round up to the brick size as this could be more efficient.
  // If mirroring, the count should be the size from *real source*.
  // TODO-MIRROR: If changing to use separate read and write wrappers,
  // this might need to change. Ditto for upsampling.
  const std::array<std::int64_t,3> bs = std::array<std::int64_t,3>
    {
     opt.chunksize[0] ? opt.chunksize[0] : surveysize[0] / opt.mirror[0],
     opt.chunksize[1] ? opt.chunksize[1] : surveysize[1] / opt.mirror[1],
     opt.chunksize[2] ? opt.chunksize[2] : surveysize[2] / opt.mirror[2]
    };

  std::vector<std::array<std::int64_t,3>> tasklist;
  {
    const std::array<std::int64_t,3> size = r->size();
    std::array<std::int64_t,3> pos;
    for (pos[0] = 0; pos[0] < size[0]; pos[0] += bs[0])
      for (pos[1] = 0; pos[1] < size[1]; pos[1] += bs[1])
        for (pos[2] = 0; pos[2] < size[2]; pos[2] += bs[2])
          tasklist.push_back(pos);
  }

  if (opt.update) {
    std::vector<std::array<std::int64_t,3>> newlist;
    for (const std::array<std::int64_t,3>& pos : tasklist) {
      std::pair<bool,double> r_isconst = r->readconst(pos, bs, 0, false);
      std::pair<bool,double> w_isconst = w->readconst(pos, bs, false);
      // Subtle issue: Technically this might turn a constant-0 brick in
      // the source into an unwritten brick in the target. Which creates
      // all kinds of risk because the replacement for a missing brick
      // is chosen when the file is read, not when it us written.
      // But the targes ws probably free of missing bricks already.
      // Not-so-subtle issue: Chunk size must not change, and I cannot
      // detect whether it has or not.
      if (r_isconst.first && w_isconst.first && r_isconst.second == w_isconst.second) {
        // Already has the correct value.
      }
      else if (!w_isconst.first) {
        // If any part of the chunk has been written with non-const data,
        // assume that the entire chunk has already been copied and might
        // be skipped. This is why it is important to (a) use the same
        // source file ehan updating, and (b) use the same chunk size.
      }
      else {
        newlist.push_back(pos);
      }
    }
    if (opt.verbose >= 2)
      std::cerr << tasklist.size() - newlist.size()
                << " of " << tasklist.size()
                << " are empty or have already been copied."
                << std::endl;
    tasklist = newlist;
  }

  if (opt.randomize)
    randomize(tasklist);

  // How many bricks to copy might be changed by the user.
  const std::int64_t total = opt.brickcount >= 0 ?
    std::min((std::int64_t)tasklist.size(), (std::int64_t)opt.brickcount) :
    (std::int64_t)tasklist.size();
  std::atomic<std::int64_t> done(0);
  if (!done.is_lock_free())
    throw std::runtime_error("Consider using plain int for atomics");

  if (opt.verbose >= 2) {
    static auto ssize = [](SampleDataType dt)
                        {
                          return static_cast<std::int64_t>
                            (dt == SampleDataType::int8 ? 1 :
                             dt == SampleDataType::int16? 2 :
                             dt == SampleDataType::float32 ? 4:
                             0);
                        };
    static auto tsize = [](const std::array<std::int64_t,3>& a)
                        {
                          return (a[0]*a[1]*a[2]) / (1024.0*1024.0);
                        };
    std::stringstream outstring;
    outstring << "Copying "
              << (opt.input.empty() ? "Synthetic data" : opt.input)
              << " size "
              << r->size()[0] << "x"
              << r->size()[1] << "x"
              << r->size()[2] << "x"
              << ssize(r->datatype())
              << "\n"
              << "to "
              << (opt.output.empty() ? "Discard data" : opt.output)
              << " size "
              << w->size()[0] << "x"
              << w->size()[1] << "x"
              << w->size()[2] << "x"
              << ssize(w->datatype())
              << "\n"
              << "Chunk size "
              << bs[0] << "x"
              << bs[1] << "x"
              << bs[2]
              << " (" << tsize(bs) * ssize(w->datatype()) << " MB)"
              //<< " LOD algorithms "
              << "\n"
              << "Will copy " << total << " of " << tasklist.size() << " chunks"
              << " (" << tsize(r->size()) * ssize(r->datatype()) << " MB)"
              << (opt.randomize ? " (in random order)" : "")
              << "\n";
    if (opt.noisefactor)
      outstring << "Add 1/" << opt.noisefactor << " noise\n";
    std::cerr << outstring.str() << std::flush;
  }

  // Size of the thread-local copy buffer, in bytes.
  // Use malloc() instead of new[] because the return from malloc is
  // guaranteed safe to cast to any other type. Not sure about new[].
  const std::size_t bufbytes = std::max(sizeof(float),sizeof(std::int16_t))*bs[0]*bs[1]*bs[2];

  // Data type to be used when copying. "native" is more efficient.
  if (opt.native && r->datatype() != w->datatype())
    throw std::runtime_error("Must read/write as float if datatypes don't match.");
  const SampleDataType dt = opt.native ? r->datatype() : SampleDataType::float32;
  MTGuardWithProgress guard(std::ref(p1), total);
  SimpleTimerEx rwt(rwtimer);
#pragma omp parallel num_threads(opt.threads) if (opt.threads > 1)
  {
    if (opt.verbose >= 4 && omp_get_thread_num() == 0) {
      std::stringstream outstring;
      outstring << ("Copy " +
                    std::to_string(total) + " chunks in " +
                    std::to_string(omp_get_num_threads()) +
                    (omp_get_num_threads() == 1 ? " thread.\n" : " threads.\n"));
      std::cerr << outstring.str() << std::flush;
    }
    std::shared_ptr<void> buf(malloc(bufbytes), [](void *d){::free(d);});
    std::int64_t last_file_check_src{0};
    std::int64_t last_file_check_dst{0};

    if (opt.ordered_write) {
#pragma omp for ordered schedule(dynamic,1)
      for (std::int64_t task = 0; task < total; ++task) {
        guard.run([&]()
                {
                  waitForTokenInCopyLoop(opt.src_token, &last_file_check_src);
                  readchunk(r, w, tasklist[task], bs, surveysize,
                            buf.get(), dt, rtimer, opt.noisefactor);
                });
#pragma omp ordered
        guard.run([&]()
                {
                  waitForTokenInCopyLoop(opt.dst_token, &last_file_check_dst);
                  writechunk(r, w, tasklist[task], bs, surveysize,
                            buf.get(), dt, wtimer);
                });
        guard.progress();
      }
    }
    else {
#pragma omp for schedule(dynamic,1)
      for (std::int64_t task = 0; task < total; ++task) {
        guard.run([&]()
                {
                  waitForTokenInCopyLoop(opt.src_token, &last_file_check_src);
                  readchunk(r, w, tasklist[task], bs, surveysize,
                            buf.get(), dt, rtimer, opt.noisefactor);
                });
        guard.run([&]()
                {
                  waitForTokenInCopyLoop(opt.dst_token, &last_file_check_dst);
                  writechunk(r, w, tasklist[task], bs, surveysize,
                            buf.get(), dt, wtimer);
                });
        guard.progress();
      }
    }
  }
  rwt.stop();
  guard.finished();

  // In case finalize or close hangs or crashes, output the timing results
  // now instead of having them printed automatically in main(). If the code
  // threw an exception we won't get here. But in that case the normal
  // mechanism of printing from the timer's destructor is invoked.
  // PS, don't use this code as a quick introduction to how simple my
  // Timer class is. You might get the wrong impression.

  rtimer.print();
  wtimer.print();
  rwt.done(); // Pretend it went out of scope.
  rwtimer.print(); // Ditto for the printing timer that owns it.

  // In case Timer logging is enabled in OpenZGY then get the
  // read-related timers to output now. close the Some timers are
  // printed on close() and some when the writer goes out of scope but
  // that is just an implementation detail. To keep them together I
  // will destruct the writer immediately afyter closing it.
  r->close();
  r.reset();

  if (opt.nolod) {
    SimpleTimerEx ft(ftimer);
    w->close_incomplete();
  }
  else {
#if HAVE_DROPCACHE
    // dropcache simulates having files so huge that the Linux
    // buffer cache won't keep the entire output file in memory.
    // so there will be a significant I/O cost reading back LOD0.
    if (opt.dropcache)
      SYSTEM("/usr/local/bin/dropcache");
#endif
    // Don't report timing for finalizing a mocked output file.
    // Yes it does actually have a (tiny) cost but the user won't
    // expect to see finalize reported at all when discarding the output.
    SimpleTimerEx ft(ftimer);
    real_w->finalize(opt.algorithm, std::ref(p2),
                opt.finalize.first, opt.finalize.second);
    w->close();
    // If Timer logging is also enabled inside OpenZGY there will now
    // be multiple lines of output when the output file is finalized
    // and closed. Unless we are writing to the mocked dummy output
    // file. All these times are related to writing. Even the
    // File::read() timer which now reports reads done while making
    // LOD bricks. Some timers are printed on close() and some when
    // the writer goes out of scope but that is just an implementation
    // detail. To keep them together I will destruct the writer
    // immediately afyter closing it.
    w.reset();
  }
#ifndef _WIN32
  // This is only for performance measurements and I am only enabling
  // it for Linux as I don't know the Windows equivalent. It only
  // makes sense when something was written to a local disk. If there
  // are any dirty buffers after a read+discard or after a write to
  // cloud then they don't belong to us. NOTE: Technically this ought
  // to have been a fsync(fd) to only flush our own data. But that
  // means exposing fsync from FileADT and that isn't worth the
  // trouble. A machine running benchmarks shouldn't have other
  // activity anyway.
  if (!opt.output.empty() && opt.output.substr(0,5) != "sd://") {
    SimpleTimerEx stim(stimer);
    sync();
  }
#endif
  ftimer.print();
  stimer.print();
}

int main(int argc, char **argv)
{
  int verbose = 2;
  try {
    Options options(argc, argv);
    verbose = options.verbose;
#if HAVE_DROPCACHE
    if (options.dropcache)
      SYSTEM("/usr/local/bin/dropcache");
#endif
    signals(options);
    openmp_config(options);
    SummaryPrintingTimerEx stimer("Tool.sync");
    SummaryPrintingTimerEx ftimer("Tool.finalize");
    SummaryPrintingTimerEx rwtimer("Tool.read+write");
    SummaryPrintingTimerEx wtimer("Tool.write");
    SummaryPrintingTimerEx rtimer("Tool.read");
    SummaryPrintingTimerEx ttimer("Tool.TOTAL");
    SimpleTimerEx t1(ttimer);
    copy(options, rtimer, wtimer, rwtimer, ftimer, stimer);
  }
  catch (const std::exception& ex) {
    std::string myname(argc >= 1 ? argv[0] : "zgycopyc");
    if (myname.find_last_of("/\\") != std::string::npos)
      myname = myname.substr(myname.find_last_of("/\\")+1);
    std::cerr << myname << ": " << ex.what() << std::endl;
    exit(1);
  }
  std::cerr << std::flush;
  std::cout << std::flush;
  if (verbose >= 2)
    std::cerr << "Ok, normal termination of zgycopy." << std::endl;
  exit(0);
}
