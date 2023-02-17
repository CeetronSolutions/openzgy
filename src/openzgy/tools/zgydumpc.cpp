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
 * \file zgydumpc.cpp
 * \brief Show metadata for a ZGY file.
 */

#include "../api.h"
#include "../impl/environment.h"
#include "../iocontext.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <chrono>
#include <random>
#include <malloc.h>
#include <cstring>

#ifndef _WIN32
#include <unistd.h>
#include <getopt.h>
#define HAVE_GETOPT
#else
#include <io.h>
#undef HAVE_GETOPT
#endif

using InternalZGY::Environment;
using OpenZGY::SeismicStoreIOContext;

/**
 * Convenience to hard code credentials for testing. Returns an IOContext.
 * Picking up sdurl/sdapikey from the environment is redundant since
 * the library already does this as a fallback.
 */
SeismicStoreIOContext getContext()
{
    return SeismicStoreIOContext()
        .sdurl(Environment::getStringEnv("OPENZGY_SDURL"))
        .sdapikey(Environment::getStringEnv("OPENZGY_SDAPIKEY"))
        .sdtoken(Environment::getStringEnv("OPENZGY_TOKEN") != "" ?
                 Environment::getStringEnv("OPENZGY_TOKEN") :
                 "FILE:carbon.slbapp.com", "");
}

/*=========================================================================*/
/*   OPTION PROCESSING   ==================================================*/
/*=========================================================================*/

/**
 * Option processing based on the old ZgyCopy.
 */
class Options
{
public:
  std::string myname;
  int verbose;
  bool histogram;
  bool offsets;
  bool sorted_offsets;
  bool only_lod0_info;
  std::vector<std::string> inputs;
  Options(int argc, char **argv)
    : myname(argc >= 1 ? argv[0] : "zgydumpc")
    , verbose(1)
    , histogram(false)
    , offsets(false)
    , sorted_offsets(false)
    , only_lod0_info(false)
  {
    if (myname.find_last_of("/\\") != std::string::npos)
      myname = myname.substr(myname.find_last_of("/\\")+1);
    if (argc > 1 && argv != nullptr) {
      parse(argc, argv);
      check();
      if (verbose > 1)
        show(std::cerr);
    }
    else {
      help(myname);
      exit(1);
    }
  }

  static void help(const std::string& myname)
  {
    static std::vector<std::string> help_options
      {
       "-h, --help:           Output this help text.",
       "-v, --verbose:        Verbose output. May be repeated.",
       "-q, --quiet:          Less verbose output. May be repeated.",
       "-H, --histogram:      Show the 256-bin histogram.",
       "-O, --offsets:        Show offset of each brick.",
       "-S, --sorted_offsets: Sort by file offset.",
       "-0, --only-lod0-info: Ignore low resolution.",
      };
    std::cerr << "Usage: " << myname << " options...\n";
    for (const std::string& s : help_options)
      std::cerr << "    " << s << "\n";
  }

  void show(std::ostream& os) const
  {
    os << myname << " "
       << (verbose<=0 ? std::string(" -q") :
           verbose==1 ? std::string():
           " -" + std::string(verbose-1, 'v'))
       << (histogram      ? " --histogram" : "")
       << (offsets        ? " --offsets" : "")
       << (sorted_offsets ? " --sorted-offsets" : "")
       << (only_lod0_info ? " --only-lod0-info" : "");
    for (const std::string& s : inputs)
      os << " '" << s << "'";
    os << std::endl;
  }

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

  static const char* short_options()
  {
    return "hvqHOS0";
  }

  static const struct option *long_options()
  {
    static const struct option result[] = {
       {"help",           no_argument, 0,  'h' },
       {"verbose",        no_argument, 0,  'v' },
       {"quiet",          no_argument, 0,  'q' },
       {"histogram",      no_argument, 0,  'H' },
       {"offsets",        no_argument, 0,  'O' },
       {"sorted_offsets", no_argument, 0,  'S' },
       {"only-lod0-info", no_argument, 0,  '0' },
       {0,                0,           0,  0 }
    };
    return result;
  }

  void setopt(int ch, const char *optarg)
  {
    switch (ch) {
    case 'h':
      help(myname);
      exit(1);
      break;

    case 'v': ++verbose; break;
    case 'q': --verbose; break;
    case 'H': histogram      = true; break;
    case 'O': offsets        = true; break;
    case 'S': sorted_offsets = true; break;
    case '0': only_lod0_info = true; break;

    default:
      help(myname);
      throw std::runtime_error("command line: unknown option. Try --help.");
    }
  }

#ifdef HAVE_GETOPT

  void parse(int argc, char** argv)
  {
    int ch;
    while ((ch = getopt_long(argc, argv, short_options(), long_options(), nullptr)) >= 0) {
      setopt(ch, optarg);
    }
    while (optind < argc)
      inputs.push_back(std::string(argv[optind++]));
  }

#else

  void parse(int argc, char** argv)
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
        if (arg[0] != '\0' && arg[0] != '-') {
          inputs.push_back(std::string(argv[ii]));
          continue;
        }
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

  void check()
  {
    if (inputs.empty()) {
      help(myname);
      throw std::runtime_error("No inputs provided.");
    }
  }
};

/*=========================================================================*/
/*   END OPTION PROCESSING   ==============================================*/
/*=========================================================================*/

void
dump_basic(std::shared_ptr<OpenZGY::IZgyReader> r, const std::string& filename, std::ostream& os)
{
  using namespace OpenZGY;
  using namespace OpenZGY::Formatters;

  auto oldprec = os.precision();
  auto oldflags = os.flags();

  auto ocp = [r](int ix) {
               // Original formats %5d %9g %11.2f
               std::stringstream ss;
               ss << "["
                  << std::setw(5) << r->indexcorners()[ix][0] << ", "
                  << std::setw(5) << r->indexcorners()[ix][1] << "] {"
                  << std::setw(9) << r->annotcorners()[ix][0] << ", "
                  << std::setw(9) << r->annotcorners()[ix][1] << "} ("
                  << std::fixed << std::setprecision(2)
                  << std::setw(11) << r->corners()[ix][0] << ", "
                  << std::setw(11) << r->corners()[ix][1] << ")";
               return ss.str();
             };

  const SampleStatistics stat = r->statistics();
  const SampleHistogram  hist = r->histogram();
  std::shared_ptr<const FileStatistics> filestats = r->filestats();

  // Note that file size in bytes, ZGY version, and projection system
  // are not available in the API. They might not be useful anyway.

  os << "File name                      = '" << filename << "'\n"
     << "File size (bytes)              = " << filestats->fileSize() << "\n"
     << "File format and version        = " << r->datatype() << " ZGY version " << filestats->fileVersion() << "\n"
     //<< "Data identifier                = " << r->dataid() << "\n"
     << "Current data Version           = " << r->verid() << "\n"
     //<< "Previous data version          = " << r->previd() << "\n"
     << "Brick size I,J,K               = " << "(" << r->bricksize()[0] << ", " << r->bricksize()[1] << ", " << r->bricksize()[2] << ")\n"
     << "Number of bricks I,J,K         = " << "(" << r->brickcount()[0][0] << ", " << r->brickcount()[0][1] << ", " << r->brickcount()[0][2] << ")\n"
     << "Number of LODs                 = " << r->nlods() << "\n"
     << "Coding range min/max           = " << std::setprecision(6) << r->datarange()[0] << " " << r->datarange()[1] << " (raw: " << r->raw_datarange()[0] << " " << r->raw_datarange()[1] << ") " << r->size()[0] * r->size()[1] * r->size()[2] << "\n"
     << "Statistical min/max/count      = " << std::setprecision(6) << stat.min << " " << stat.max << " " << stat.cnt << "\n"
     << "Histogram range min/max/count  = " << std::setprecision(6) << hist.minvalue << " " << hist.maxvalue << " " << hist.samplecount << "\n"
     << "Inline start/increment/count   = " << r->annotstart()[0] << " " << r->annotinc()[0] << " " << r->size()[0] << "\n"
     << "Xline  start/increment/count   = " << r->annotstart()[1] << " " << r->annotinc()[1] << " " << r->size()[1] << "\n"
     << "Sample start/increment/count   = " << r->zstart() << " " << r->zinc() << " " << r->size()[2] << "\n"
     << "Horizontal projection system   = " << "?\n" // {r._accessor._metadata._ih._hprjsys}
     << "Horizontal dim/factor/name     = " << r->hunitdim() << " " << r->hunitfactor() << " '" << r->hunitname() << "'\n"
     << "Vertical dim/factor/name       = " << r->zunitdim() << " " << r->zunitfactor() << " '" << r->zunitname() << "'\n"
     << "Ordered Corner Points Legend   = " << "[  <i>,   <j>] { <inline>,   <xline>} (  <easting>,  <northing>)\n"
     << "Ordered Corner Point 1         = " << ocp(0) << "\n"
     << "Ordered Corner Point 2         = " << ocp(1) << "\n"
     << "Ordered Corner Point 3         = " << ocp(2) << "\n"
     << "Ordered Corner Point 4         = " << ocp(3) << "\n"
     << std::flush;

  os.flags(oldflags);
  os.precision(oldprec);
}

void
dump_histogram(std::shared_ptr<OpenZGY::IZgyReader> r, std::ostream& os)
{
  int ii = 0;
  for (std::int64_t c : r->histogram().bins)
    os << "Histogram bin "
       << std::setw(3) << ii++
       << "              = "
       << std::setw(11) << c << "\n";
}

void
dump_summary_brick_offsets(std::shared_ptr<OpenZGY::IZgyReader> r, std::ostream& os)
{
  // TBD -- need access to internals
  // See all_brick(), all_alpha(), summary_brick_offsets() in the Python version.
  r->filestats()->dump(os);
}

void
dump_summary_normal_size(std::shared_ptr<OpenZGY::IZgyReader> r, bool header, std::ostream& os)
{
  std::int64_t normal = 0;  // TBD -- need access to internals
  std::int64_t bytespersample = 0;
  switch (r->datatype()) {
  case OpenZGY::SampleDataType::int8:    bytespersample = 1; break;
  case OpenZGY::SampleDataType::int16:   bytespersample = 2; break;
  case OpenZGY::SampleDataType::float32: bytespersample = 4; break;
  case OpenZGY::SampleDataType::unknown: bytespersample = 0; break;
  }
  std::int64_t bytesperbrick =
    bytespersample *
    r->bricksize()[0] *
    r->bricksize()[1] *
    r->bricksize()[2];
  std::int64_t colsize =
    (r->size()[2] + r->bricksize()[2] - 1) / r->bricksize()[2];

  if (header) {
    os << "Normal LOD0 bricks & col size  = LOD0: "
       << normal << " "
       << (normal * bytesperbrick) / (1024*1024) << " MB "
       << "column "
       << colsize << " "
       << (colsize * bytesperbrick) / (1024*1024) << " MB "
       << "brick " << bytesperbrick << "\n";
  }
  else {
    os << normal << " "
       << (normal * bytesperbrick) / (1024*1024) << " "
       << colsize << " "
       << (colsize * bytesperbrick) / (1024*1024) << " "
       << bytesperbrick << "\n";
  }
}

void
dump_combined_offsets(std::shared_ptr<OpenZGY::IZgyReader> r, std::ostream& os)
{
  // TBD -- need access to internals
  os << "BRICK and ALPHA offsets sorted by address: (NOT IMPLEMENTED)\n";
}

void
dump_brick_offsets(std::shared_ptr<OpenZGY::IZgyReader> r, std::ostream& os)
{
  // TBD -- need access to internals
  os << "BRICK offsets: (NOT IMPLEMENTED)\n";
}

void
dump_alpha_offsets(std::shared_ptr<OpenZGY::IZgyReader> r, std::ostream& os)
{
  // TBD -- need access to internals
  os << "ALPHA offsets: (NOT IMPLEMENTED)\n";
}

void
run(const std::string filename, const Options& opt, std::ostream& os)
{
  SeismicStoreIOContext context(getContext());
  std::shared_ptr<OpenZGY::IZgyReader> r = OpenZGY::IZgyReader::open(filename, &context);
  if (opt.only_lod0_info) {
    dump_summary_normal_size(r, false, os);
    return;
  }
  dump_basic(r, filename, os);
  dump_summary_brick_offsets(r, os);
  dump_summary_normal_size(r, true, os);
  if (opt.histogram) {
    dump_histogram(r, os);
  }
  if (opt.sorted_offsets) {
    dump_combined_offsets(r, os);
  }
  else if (opt.offsets) {
    dump_brick_offsets(r, os);
    dump_alpha_offsets(r, os);
  }
  r->close();
}

// Features from the old C++ zgydump
//    -b --brief BriefInfo()
//    (default)  FullInfo()
//       In new code, --histogram must be requested explicitly.
//       There is no --brief since the two outputs then end up
//       rather similar.
//    -p --performance PerfInfo()
//       normal/empty/const counts now printed unconditionally.
//       Can also consider merging test_zgydump+ and test_isoptimal.
//    -o --offset Offsets()
//       Better than the original, handles const(value) and compressed(size).
//    -s --slice  Slice()
//       Covered by test_show. Might as well keep this as a separate exe.
//    -a --alpha Alpha()
//       Not implemented, as the alpha tiles are deprecated.

int main(int argc, char **argv)
{
  using namespace OpenZGY;
  try {
    Options options(argc, argv);
    for (const std::string& filename : options.inputs) {
      run(filename, options, std::cout);
    }
  }
  catch (const std::exception& ex) {
    std::string myname(argc >= 1 ? argv[0] : "zgydumpc");
    if (myname.find_last_of("/\\") != std::string::npos)
      myname = myname.substr(myname.find_last_of("/\\")+1);
    std::cerr << myname << ": " << ex.what() << std::endl;
    exit(1);
  }
}
