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


#include "../api.h"
#include "../iocontext.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <stdlib.h>
#include <atomic>
#include <cstring>
#include <cmath>

typedef std::array<std::array<double,2>,4> corners_t;
int verbose{2};

std::string
format_corners(const corners_t& index, const corners_t& annot, const corners_t& world)
{
  std::stringstream ss;
  for (int ix=0; ix<4; ++ix) {
    ss << "Corner Point " << ix << " = ["
       << std::fixed << std::setprecision(0)
       << std::setw(5) << index[ix][0] << ", "
       << std::setw(5) << index[ix][1] << "] {"
       << std::fixed << std::setprecision(1)
       << std::setw(8) << annot[ix][0] << ", "
       << std::setw(8) << annot[ix][1] << "} ("
       << std::fixed << std::setprecision(2)
       << std::setw(11) << world[ix][0] << ", "
       << std::setw(11) << world[ix][1] << ")\n";
  }
  return ss.str();
}

std::string
format_corners(const OpenZGY::IZgyReader& r)
{
  return format_corners(r.indexcorners(), r.annotcorners(), r.corners());
}

/**
 * Return largest difference between any of the 4 corners.
 * The corners might represent index, annot, or world values.
 */
double
abs_error(const corners_t& c1, const corners_t& c2)
{
  double max_error{0};
  for (int ii=0; ii<4; ++ii) {
    double err = std::hypot(c1[ii][0] - c2[ii][0], c1[ii][1] - c2[ii][1]);
    max_error = std::max(max_error, err);
  }
  return max_error;
}

void
check_single_value(double measured, double allowed, const std::string& what)
{
  if (allowed >= 0 && measured > allowed)
    throw std::runtime_error(what + " " + std::to_string(measured) +
                             " exceeds allowable " + std::to_string(allowed));
}

void
check_size_error(const OpenZGY::IZgyReader& r1, const OpenZGY::IZgyReader& r2)
{
  std::stringstream ss;
  const std::array<std::int64_t,3> a(r1.size());
  const std::array<std::int64_t,3> b(r2.size());
  ss << "Size A: (" << a[0] << ", " << a[1] << ", " << a[2]
     << ")    B: (" << b[0] << ", " << b[1] << ", " << b[2] << ")";
  if (verbose >= 2) {
    std::cout << ss.str() << std::endl;
  }
  if (a[0] != b[0] || a[1] != b[1] || a[2] != b[2])
    throw std::runtime_error("Survey size mismatch " + ss.str());
}

/**
 * Print and check for errors in annotation- and world corners.
 * And index corners for completeness although those are computed
 * from size and are very unlikely to be wrong.
 */
void
check_corners_error(
     const OpenZGY::IZgyReader& r1, const OpenZGY::IZgyReader& r2,
     double epsilon_index, double epsilon_annot, double epsilon_world)
{
  if (verbose >= 2) {
    std::cout << "A Lattice\n" << format_corners(r1) << std::flush;
    std::cout << "B Lattice\n" << format_corners(r2) << std::flush;
  }
  double abs_error_index = abs_error(r1.indexcorners(), r2.indexcorners());
  double abs_error_annot = abs_error(r1.annotcorners(), r2.annotcorners());
  double abs_error_world = abs_error(r1.corners(), r2.corners());
  if (verbose >= 2) {
    std::cout << "Max error:"
              << " index " << abs_error_index
              << " annot " << abs_error_annot
              << " world " << abs_error_world
              << "\n" << std::flush;
  }
  check_single_value(abs_error_index, epsilon_index, "absolute index error");
  check_single_value(abs_error_annot, epsilon_annot, "absolute annot error");
  check_single_value(abs_error_world, epsilon_world, "absolute world error");
}

/**
 * Read both files and compare sample values.
 */
void
check_sample_error(
     const OpenZGY::IZgyReader& r1, const OpenZGY::IZgyReader& r2,
     double epsilon_value)
{
  double sum{0}, ssq{0}, err{0}, range{0};
  OpenZGY::ProgressWithDots p1(verbose >= 1 ? 51 : 0);
  const std::array<std::int64_t,3> size = r1.size();
  const std::array<std::int64_t,3> column{64, 64, size[2]};
  const std::int64_t samples = column[0] * column[1] * column[2];
  std::int64_t done{0};
  const std::int64_t total{((size[0]+63)/64) * ((size[1]+63)/64)};
  std::shared_ptr<float> buf1(new float[samples]);
  std::shared_ptr<float> buf2(new float[samples]);
  for (std::int64_t ii=0; ii < size[0]; ii += column[0]) {
    for (std::int64_t jj=0; jj < size[1]; jj += column[1]) {
      const std::int64_t ni = std::min(column[0], size[0] - ii);
      const std::int64_t nj = std::min(column[1], size[1] - jj);
      r1.read(std::array<std::int64_t,3>{ii, jj, 0},
              std::array<std::int64_t,3>{ni, nj, column[2]},
              buf1.get());
      r2.read(std::array<std::int64_t,3>{ii, jj, 0},
              std::array<std::int64_t,3>{ni, nj, column[2]},
              buf2.get());
      const float *ptr1 = buf1.get();
      const float *ptr2 = buf2.get();
      const float *const ptr1end = ptr1 + (ni * nj * column[2]);
      for (; ptr1 < ptr1end; ++ptr1, ++ptr2) {
        const double diff = std::abs(*ptr1 - *ptr2);
        range = std::max(range, (double)std::abs(*ptr1));
        range = std::max(range, (double)std::abs(*ptr2));
        err   = std::max(err, diff);
        sum  += diff;
        ssq  += (diff * diff);
        static int dejavu{0};
        if (verbose >= 3 && diff > 20 && ++dejavu < 10)
          std::cout << "Error:"
                    << " a: " << *ptr1
                    << " b: " << *ptr2
                    << " at (" << ii << "," << jj << ")"
                    << " offset " << (ptr1 - buf1.get()) << "\n";
      }
      if (ii==0 && jj==0 && verbose >= 3) {
        std::cout << "1st: "
                  << buf1.get()[0] << ", "
                  << buf1.get()[1] << ", "
                  << buf1.get()[2] << ", "
                  << buf1.get()[3] << "\n"
                  << "2nd: "
                  << buf2.get()[0] << ", "
                  << buf2.get()[1] << ", "
                  << buf2.get()[2] << ", "
                  << buf2.get()[3] << "\n";
      }
      p1(++done, total);
    }
  }
  std::int64_t total_samples = size[0] * size[1] * size[2];
  if (verbose >= 2) {
    std::cout << "differences in value:"
              << " max " << err
              << " avg " << (sum / total_samples)
              << " rms " << std::sqrt(ssq / total_samples)
              << " range " << range
              << " relative-error " << err / range
              << std::endl;
  }
  check_single_value(err/range, epsilon_value, "relative sample error");
}

/**
 * Rudimentary compare of two ZGY files, read as float regardless
 * of the actual value type. Computes max absolute error, average
 * error, and root-mean-square of errors.
 *
 * Raises an exception if the two files differ by more than the
 * specified allowable amount.
 *
 *   - epsilon_value: Max relative error in sample values.
 *   - epsilon_annot: Max abolute error in annotation.
 *   - epsilon_world: Max abolute error in meters of world coordinates.
 *
 * Could have done much more...
 * - Speed up by using an OpenMP loop.
 * - Larger read chunks and allow OpenGY to parallelize.
 * - Optionally output cube1 - cube2 as a ZGY file.
 * - Compare annotation.
 * - Compare lattice.
 */
void
compare(const std::string& file1, const std::string& file2,
        double epsilon_value,
        double epsilon_index,
        double epsilon_annot,
        double epsilon_world)
{
  OpenZGY::SeismicStoreIOContext rcontext;
  std::shared_ptr<OpenZGY::IZgyReader> r1 = OpenZGY::IZgyReader::open(file1, &rcontext);
  std::shared_ptr<OpenZGY::IZgyReader> r2 = OpenZGY::IZgyReader::open(file2, &rcontext);
  check_size_error(*r1, *r2);
  check_corners_error(*r1, *r2, epsilon_index, epsilon_annot, epsilon_world);
  check_sample_error(*r1, *r2, epsilon_value);
}

int main(int argc, const char **argv)
{
  try {
    if (argc < 3)
      throw std::runtime_error("Usage: zgydiff file1 file2 [epsilon] [-v|-q]");
    double epsilon_value{-1};
    double epsilon_index{-1};
    double epsilon_annot{-1};
    double epsilon_world{-1};
    if (argc > 3) {
      char *tail{nullptr};
      epsilon_value = strtod(argv[3], &tail);
      if (*tail)
        throw std::runtime_error("The third argument is not a number.");
    }
    if (epsilon_value >= 0) {
      // The limits for metadata errors are hard coded, but checks are
      // only enabled if the values are checked.
      epsilon_index = 0.0001;
      epsilon_annot = 0.01;
      epsilon_world = 0.1;
    }
    for (int ii=4; ii<argc; ++ii) {
      if (0 == strcmp(argv[ii], "-v"))
        ++verbose;
      else if(0 == strcmp(argv[ii], "-q"))
        --verbose;
      else
        throw std::runtime_error("Invalid command line option");
    }
    compare(argv[1], argv[2], epsilon_value, epsilon_index, epsilon_annot, epsilon_world);
    return 0;
  }
  catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }
}
