// Copyright 2017-2020, Schlumberger
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

#include <openzgy/api.h>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

/**
 * \file example.h
 * \brief Example program using the C++ API.
 * \details To build this e.g. on Ububtu focal:
 * Compile: \code g++ -std=c++11 -Iinclude -oexample -x c++ example.h -Lfocal-gcc93 -Wl,-rpath-link=focal-gcc93 -lopenzgy \endcode
 * Run: \code env LD_LIBRARY_PATH=focal-gcc93 ./example fromfile.zgy tofile.zgy \endcode
 */

void copy(const std::string& srcname, const std::string& dstname)
{
  using namespace OpenZGY;
  ProgressWithDots p1, p2;
  std::shared_ptr<IZgyReader> r = IZgyReader::open(srcname);
  std::shared_ptr<IZgyWriter> w = IZgyWriter::open(ZgyWriterArgs().metafrom(r).filename(dstname));
  const std::array<std::int64_t,3> size = r->size();
  const std::array<std::int64_t,3> brick = r->bricksize();
  const std::array<std::int64_t,3> bs{brick[0], brick[1], size[2]};
  const std::int64_t total = ((size[0] + bs[0] - 1) / bs[0]) *
                             ((size[1] + bs[1] - 1) / bs[1]);
  std::int64_t done{0};
  std::unique_ptr<float> buf(new float[bs[0]*bs[1]*bs[2]]);
  std::array<std::int64_t,3> pos;
  for (pos[0] = 0; pos[0] < size[0]; pos[0] += bs[0]) {
    for (pos[1] = 0; pos[1] < size[1]; pos[1] += bs[1]) {
      for (pos[2] = 0; pos[2] < size[2]; pos[2] += bs[2]) {
        r->read(pos, bs, buf.get(), 0);
        w->write(pos, bs, buf.get());
        p1(++done, total);
      }
    }
  }
  w->finalize(std::vector<DecimationType>(), std::ref(p2));
  w->close();
}

int main(int argc, const char **argv)
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " infile outfile" << std::endl;
    exit(1);
  }
  try {
    copy(argv[1], argv[2]);
  }
  catch (const std::exception& ex) {
    std::cerr << argv[0] << ": " << ex.what() << std::endl;
    exit(1);
  }
}
