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

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "../api.h"
#include "../iocontext.h"
#include "../impl/environment.h"

using namespace OpenZGY;
using InternalZGY::Environment;

const std::shared_ptr<IOContext> getContext(){
  return Environment::getStringEnv("OPENZGY_SDURL") != ""
    && Environment::getStringEnv("OPENZGY_SDAPIKEY") != ""
    && Environment::getStringEnv("OPENZGY_TOKEN") != ""
    ? std::shared_ptr<IOContext>(new SeismicStoreIOContext()) : std::shared_ptr<IOContext>();
}

template<typename T>
T coerceToRange(T value, T min, T max) {
  return std::min(max, std::max(min, value));
}

template<class T, class C> 
void writeToAsciiFile(std::string filename, const std::shared_ptr<T>& buffer, int64_t bufferSize) {
  std::fstream file(filename, std::ios::out);
  ProgressWithDots progress(42, std::cout);
  std::cout << "Writing output " << std::endl;
  for (int i = 0; i < bufferSize; i++) {
    progress(i, bufferSize - 1);
    file << (C)buffer.get()[i] << " ";
  }

  file.close();
}

double elapsedTime(const std::chrono::high_resolution_clock::time_point& startTime) {
  auto endTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double, std::micro>(endTime - startTime);
  return duration.count() / 1000000.0F;
}

int sampleSize(const std::shared_ptr<IZgyReader>& reader) {
  return 1 << (static_cast<int>(reader->datatype()) - static_cast<int>(SampleDataType::int8));
}

void printBenchmark(const int64_t& bufferSize, std::shared_ptr<OpenZGY::IZgyReader>& reader, double elapsedSeconds)
{
  std::cout
    << "  MBytes read: " << bufferSize * sampleSize(reader) / (1048576.0) << std::endl
    << "  Reading time: " << elapsedSeconds << "s" << std::endl
    << "  Throughput: " << (bufferSize * sampleSize(reader) / (1048576.0)) / elapsedSeconds << "MB/s" << std::endl;
}

int main(int argc, const char** argv){
  if (argc < 9) {
    std::cerr << "Usage: " << argv[0] << " zgy_file i0 j0 k0 i1 j1 k1 lod [output_ascii_file]" << std::endl;
    exit(1);
  }
  try {
    bool genOutput = argc == 10;
    std::cout << "Loading...";
    std::shared_ptr<IZgyReader> reader = IZgyReader::open(argv[1], getContext().get());
    std::cout << "Done" << std::endl;

    int64_t i0 = coerceToRange<int64_t>(strtoll(argv[2], nullptr, 0), 0, reader->size()[0] - 1);
    int64_t j0 = coerceToRange<int64_t>(strtoll(argv[3], nullptr, 0), 0, reader->size()[1] - 1);
    int64_t k0 = coerceToRange<int64_t>(strtoll(argv[4], nullptr, 0), 0, reader->size()[2] - 1);
    int64_t i1 = coerceToRange<int64_t>(strtoll(argv[5], nullptr, 0), 0, reader->size()[0] - 1);
    int64_t j1 = coerceToRange<int64_t>(strtoll(argv[6], nullptr, 0), 0, reader->size()[1] - 1);
    int64_t k1 = coerceToRange<int64_t>(strtoll(argv[7], nullptr, 0), 0, reader->size()[2] - 1);
    int lod = coerceToRange(atoi(argv[8]), 0, reader->nlods() - 1);

    if (i0 > i1) std::swap(i0, i1);
    if (j0 > j1) std::swap(j0, j1);
    if (k0 > k1) std::swap(k0, k1);
    std::array<std::int64_t, 3> start { i0, j0, k0 };
    std::array<std::int64_t, 3> size { i1 - i0, j1 - j0, k1 - k0 };

    auto startTime = std::chrono::high_resolution_clock::now();

    int64_t bufferSize = size[0] * size[1] * size[2];
    std::cout << "Reading requested area..." << std::endl;

    if (reader->datatype() == SampleDataType::int8){
      std::shared_ptr<std::int8_t> buf8(new std::int8_t[bufferSize]);
      reader->read(start, size, buf8.get(), lod);
      printBenchmark(bufferSize, reader, elapsedTime(startTime));
      if (genOutput) writeToAsciiFile<std::int8_t, int>(argv[9], buf8, bufferSize);
    }
    else if (reader->datatype() == SampleDataType::int16) {
      std::shared_ptr<std::int16_t> buf16(new std::int16_t[bufferSize]);
      reader->read(start, size, buf16.get(), lod);
      printBenchmark(bufferSize, reader, elapsedTime(startTime));
      if (genOutput) writeToAsciiFile<std::int16_t, int>(argv[9], buf16, bufferSize);
    }
    else if (reader->datatype() == SampleDataType::float32) {
      std::shared_ptr<float> buf32(new float[bufferSize]);
      reader->read(start, size, buf32.get(), lod);
      printBenchmark(bufferSize, reader, elapsedTime(startTime));
      if (genOutput) writeToAsciiFile<float, float>(argv[9], buf32, bufferSize);
    }
  }
  catch (const std::exception& ex) {
    std::cerr << argv[0] << ": " << ex.what() << std::endl;
    exit(1);
  }
}

