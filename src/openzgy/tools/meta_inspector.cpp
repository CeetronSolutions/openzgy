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

std::string toString(const SampleDataType& dataType) {
  switch (dataType) {
    case SampleDataType::int8: return "int8";
    case SampleDataType::int16: return "int16";
    case SampleDataType::float32: return "float32";
    default: return "unknown";
  }
}

std::string toString(const UnitDimension& dimension) {
  switch (dimension) {
    case UnitDimension::arcangle: return "arcangle";
    case UnitDimension::time: return "time";
    case UnitDimension::length: return "length";
    default: return "unknown";
  }
}

std::string toString(const std::vector<IZgyMeta::size3i_t>& sizes) {
  std::ostringstream out;
  for (IZgyMeta::size3i_t size : sizes) out << "[" << size[0] << ", " << size[1] << ", " << size[2] << "] ";
  return out.str();
}

std::string toString(const std::vector<std::int64_t>& bins) {
  std::ostringstream out; 
  for (std::int64_t bin : bins) out << bin << " ";
  return out.str();
}

std::string showMeta(const std::shared_ptr<OpenZGY::IZgyReader>& reader){
  std::ostringstream out;
  out << "Metadata\n"
    << "  Size: [" << reader->size()[0] << ", " << reader->size()[1] << ", " << reader->size()[2] << "]\n"
    << "  Brick size: [" << reader->bricksize()[0] << ", " << reader->bricksize()[1] << ", " << reader->bricksize()[2] << "]\n"
    << "  Data type: " << toString(reader->datatype()) << "\n"
    << "  Data range: " << reader->datarange()[0] << " <-> " << reader->datarange()[1] << "\n"
    << "  Z unit dimension: " << toString(reader->zunitdim()) << "\n"
    << "  Z unit name: " << reader->zunitname() << "\n"
    << "  Z unit factor: " << reader->zunitfactor() << "\n"
    << "  Z slice start: " << reader->zstart() << "\n"
    << "  Z slice increment: " << reader->zinc() << "\n"
    << "  XY unit dimension: " << toString(reader->hunitdim()) << "\n"
    << "  XY unit name: " << reader->hunitname() << "\n"
    << "  XY unit factor: " << reader->hunitfactor() << "\n"
    << "  Inline start: " << reader->annotstart()[0] << "\n"
    << "  Inline increment: " << reader->annotinc()[0] << "\n"
    << "  Crossline start: " << reader->annotstart()[1] << "\n"
    << "  Crossline increment: " << reader->annotinc()[1] << "\n"
    << "  World corner 0: " << reader->corners()[0][0] << ", " << reader->corners()[0][1] << "\n"
    << "  World corner 1: " << reader->corners()[1][0] << ", " << reader->corners()[1][1] << "\n"
    << "  World corner 2: " << reader->corners()[2][0] << ", " << reader->corners()[2][1] << "\n"
    << "  World corner 3: " << reader->corners()[3][0] << ", " << reader->corners()[3][1] << "\n"
    << "  Index corner 0: " << reader->indexcorners()[0][0] << ", " << reader->indexcorners()[0][1] << "\n"
    << "  Index corner 1: " << reader->indexcorners()[1][0] << ", " << reader->indexcorners()[1][1] << "\n"
    << "  Index corner 2: " << reader->indexcorners()[2][0] << ", " << reader->indexcorners()[2][1] << "\n"
    << "  Index corner 3: " << reader->indexcorners()[3][0] << ", " << reader->indexcorners()[3][1] << "\n"
    << "  Annotation corner 0: " << reader->annotcorners()[0][0] << ", " << reader->annotcorners()[0][1] << "\n"
    << "  Annotation corner 1: " << reader->annotcorners()[1][0] << ", " << reader->annotcorners()[1][1] << "\n"
    << "  Annotation corner 2: " << reader->annotcorners()[2][0] << ", " << reader->annotcorners()[2][1] << "\n"
    << "  Annotation corner 3: " << reader->annotcorners()[3][0] << ", " << reader->annotcorners()[3][1] << "\n"
    << "  Number of LODs: " << reader->nlods() << "\n"
    << "  Bricks per LOD: " << toString(reader->brickcount()) << "\n"    ;

  return out.str();
}

std::string showHistogram(const std::shared_ptr<OpenZGY::IZgyReader>& reader) {
  std::ostringstream out;
  const SampleHistogram& sh = reader->histogram();
  out << "Histogram\n"
    << "  Center of the first bin: " << sh.minvalue << "\n"
    << "  Center of the last bin: " << sh.maxvalue << "\n"
    << "  Sum of counts of all bins: " << sh.samplecount << "\n"
    << "  Counts by bin: " << toString(sh.bins) << "\n";
  return out.str();
}

std::string showStatistics(const std::shared_ptr<OpenZGY::IZgyReader>& reader) {
  std::ostringstream out;

  const SampleStatistics& ss = reader->statistics();
  out << "Statistics\n"
    << "  Min: " << ss.min << "\n"
    << "  Max: " << ss.max << "\n"
    << "  Sum: " << ss.sum << "\n"
    << "  Sum of squares: " << ss.ssq << "\n";
  return out.str();
}

int main(int argc, const char** argv)
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " zgyfile" << std::endl;
    exit(1);
  }
  try {
    std::shared_ptr<IZgyReader> reader = IZgyReader::open(argv[1], getContext().get());
    std::cout << showMeta(reader) << std::endl;
    std::cout << showStatistics(reader) << std::endl;
    std::cout << showHistogram(reader) << std::endl;
  }
  catch (const std::exception& ex) {
    std::cerr << argv[0] << ": " << ex.what() << std::endl;
    exit(1);
  }
}
