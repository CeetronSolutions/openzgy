// Copyright 2017-2023, SLB
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

#include "test_all.h"
#include "test_utils.h"
#include "mock.h"
#include "../writerargs.h"
#include "../impl/writerargsimpl.h"
#include "../impl/lodalgo.h"
#include "../impl/guid.h"
#include "../impl/enum_mapper.h"
#include "../impl/lodalgo.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <set>

using namespace OpenZGY;
using Test_Utils::must_throw;

namespace Test {
  class PrivateWriterArgsV1 : public ZgyWriterArgs
  {
  };

  class PrivateWriterArgsV2 : public ZgyWriterArgsV2
  {
  };

  // Might need to be declared as friend somewhere.
  class TestWriterArgs {
  public:
    static void test_v2_defaults();
    static void test_v2_overloads();
    static void test_v2_setters();
    static void test_v2_lodmode();
    static void test_v2_errors();
    static void test_v2_chain();
    static void test_v2_deepcopy();
    static void test_v2_guids();
    static void test_v2_setguids();
    static void test_v1_v2();
    static int which_overload(const ZgyWriterArgs&) { return 1; };
    static int which_overload(const ZgyWriterArgsV2&) { return 2; };
    static void test_enums();
  };
}

void Test::TestWriterArgs::test_v2_defaults()
{
  ZgyWriterArgsV2 args;
  TEST_EQUAL(args.impl().getInternalLodMode().incr.force, 1);
  TEST_EQUAL(args.impl().getInternalLodMode().incr.level, -1);
  TEST_EQUAL(args.impl().getInternalLodMode().last.force, 2);
  TEST_EQUAL(args.impl().getInternalLodMode().last.level, -1);
  TEST_CHECK(args.impl().getHistoRange()[0] > args.impl().getHistoRange()[1]);
  TEST_CHECK(args.impl().getDecimation().empty());
  TEST_CHECK(args.impl().getLodAlgorithm().empty());
  TEST_CHECK(args.impl().getDataId().empty());
  TEST_CHECK(args.impl().getVerId().empty());
  TEST_CHECK(args.impl().getPrevId().empty());
}

void Test::TestWriterArgs::test_v2_overloads()
{
  // Because the existing ZgyWriterArgs has no vtable, I will need
  // some tricky overload resolution where the compiler must choose
  // between ZgyWriterArgsV2 and its parent ZgyWriterArgs. Check that
  // this works.
  ZgyWriterArgs argsV1;
  ZgyWriterArgsV2 argsV2;
  PrivateWriterArgsV1 privateV1;
  PrivateWriterArgsV2 privateV2;

  TEST_EQUAL(which_overload(argsV1), 1);
  TEST_EQUAL(which_overload(argsV2), 2);
  TEST_EQUAL(which_overload(privateV1), 1);
  TEST_EQUAL(which_overload(privateV2), 2);
}

void Test::TestWriterArgs::test_v2_setters()
{
  ZgyWriterArgsV2 args;
  args
    .historange(-12, +42)
    .decimation(std::vector<OpenZGY::DecimationType>
                {OpenZGY::DecimationType::Maximum,
                 OpenZGY::DecimationType::Minimum});
  TEST_EQUAL(args.impl().getHistoRange()[0], -12);
  TEST_EQUAL(args.impl().getHistoRange()[1], 42);
  // Note that the setter in the API will use DecimationType;
  // I have not decided whether the data members and/or the getter
  // should convert those to the internal LodAlgorithm type.
  const auto decimation = args.impl().getDecimation();
  if (TEST_EQUAL(decimation.size(), 2)) {
    TEST_EQUAL((int)decimation[0], (int)OpenZGY::DecimationType::Maximum);
    TEST_EQUAL((int)decimation[1], (int)OpenZGY::DecimationType::Minimum);
  }
  //const auto algorithm = args.impl().getLodAlgorithm();
  //if (TEST_EQUAL(algorithm.size(), 2)) {
  //  TEST_EQUAL((int)algorithm[0], (int)InternalZGY::LodAlgorithm::Maximum);
  //  TEST_EQUAL((int)algorithm[1], (int)InternalZGY::LodAlgorithm::Minimum);
  //}
  if (verbose())
    args.dump(std::cout);
}

void Test::TestWriterArgs::test_v2_lodmode()
{
  auto modeEarly   = ZgyWriterArgsV2().lodmode(OpenZGY::LodMode::Early);
  auto modeEarly1  = ZgyWriterArgsV2().lodmode(OpenZGY::LodMode::Early1);
  auto modeNever   = ZgyWriterArgsV2().lodmode(OpenZGY::LodMode::Never);
  auto modeRebuild = ZgyWriterArgsV2().lodmode(OpenZGY::LodMode::Rebuild);
  auto modeDefault = ZgyWriterArgsV2().lodmode(OpenZGY::LodMode::Default);

  TEST_EQUAL(modeEarly.impl().getInternalLodMode().incr.force, 1);
  TEST_EQUAL(modeEarly.impl().getInternalLodMode().incr.level, -1);
  TEST_EQUAL(modeEarly.impl().getInternalLodMode().last.force, 2);
  TEST_EQUAL(modeEarly.impl().getInternalLodMode().last.level, -1);

  TEST_EQUAL(modeEarly1.impl().getInternalLodMode().incr.force, 1);
  TEST_EQUAL(modeEarly1.impl().getInternalLodMode().incr.level, 1);
  TEST_EQUAL(modeEarly1.impl().getInternalLodMode().last.force, 2);
  TEST_EQUAL(modeEarly1.impl().getInternalLodMode().last.level, -1);

  TEST_EQUAL(modeNever.impl().getInternalLodMode().incr.force, 0);
  TEST_EQUAL(modeNever.impl().getInternalLodMode().last.force, 0);

  TEST_EQUAL(modeRebuild.impl().getInternalLodMode().incr.force, 0);
  TEST_EQUAL(modeRebuild.impl().getInternalLodMode().last.force, 3);
  TEST_EQUAL(modeRebuild.impl().getInternalLodMode().last.level, -1);

  TEST_EQUAL(modeDefault.impl().getInternalLodMode().incr.force, 1);
  TEST_EQUAL(modeDefault.impl().getInternalLodMode().incr.level, -1);
  TEST_EQUAL(modeDefault.impl().getInternalLodMode().last.force, 2);
  TEST_EQUAL(modeDefault.impl().getInternalLodMode().last.level, -1);
}

void Test::TestWriterArgs::test_v2_errors()
{
  // If I remove lodIncrForce from the public API
  // then most of these tests go away. They might still
  // be settable by environment variables but I
  // probably don't need to test that.
  TEST_CHECK(must_throw("don't mix", [&](){
      auto args = ZgyWriterArgsV2()
        .lodmode(OpenZGY::LodMode::Never)
        .lodIncrForce(0);
      TEST_EQUAL(which_overload(args), 2);
  }));
  TEST_CHECK(must_throw("don't mix", [&](){
      auto args = ZgyWriterArgsV2()
        .lodIncrForce(0)
        .lodmode(OpenZGY::LodMode::Never);
      TEST_EQUAL(which_overload(args), 2);
  }));
  TEST_CHECK(must_throw("Invalid value", [&](){
      auto args = ZgyWriterArgsV2().lodIncrForce(2);
      TEST_EQUAL(which_overload(args), 2);
  }));
  TEST_CHECK(must_throw("Invalid value", [&](){
      auto args = ZgyWriterArgsV2().lodLastForce(1);
      TEST_EQUAL(which_overload(args), 2);
  }));
  TEST_CHECK(must_throw("Invalid value", [&](){
      auto args = ZgyWriterArgsV2().lodmode(static_cast<OpenZGY::LodMode>(999));
      TEST_EQUAL(which_overload(args), 2);
  }));
}

void Test::TestWriterArgs::test_v2_chain()
{
  // V2 settings cannot come after the original settings. This is
  // annoying but not serious because there will be a compilation error.
  // (fixed now)
  TEST_EQUAL(which_overload(ZgyWriterArgsV2().filename("foo").historange(-12, +42)), 2);

  // This is a more insiduous problem: The chain evaluates to type
  // ZgyWriterArgs, not V2. The called function has no way of knowing
  // that its input parameter can be cast to ZgyWriterArgsV2. The
  // settings from V1 WILL BE SILENTLY IGNORED. (fixed now)
  TEST_EQUAL(which_overload(ZgyWriterArgsV2().historange(-12, +42).filename("foo")), 2);

  // This is what the application currently needs to do.
  // Looks a lot less elegant.
  // (no longer required)
  ZgyWriterArgsV2 args;
  args.historange(-12, +42).filename("foo");
  TEST_EQUAL(which_overload(args), 2);

  // Done to fix this: Override all the setters in ZgyWriterArgs
  // in the ZgyWriterArgsV2 class so they won't return the base class.
}

namespace {
  std::string
  lookup(const std::string& data, const std::string& lookfor)
  {
    std::string result("(key not found)");
    std::size_t pos = data.find(lookfor);
    if (pos != std::string::npos) {
      std::size_t end = data.find("\n", pos);
      if (end == std::string::npos)
        result = data.substr(pos);
      else
        result = data.substr(pos, end-pos);
      result = result.substr(lookfor.size()); // strip key
    }
    std::size_t trim = result.find_first_not_of(" \t");
    if (trim != 0 && trim != std::string::npos)
      result = result.substr(trim); // strip leading spaces
    //std::cout << "Lookup \"" << lookfor << "\" -> \"" << result << "\"\n";
    return result;
  }
}

void Test::TestWriterArgs::test_v2_deepcopy()
{
  if (verbose())
    std::cout << std::endl;

  // Set a couple of args, both old and V2
  ZgyWriterArgsV2 args;
  args
    .historange(-12, +42)
    .decimation(std::vector<OpenZGY::DecimationType>
                {OpenZGY::DecimationType::Maximum,
                 OpenZGY::DecimationType::Minimum})
    .filename("one");

  // The copy should match exactly.
  ZgyWriterArgsV2 copy(args);
  std::stringstream args1, copy1;
  args.dump(args1);
  copy.dump(copy1);
  TEST_EQUAL(copy1.str(), args1.str());

  // Also verify the values that were set.
  TEST_EQUAL(lookup(args1.str(), "historange:"), "-12 to 42");
  TEST_EQUAL(lookup(args1.str(), "decimation:"), "105 104");
  TEST_EQUAL(lookup(args1.str(), "filename:"), "\"one\"");
  TEST_EQUAL(lookup(args1.str(), "zstart/inc:"), "0 / 0");

  // Change some values both in the old and the new parts.
  args
    .historange(-120, +420)
    .filename("two");
  copy
    .decimation(std::vector<OpenZGY::DecimationType>
                {OpenZGY::DecimationType::Average})
    .zstart(-999.25);
  std::stringstream args2, copy2;
  args.dump(args2);
  copy.dump(copy2);

  // Did the new values get set?
  TEST_EQUAL(lookup(args2.str(), "historange:"), "-120 to 420");
  TEST_EQUAL(lookup(copy2.str(), "decimation:"), "102");
  TEST_EQUAL(lookup(args2.str(), "filename:"), "\"two\"");
  TEST_EQUAL(lookup(copy2.str(), "zstart/inc:"), "-999.25 / 0");

  // Did values leak across instances?
  TEST_EQUAL(lookup(copy2.str(), "historange:"), "-12 to 42");
  TEST_EQUAL(lookup(args2.str(), "decimation:"), "105 104");
  TEST_EQUAL(lookup(copy2.str(), "filename:"), "\"one\"");
  TEST_EQUAL(lookup(args2.str(), "zstart/inc:"), "0 / 0");

  if (verbose()) {
    std::cout << "Initial args:\n" << args1.str()
              << "Initial copy:\n" << copy1.str()
              << "Updated args:\n" << args2.str()
              << "Updated copy:\n" << copy2.str()
              << "END" << std::endl;
  }
}

#ifndef _WIN32
/**
 * Unit test for class GUID. Put here because impl/guid.cpp doesn't
 * warrant its own test. And ZgyWriterArgs is one of the users.
 *
 * NOTE: the unit test may need to be commented out, because the
 * class being tested is aggressive about data hiding. Specifically,
 * the constructors are not tagged OPENZGY_API or even OPENZGY_TEST_API.
 * It is really important that application code won't be able to access
 * anything that deals with the internal representation and not string.
 *
 * Hiding isn't implemented on the Linux end. So it works there.
 */
void Test::TestWriterArgs::test_v2_guids()
{
  using InternalZGY::GUID;
  const GUID empty1(nullptr);
  const GUID empty2(GUID::guid_bytes_t{});
  const GUID normal("12345678-2345-6789-abcd-00ffABCDEF01");
  const GUID empty3("00000000-0000-0000-0000-000000000000");
  const GUID empty4("");
  const GUID random;
  const GUID r_normal(normal.toString());
  const GUID r_random(random.toString());
  GUID::guid_bytes_t raw;
  normal.copyTo(raw.data(), raw.size());

  const GUID check1(raw);
  TEST_EQUAL(empty1.toString(), "00000000-0000-0000-0000-000000000000");
  TEST_EQUAL(empty2.toString(), "00000000-0000-0000-0000-000000000000");
  TEST_EQUAL(empty3.toString(), "00000000-0000-0000-0000-000000000000");
  TEST_EQUAL(empty4.toString(), "00000000-0000-0000-0000-000000000000");
  TEST_EQUAL(normal.toString(), "12345678-2345-6789-abcd-00ffabcdef01");
  TEST_EQUAL(r_normal.toString(), "12345678-2345-6789-abcd-00ffabcdef01");
  TEST_CHECK(r_random.toString() != normal.toString());
  TEST_CHECK(r_random.toString() != empty1.toString());
  TEST_EQUAL(check1.toString(), "12345678-2345-6789-abcd-00ffabcdef01");
}
#endif

/**
 * Unit test for ZgyWriterArgs able to pick up guids from another file.
 * See also api.setguids for whether the guods will then get propagated
 * to the file.
 */
void Test::TestWriterArgs::test_v2_setguids()
{
  if (verbose())
    std::cout << std::endl;

  std::array<std::int64_t,3> size{7,9,13};
  std::shared_ptr<IZgyReader> r = Test::ZgyReaderMock::mock(size, SampleDataType::float32);
  TEST_EQUAL(r->verid(), "0000000a-000b-000c-000d-222222222222");

  // Check that the argument package picks up the guid.
  // Whether it actually is used needs to be checked
  // in a higher level test that is allowed to create data.
  ZgyWriterArgsV2 args = ZgyWriterArgsV2().guidsfrom(r);
  if (verbose())
    args.dump(std::cout);
  TEST_EQUAL(args.impl().getDataId(), "");
  TEST_EQUAL(args.impl().getVerId(), r->verid());
  TEST_EQUAL(args.impl().getPrevId(), "");

  // Check that merge() propagates the guid if it was
  // explicity set (by guidsfrom).
  ZgyWriterArgsV2 a2 = ZgyWriterArgsV2().merge(args);
  TEST_EQUAL(a2.impl().getDataId(), "");
  TEST_EQUAL(a2.impl().getVerId(), r->verid());
  TEST_EQUAL(a2.impl().getPrevId(), "");

  // Check that merge() will not propagate the (empty)
  // guids and thus overwrite guids already set.
  ZgyWriterArgsV2 a3 = args.merge(ZgyWriterArgsV2());

  TEST_EQUAL(a2.impl().getDataId(), "");
  TEST_EQUAL(a2.impl().getVerId(), r->verid());
  TEST_EQUAL(a2.impl().getPrevId(), "");
}

void Test::TestWriterArgs::test_v1_v2()
{
  if (verbose())
    std::cout << std::endl;

  // Completely bogus, not even remotely orthogonal.
  const ZgyWriterArgs::corners_t corners{{{1,2},{3,4},{5,6},{7,8}}};

  ZgyWriterArgs v1;
  v1
    .filename("yahoo.zgy")
    .size(111,222,333)
    .bricksize(16,32,128)
    .datatype(SampleDataType::int16)
    .datarange(-1001, +5432)
    .hunit(UnitDimension::length, "toes", 0.19)
    .zunit(UnitDimension::time, "days", 998.0)
    .ilstart(202)
    .ilinc(2)
    .xlstart(303)
    .xlinc(3)
    .zstart(404)
    .zinc(4)
    .corners(corners);

  ZgyWriterArgsV2 v2;
  v2
    .filename("yahoo.zgy")
    .size(111,222,333)
    .bricksize(16,32,128)
    .datatype(SampleDataType::int16)
    .datarange(-1001, +5432)
    .hunit(UnitDimension::length, "toes", 0.19)
    .zunit(UnitDimension::time, "days", 998.0)
    .ilstart(202)
    .ilinc(2)
    .xlstart(303)
    .xlinc(3)
    .zstart(404)
    .zinc(4)
    .corners(corners)
    .decimation(std::vector<DecimationType>
                {DecimationType::MostFrequent,DecimationType::Average})
    .lodmode(LodMode::Early1);

  ZgyWriterArgsV3 v3;
  v3
    .filename("yahoo.zgy")
    .size(111,222,333)
    .bricksize(16,32,128)
    .datatype(SampleDataType::int16)
    .datarange(-1001, +5432)
    .hunit(UnitDimension::length, "toes", 0.19)
    .zunit(UnitDimension::time, "days", 998.0)
    .ilstart(202)
    .ilinc(2)
    .xlstart(303)
    .xlinc(3)
    .zstart(404)
    .zinc(4)
    .corners(corners)
    .decimation(std::vector<DecimationType>
                {DecimationType::MostFrequent,DecimationType::Average})
    .lodmode(LodMode::Early1);

  std::stringstream dump1;
  v1.dump(dump1);
  std::stringstream dump2;
  v2.dump(dump2);
  std::stringstream dump3;
  v3.dump(dump3);

  if (verbose()) {
    std::cout << "V1 arguments (" << dump1.str().size() << " bytes)\n"
              << dump1.str()
              << "\nV2 arguments (" << dump2.str().size() << " bytes)\n"
              << dump2.str()
              << "\nV3 arguments (" << dump3.str().size() << " bytes)\n"
              << dump3.str()
              << std::flush;
  }

  std::string check1 = dump1.str();
  std::string check2 = dump2.str().substr(0, check1.size());
  // The V2 struct has more members, but the first parts should match.
  TEST_EQUAL(check1, check2);

  // V3 has nothing new.
  TEST_EQUAL(dump3.str(), dump2.str());

  // Convert a V1 instance to V3 explicitly.
  ZgyWriterArgsV3 v31(v1, ZgyWriterArgsV3::FromOldWriterArgs{});
  std::stringstream dump31;
  v31.dump(dump31);
  TEST_EQUAL(dump31.str().substr(0, dump1.str().size()), dump1.str());

  // Convert a V2 instance to V3 explicitly.
  ZgyWriterArgsV3 v32(v2, ZgyWriterArgsV3::FromOldWriterArgsV2{});
  std::stringstream dump32;
  v32.dump(dump32);
  TEST_EQUAL(dump32.str(), dump3.str());
}

void Test::TestWriterArgs::test_enums()
{
  using OpenZGY::Impl::EnumMapper;
  using OpenZGY::SampleDataType;
  using OpenZGY::UnitDimension;
  using OpenZGY::DecimationType;
  using InternalZGY::RawDataType;
  using InternalZGY::RawHorizontalDimension;
  using InternalZGY::RawVerticalDimension;
  using InternalZGY::LodAlgorithm;

  std::vector<SampleDataType> known_dt
    {//SampleDataType::unknown,
     SampleDataType::int8,
     SampleDataType::int16,
     SampleDataType::float32};
  for (SampleDataType dt : known_dt) {
    const RawDataType mapped =
      EnumMapper::mapSampleDataTypeToRawDataType(dt);
    const SampleDataType check =
      EnumMapper::mapRawDataTypeToSampleDataType(mapped);
    TEST_EQUAL((int)check, (int)dt);
    TEST_CHECK((int)dt != (int)mapped);
  }

  std::vector<UnitDimension> known_dim
    {UnitDimension::unknown,
     UnitDimension::time,
     UnitDimension::length,
     UnitDimension::arcangle};
  for (UnitDimension dim : known_dim) {
    if (dim != UnitDimension::time) {
      const RawHorizontalDimension mapped_h =
        EnumMapper::mapUnitDimensionToRawHorizontalDimension(dim);
      const UnitDimension check_h =
        EnumMapper::mapRawHorizontalDimensionToUnitDimension(mapped_h);
      TEST_EQUAL((int)check_h, (int)dim);
      TEST_CHECK((int)dim != (int)mapped_h);
    }
    else {
      TEST_CHECK(must_throw("time is not allowed", [&](){
        (void)EnumMapper::mapUnitDimensionToRawHorizontalDimension(dim);
      }));
    }
    if (dim != UnitDimension::arcangle) {
      const RawVerticalDimension mapped_v =
        EnumMapper::mapUnitDimensionToRawVerticalDimension(dim);
      const UnitDimension check_v =
        EnumMapper::mapRawVerticalDimensionToUnitDimension(mapped_v);
      TEST_EQUAL((int)check_v, (int)dim);
      TEST_CHECK((int)dim != (int)mapped_v);
    }
    else {
      TEST_CHECK(must_throw("arcangle is not allowed", [&](){
        (void)EnumMapper::mapUnitDimensionToRawVerticalDimension(dim);
      }));
    }
  }

  std::vector<DecimationType> known_decim
    {DecimationType::LowPass,
     DecimationType::LowPassNew,
     DecimationType::WeightedAverage,
     DecimationType::Average,
     DecimationType::Median,
     DecimationType::Minimum,
     DecimationType::Maximum,
     DecimationType::Decimate,
     DecimationType::DecimateSkipNaN,
     DecimationType::AllZero,
     DecimationType::MostFrequent,
     DecimationType::MostFrequentNon0,
     DecimationType::AverageNon0};
  const std::vector<LodAlgorithm> all_algo =
    EnumMapper::mapDecimationTypeToLodAlgorithm(known_decim);
  TEST_EQUAL(all_algo.size(), known_decim.size());
  std::set<LodAlgorithm> all_algo_set;
  for (auto it : all_algo)
    all_algo_set.insert(it);
  TEST_EQUAL(all_algo_set.size(), all_algo.size()); // All unique?
}

namespace {
  class Register
  {
  public:
    Register()
    {
      register_test("writerargs.v2.defaults",  Test::TestWriterArgs::test_v2_defaults);
      register_test("writerargs.v2.overloads", Test::TestWriterArgs::test_v2_overloads);
      register_test("writerargs.v2.setters",   Test::TestWriterArgs::test_v2_setters);
      register_test("writerargs.v2.lodmode",   Test::TestWriterArgs::test_v2_lodmode);
      register_test("writerargs.v2.errors",    Test::TestWriterArgs::test_v2_errors);
      register_test("writerargs.v2.chain",     Test::TestWriterArgs::test_v2_chain);
      register_test("writerargs.v2.deepcopy",  Test::TestWriterArgs::test_v2_deepcopy);
#ifndef _WIN32
      register_test("writerargs.v2.guids",     Test::TestWriterArgs::test_v2_guids);
#endif
      register_test("writerargs.v2.setguids",  Test::TestWriterArgs::test_v2_setguids);
      register_test("writerargs.v1-v2",        Test::TestWriterArgs::test_v1_v2);
      register_test("writerargs.enums",        Test::TestWriterArgs::test_enums);
    }
  } dummy;
} // namespace for registration
