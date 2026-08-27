// Copyright 2017-2024, Schlumberger
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

#include "mock.h"

#include "../writerargs.h"
#include "../impl/writerargsimpl.h"
#include "../impl/enum_mapper.h"

Test::ZgyMetaMock::~ZgyMetaMock()
{
}

Test::ZgyReaderMock::~ZgyReaderMock()
{
}

Test::ZgyWriterMock::~ZgyWriterMock()
{
}

std::shared_ptr<OpenZGY::IZgyWriter>
Test::ZgyWriterMock::mock(const OpenZGY::ZgyWriterArgsV3& args)
{
  using OpenZGY::Impl::EnumMapper;
  const OpenZGY::SampleDataType dt =
    EnumMapper::mapRawDataTypeToSampleDataType(args.impl().datatype_);
  const auto size = args.impl().size_;
  return mock(size, dt);
}
