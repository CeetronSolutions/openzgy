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

#include "file_parallelizer.h"
#include "../exception.h"

namespace InternalZGY {
#if 0
}
#endif

FileRelayBase::FileRelayBase(std::shared_ptr<IFileBase> relay)
  : _relay(relay)
{
  if (!_relay)
    throw OpenZGY::Errors::ZgyInternalError("FileRelayBase created with no target");
}

FileRelayBase::~FileRelayBase()
{
}

FileRelay::FileRelay(std::shared_ptr<IFileADT> relay)
  : FileRelayBase(relay),_relay(relay)
{
  if (!_relay)
    throw OpenZGY::Errors::ZgyInternalError("FileRelay created with no target");
}

FileRelay::~FileRelay()
{
}

} // namespace
