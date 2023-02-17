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

// Minor tweaks to the original "cutest.h",
// allowing tests in multiple source files to execute as one.
#pragma once
#ifndef CUTEST_HAS_MAIN
#define TEST_NO_MAIN
#endif
#define NOMINMAX // cutest includes Windows.h. Ugh!
#include "cutest.h"
extern void register_test(const char* name, void (*func)(void));
extern void register_sd_test(const char* name, void (*func)(void));
extern int verbose();
