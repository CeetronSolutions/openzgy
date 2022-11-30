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

#pragma once

#if defined OPENZGY_STATIC || !defined _WIN32
    #define OPENZGY_API
    #define OPENZGY_TEST_API
    #define OPENZGY_DECLARE_EXPLICIT_TEMPLATE(...)   extern template class __VA_ARGS__;
    #define OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(...) template class  __VA_ARGS__;
#else
    #ifdef OPENZGY_DLL
        #define OPENZGY_API      __declspec(dllexport) // this is the public API
        #define OPENZGY_TEST_API __declspec(dllexport) // exported only for unit tests
        #define OPENZGY_DECLARE_EXPLICIT_TEMPLATE(...)
        #define OPENZGY_IMPLEMENT_EXPLICIT_TEMPLATE(...) template class __declspec(dllexport) __VA_ARGS__;
    #else
        #define OPENZGY_API      __declspec(dllimport)
        #define OPENZGY_TEST_API __declspec(dllimport)
        #define OPENZGY_DECLARE_EXPLICIT_TEMPLATE(...) template class __declspec(dllimport) __VA_ARGS__;
    #endif
#endif
