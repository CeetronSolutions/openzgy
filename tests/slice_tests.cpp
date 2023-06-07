/////////////////////////////////////////////////////////////////////////////////
//
// Copyright 2023 Equinor ASA
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http ://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
/////////////////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"

#include <algorithm>
#include <string>
#include <cmath>

#include "zgyaccess/seismicslice.h"

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(slice_tests, testSlice)
{
    ZGYAccess::SeismicSliceData slice(200, 133);

    ASSERT_EQ(slice.size(), 200*133);
    ASSERT_EQ(slice.width(), 200);
    ASSERT_EQ(slice.depth(), 133);
    ASSERT_FALSE(slice.isEmpty());

    ASSERT_NE(slice.values(), nullptr);

    float* pData = slice.values();

    for (int i = 0; i < slice.size(); i++, pData++)
    {
        *pData = 1.0f * i;
    }

    ASSERT_EQ(slice.valueAt(10, 10), 1340.0);

    slice.reset();
    ASSERT_EQ(slice.size(), 0);
    ASSERT_EQ(slice.width(), 0);
    ASSERT_EQ(slice.depth(), 0);
    ASSERT_EQ(slice.values(), nullptr);

    ASSERT_TRUE(slice.isEmpty());

}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(slice_tests, testMute) {
    ZGYAccess::SeismicSliceData slice(60, 90);

    float *pData = slice.values();

    for (int i = 0; i < slice.size(); i++, pData++) 
    {
        *pData = 1.0f * i;
    }

    slice.mute(100.0);

    pData = slice.values();

    for (int i = 0; i < slice.size(); i++, pData++)
    {
        ASSERT_TRUE((*pData == 0.0f) || (*pData >= 100.0f));
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(slice_tests, testLimiter) {
    ZGYAccess::SeismicSliceData slice(120, 77);

    float* pData = slice.values();

    int offset = slice.size() / 2;

    for (int i = 0; i < slice.size(); i++, pData++)
    {
        *pData = 1.0f * (i - offset);
    }

    slice.limitTo(-40.f, 50.f);

    pData = slice.values();

    for (int i = 0; i < slice.size(); i++, pData++)
    {
        ASSERT_TRUE((*pData >= -40.0f) || (*pData <= 50.0f));
    }
}
