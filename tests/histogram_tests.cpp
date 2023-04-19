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

#include "zgyaccess/zgy_histogram.h"

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(histogram_tests, testHistogram)
{

    std::vector<float> testdata = { -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 1.0, 1.0 };

    auto hist = ZGYAccess::HistogramGenerator::getHistogram(testdata, 5, -1.1f, 1.1f);

    ASSERT_EQ(hist->Xvalues.size(), 5);
    ASSERT_EQ(hist->Yvalues.size(), 5);

    ASSERT_EQ(hist->Yvalues[0], 2);
    ASSERT_EQ(hist->Yvalues[1], 0);
    ASSERT_EQ(hist->Yvalues[2], 5);
    ASSERT_EQ(hist->Yvalues[3], 1);
    ASSERT_EQ(hist->Yvalues[4], 2);
}
