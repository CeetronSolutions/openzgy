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
#include <fstream>
#include <string>
#include <variant>

#include "zgyaccess/zgyreader.h"
#include "testdatafolder.h"

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testOpenFile)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto metadata = reader.metaData();

    ASSERT_EQ(metadata.size(), 21);


    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testFailOpenFile)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_FALSE(reader.open(std::string(TEST_DATA_DIR) + "does_not_exist.zgy"));

    ASSERT_EQ(reader.inlineSize(), 0);
    ASSERT_EQ(reader.xlineSize(), 0);
    ASSERT_EQ(reader.zSize(), 0);

}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testSeismicSize)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    ASSERT_EQ(reader.inlineSize(), 112);
    ASSERT_EQ(reader.xlineSize(), 64);
    ASSERT_EQ(reader.zSize(), 176);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testHistogram)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    ZGYAccess::HistogramData* hist = reader.histogram();

    ASSERT_EQ(hist->Xvalues.size(), 256);
    ASSERT_EQ(hist->Xvalues.size(), hist->Yvalues.size());

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testMinMaxDataValue)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [minval, maxval] = reader.dataRange();

    ASSERT_DOUBLE_EQ(minval, -28);
    ASSERT_DOUBLE_EQ(maxval, 227);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testSeismicOutline)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    ZGYAccess::Outline outline = reader.seismicWorldOutline();

    ASSERT_TRUE(outline.isValid());
    ASSERT_EQ(outline.points().size(), 4);

    ASSERT_DOUBLE_EQ(outline.points()[0].x(), 1000.0);
    ASSERT_DOUBLE_EQ(outline.points()[0].y(), 1000.0);

    ASSERT_TRUE(outline.points()[3] == ZGYAccess::Point2d(3775.0, 2890.0));

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testZRange)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [zmin, zmax] = reader.zRange();
    double ztep = reader.zStep();

    ASSERT_DOUBLE_EQ(ztep, 4.125);
    ASSERT_DOUBLE_EQ(zmin, 2500);
    ASSERT_DOUBLE_EQ(zmax, 3226);


    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testWorldCoord)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [inlineFrom, inlineTo] = reader.inlineRange();
    int inlineStep = reader.inlineStep();
    auto [xlineFrom, xlineTo] = reader.xlineRange();
    int xlineStep = reader.xlineStep();

    auto [wX, wY] = reader.toWorldCoordinate(inlineFrom, xlineFrom);
    auto [w2X, w2Y] = reader.toWorldCoordinate(inlineTo, xlineFrom);

    ASSERT_EQ(inlineStep, 5);
    ASSERT_EQ(xlineStep, 2);

    ASSERT_EQ(inlineFrom, 1234);
    ASSERT_EQ(inlineTo, 1794);

    ASSERT_EQ(xlineFrom, 5678);
    ASSERT_EQ(xlineTo, 5806);

    ASSERT_DOUBLE_EQ(wX, 1000.0);
    ASSERT_DOUBLE_EQ(w2X, 3800.0);
    ASSERT_DOUBLE_EQ(wY, 1000.0);
    ASSERT_DOUBLE_EQ(w2Y, 1000.0);


    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadInlineSlice)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.inlineSlice(50);

    double value = 0.0;

    float* values = data->values();

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadInlineSliceSubset)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.inlineSlice(50, 15, 75);

    double value = 0.0;

    float* values = data->values();

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadInlineSliceFail)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.inlineSlice(500);

    float* values = data->values();

    ASSERT_EQ(data->size(), 0);
    ASSERT_EQ(values, nullptr);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadXlineSlice)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.xlineSlice(20);

    double value = 0.0;

    float* values = data->values();

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadXlineSliceSubset)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.xlineSlice(30, 40, 5);

    double value = 0.0;

    float* values = data->values();

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadZSlice)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.zSlice(30);

    double value = 0.0;

    float* values = data->values();

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadZTrace)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.zTrace(20,22);

    double value = 0.0;

    float* values = data->values();

    ASSERT_TRUE(values != nullptr);

    ASSERT_EQ(data->size(), reader.zSize());

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testReadZTraceSubset)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto data = reader.zTrace(20, 22, 10, 20);

    double value = 0.0;

    float* values = data->values();

    ASSERT_TRUE(values != nullptr);

    ASSERT_EQ(data->size(), 20);

    for (int i = 0; i < data->size(); i++)
    {
        if (values[i] != 0.0) value = values[i];
    }

    ASSERT_NE(value, 0.0);

    reader.close();
}


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testWorldtoIndexCoord)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [inlineFrom, inlineTo] = reader.inlineRange();
    auto [xlineFrom, xlineTo] = reader.xlineRange();

    auto [wX, wY] = reader.toWorldCoordinate(inlineFrom, xlineFrom);
    auto [w2X, w2Y] = reader.toWorldCoordinate(inlineTo, xlineTo);

    auto [convInlineFrom, convXlineFrom] = reader.toInlineXline(wX, wY);
    auto [convInlineTo, convXlineTo] = reader.toInlineXline(w2X, w2Y);

    ASSERT_EQ(inlineFrom, convInlineFrom);
    ASSERT_EQ(inlineTo, convInlineTo);
    ASSERT_EQ(xlineFrom, convXlineFrom);
    ASSERT_EQ(xlineTo, convXlineTo);

    reader.close();
}
