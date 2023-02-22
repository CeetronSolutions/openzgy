

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

    ASSERT_TRUE(reader.Open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto metadata = reader.MetaData();

    ASSERT_EQ(metadata.size(), 21);


    reader.Close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testFailOpenFile)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_FALSE(reader.Open(std::string(TEST_DATA_DIR) + "does_not_exist.zgy"));
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testHistogram)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.Open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    ZGYAccess::HistogramData* hist = reader.histogram();

    ASSERT_EQ(hist->Xvalues.size(), 256);
    ASSERT_EQ(hist->Xvalues.size(), hist->Yvalues.size());

    reader.Close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testMinMaxDataValue)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.Open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [minval, maxval] = reader.DataRange();

    ASSERT_DOUBLE_EQ(minval, -28);
    ASSERT_DOUBLE_EQ(maxval, 227);

    reader.Close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testSeismicOutline)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.Open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    ZGYAccess::Outline outline = reader.seismicOutline();

    ASSERT_TRUE(outline.isValid());
    ASSERT_EQ(outline.points().size(), 4);

    ASSERT_DOUBLE_EQ(outline.points()[0].x(), 1000.0);
    ASSERT_DOUBLE_EQ(outline.points()[0].y(), 1000.0);

    ASSERT_TRUE(outline.points()[3] == ZGYAccess::Point2d(3775.0, 2890.0));

    reader.Close();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(reader_tests, testZRange)
{
    ZGYAccess::ZGYReader reader;

    ASSERT_TRUE(reader.Open(std::string(TEST_DATA_DIR) + "Fancy-int8.zgy"));

    auto [zmin, zmax] = reader.ZRange();
    double ztep = reader.ZStep();

    ASSERT_DOUBLE_EQ(ztep, 4.125);
    ASSERT_DOUBLE_EQ(zmin, 2500);
    ASSERT_DOUBLE_EQ(zmax, 3292);


    reader.Close();
}

