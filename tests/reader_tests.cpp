

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
