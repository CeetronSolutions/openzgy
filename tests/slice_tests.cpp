

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

    ASSERT_NE(slice.values(), nullptr);

    slice.reset();
    ASSERT_EQ(slice.size(), 0);
    ASSERT_EQ(slice.width(), 0);
    ASSERT_EQ(slice.depth(), 0);
    ASSERT_EQ(slice.values(), nullptr);

}
