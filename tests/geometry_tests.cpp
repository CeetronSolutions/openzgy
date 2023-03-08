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

#include "zgyaccess/zgy_point.h"
#include "zgyaccess/zgy_line.h"
#include "zgyaccess/zgy_outline.h"

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(geometry_tests, testPoint)
{
    ZGYAccess::Point2d p1(1.0, 2.0);
    ZGYAccess::Point2d p2(3.0, 4.0);

    ASSERT_DOUBLE_EQ(p1.x(), 1.0);
    ASSERT_DOUBLE_EQ(p1.y(), 2.0);
    ASSERT_DOUBLE_EQ(p2.x(), 3.0);
    ASSERT_DOUBLE_EQ(p2.y(), 4.0);

    ASSERT_FALSE(p1 == p2);
    ASSERT_TRUE(p2 == p2);

    ASSERT_DOUBLE_EQ(p1.distanceTo(p2), std::sqrt(8.0));
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(geometry_tests, testLine)
{
    ZGYAccess::LineSegment l1(1.0, 2.0, 3.0, 4.0);

    ZGYAccess::Point2d p1(11.0, 12.0);
    ZGYAccess::Point2d p2(13.0, 14.0);

    ZGYAccess::LineSegment l2(p1, p2);

    ASSERT_TRUE(l2.p1() == p1);
    ASSERT_TRUE(l2.p2() == p2);



}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
TEST(geometry_tests, testOutline)
{

    ZGYAccess::Point2d p1(1.0, 10.0);
    ZGYAccess::Point2d p2(13.0, 11.0);
    ZGYAccess::Point2d p3(14.0, 2.0);

    ZGYAccess::Outline o;

    ASSERT_TRUE(o.isEmpty());

    o.addPoint(p1);
    ASSERT_FALSE(o.isEmpty());
    ASSERT_FALSE(o.isValid());

    o.addPoint(p2);
    ASSERT_FALSE(o.isEmpty());
    ASSERT_FALSE(o.isValid());

    o.addPoint(p3);
    ASSERT_FALSE(o.isEmpty());
    ASSERT_TRUE(o.isValid());

    o.addPoint(2.0, 1.0);
    ASSERT_FALSE(o.isEmpty());
    ASSERT_TRUE(o.isValid());

}
