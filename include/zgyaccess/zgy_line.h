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

#pragma once

#include "zgy_point.h"

namespace ZGYAccess
{

    class LineSegment
    {
    public:
        LineSegment(double x1, double y1, double x2, double y2);
        LineSegment(Point2d p1, Point2d p2);
        ~LineSegment();

        Point2d p1() const;
        Point2d p2() const;

        //Point2d intersectsWith(Line2d other, bool& noIntersection);

        double distanceTo(Point2d& p) const;

        Point2d closestPositionTo(Point2d& p, bool& onLine) const;


    private:
        Point2d m_p1;
        Point2d m_p2;
    };

}