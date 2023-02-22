
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