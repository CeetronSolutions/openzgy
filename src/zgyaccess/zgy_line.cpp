
#include "zgy_line.h"

#include <cmath>

namespace ZGYAccess
{

    LineSegment::LineSegment(Point2d p1, Point2d p2)
        : m_p1(p1)
        , m_p2(p2)
    {

    }

    LineSegment::LineSegment(double x1, double y1, double x2, double y2)
        : m_p1(x1,y1)
        , m_p2(x2,y2)
    {

    }


    LineSegment::~LineSegment()
    {
    }

    Point2d LineSegment::p1() const
    {
        return m_p1;
    }

    Point2d LineSegment::p2() const
    {
        return m_p2;
    }

    double LineSegment::distanceTo(Point2d& p) const
    {
        bool onLine = false;

        Point2d p2 = closestPositionTo(p, onLine);

        if (!onLine) return std::numeric_limits<double>::infinity();

        return p.distanceTo(p2);
    }

    Point2d LineSegment::closestPositionTo(Point2d& p, bool& onLine) const
    {
        double dx1 = p.x() - m_p1.x();
        double dx2 = m_p2.x() - m_p1.x();
        double dy1 = p.x() - m_p1.y();
        double dy2 = m_p2.y() - m_p1.y();

        double reldist = (dx1 * dx2 + dy1 * dy2) / (dx2 * dx2 + dy2 * dy2);

        onLine = ((reldist >= 0.0) && (reldist <= 1.0));

        return Point2d(m_p1.x() + dx2 * reldist, m_p1.y() + dy2 * reldist);
    }


}