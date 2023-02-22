
#pragma once

namespace ZGYAccess
{

    class Point2d
    {
    public:
        Point2d(double x, double y);
        ~Point2d();

        double x() const;
        double y() const;

        double distanceTo(Point2d& b) const;

        bool operator ==(const Point2d& b) const;

    private:
        double m_x;
        double m_y;
    };

}