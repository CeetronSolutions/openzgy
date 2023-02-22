
#include "zgy_point.h"

#include <cmath>

namespace ZGYAccess
{

Point2d::Point2d(double x, double y)
    : m_x(x)
    , m_y(y)
{

}

Point2d::~Point2d()
{
}

double Point2d::x() const
{
    return m_x;
}

double Point2d::y() const
{
    return m_y;
}

bool Point2d::operator ==(const Point2d& b) const
{
    return (b.m_x == m_x) && (b.m_y == m_y);
}

double Point2d::distanceTo(Point2d& b) const
{
    double dx = b.m_x - m_x;
    double dy = b.m_y - m_y;

    return std::sqrt(dx * dx + dy * dy);
}



}