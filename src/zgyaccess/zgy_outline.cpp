
#include "zgy_outline.h"


namespace ZGYAccess
{

    Outline::Outline()
    {

    }

    Outline::~Outline()
    {

    }

    void Outline::addPoint(Point2d p)
    {
        m_points.push_back(p);
    }

    void Outline::addPoint(double x, double y)
    {
        m_points.push_back(Point2d(x, y));
    }

    std::vector<Point2d> Outline::points() const
    {
        return m_points;
    }

    void Outline::reset()
    {
        m_points.clear();
    }

    bool Outline::isEmpty() const
    {
        return (m_points.size() == 0);
    }

    bool Outline::isValid() const
    {
        return (m_points.size() > 2);
    }


}