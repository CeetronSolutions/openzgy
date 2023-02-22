
#pragma once

#include "zgy_point.h"

#include <utility>
#include <vector>

namespace ZGYAccess
{

    class Outline
    {
    public:
        Outline();
        ~Outline();

        void addPoint(Point2d p);
        void addPoint(double x, double y);

        std::vector<Point2d> points() const;

        bool isValid() const;
        bool isEmpty() const;
        void reset();

    private:
        std::vector<Point2d> m_points;
    };

}