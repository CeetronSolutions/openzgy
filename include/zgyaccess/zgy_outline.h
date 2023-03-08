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