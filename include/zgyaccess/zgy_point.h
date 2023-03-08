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