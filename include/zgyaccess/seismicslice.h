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

#include <vector>
#include <memory>

namespace ZGYAccess
{


class SeismicSliceData
{

public:
    SeismicSliceData(int width, int depth);
    ~SeismicSliceData();

    //void addSubSliceData(int fromX, int fromY, int fromWidth, const std::vector<float> values, int toX, int toY);

    float* values();
    int size() const;
    int depth() const;
    int width() const;

    void transpose();

    void reset();

private:
    int m_width;
    int m_depth;
    std::unique_ptr<float> m_values;
};




}