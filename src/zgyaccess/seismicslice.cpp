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

#include "zgyreader.h"

#include "exception.h"
#include "api.h"

#include "seismicslice.h"

namespace ZGYAccess
{

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
SeismicSliceData::SeismicSliceData(int width, int height)
    : m_width(width)
    , m_depth(height)
{
    const int size = height * width;

    float* buffer = new float[size];
    m_values = std::unique_ptr<float>(buffer);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
SeismicSliceData::~SeismicSliceData()
{

}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int SeismicSliceData::size() const
{
    return m_width * m_depth;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int SeismicSliceData::width() const
{
    return m_width;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int SeismicSliceData::depth() const
{
    return m_depth;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float* SeismicSliceData::values()
{
    return m_values.get();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void SeismicSliceData::reset()
{
    m_width = 0;
    m_depth = 0;

    m_values.reset();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool SeismicSliceData::isEmpty() const
{
    return size() == 0;
}

}