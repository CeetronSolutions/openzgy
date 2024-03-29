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

#include "zgyaccess/zgyreader.h"

#include "exception.h"
#include "api.h"

#include "zgyaccess/seismicslice.h"

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

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
float SeismicSliceData::valueAt(int width, int depth)
{
    if ((width < m_width) && (depth < m_depth))
    {
        return m_values.get()[width * m_depth + depth];
    }

    return 0.0;
}


//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void SeismicSliceData::limitTo(float minVal, float maxVal)
{
    const int nVals = size();
    float *pData = m_values.get();
    for (int i = 0; i < nVals; i++, pData++)
    {
        const float tmp = *pData;
        if (tmp < minVal)
          *pData = minVal;
        else if (tmp > maxVal)
          *pData = maxVal;
    }
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void SeismicSliceData::mute(float threshold)
{
    const int nVals = size();
    float *pData = m_values.get();
    for (int i = 0; i < nVals; i++, pData++) 
    {
        if (std::abs(*pData) < threshold)
          *pData = 0.0;
    }
}


}
