

#include "zgyreader.h"

#include "exception.h"
#include "api.h"

#include "seismicslice.h"

namespace ZGYAccess
{

SeismicSliceData::SeismicSliceData(int width, int height)
    : m_width(width)
    , m_depth(height)
{
    const int size = height * width;

    float* buffer = new float[size];
    m_values = std::unique_ptr<float>(buffer);
}

SeismicSliceData::~SeismicSliceData()
{

}

int SeismicSliceData::size() const
{
    return m_width * m_depth;
}

float* SeismicSliceData::values()
{
    return m_values.get();
}

void SeismicSliceData::reset()
{
    m_width = 0;
    m_depth = 0;

    m_values.reset();
}


}