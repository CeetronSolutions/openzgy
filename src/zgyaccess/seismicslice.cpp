

#include "zgyreader.h"
#include "openzgy.h"

#include "seismicslice.h"

namespace ZGYAccess
{

SeismicSliceData::SeismicSliceData(int width, int height)
    : m_width(width)
    , m_depth(height)
{

}

SeismicSliceData::~SeismicSliceData()
{

}

}