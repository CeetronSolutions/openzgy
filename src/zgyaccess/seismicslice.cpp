

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

}

SeismicSliceData::~SeismicSliceData()
{

}

}