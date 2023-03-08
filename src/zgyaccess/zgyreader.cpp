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

namespace ZGYAccess
{

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
ZGYReader::ZGYReader()
{

}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
ZGYReader::~ZGYReader()
{

}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
bool ZGYReader::open(std::string filename)
{
    if (m_reader != nullptr) return false;

    try
    {
        m_reader = OpenZGY::IZgyReader::open(filename);
    }
    catch (const std::exception&)
    {
        m_reader = nullptr;
        return false;
    }

    return true;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
void ZGYReader::close()
{
    if (m_reader == nullptr) return;

    try
    {
        m_reader->close();
    }
    catch (const std::exception&)
    {
    }

    m_reader = nullptr;

    return;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::vector<std::pair<std::string, std::string>> ZGYReader::metaData()
{
    std::vector<std::pair<std::string, std::string>> retValues;

    if (m_reader == nullptr) return retValues;

    auto stats = m_reader->filestats();

    retValues.push_back(std::make_pair("ZGY version", std::to_string(stats->fileVersion())));
    if (stats->isCompressed())
    {
        retValues.push_back(std::make_pair("Compressed", "yes"));
        retValues.push_back(std::make_pair("Compression factor", std::to_string(stats->compressionFactor())));
    }
    retValues.push_back(
        std::make_pair("File size", std::to_string(1.0 * stats->fileSize() / 1024.0 / 1024.0) + " MBytes"));
    retValues.push_back(std::make_pair("Header size", std::to_string(stats->headerSize()) + " Bytes"));

    std::string tmp;

    for (auto& a : m_reader->brickcount())
    {
        tmp.append(sizeToString(a) + " ");
    }
    retValues.push_back(std::make_pair("Brick count", tmp));
    retValues.push_back(std::make_pair("Levels of details", std::to_string(m_reader->nlods())));

    const auto& bricksize = m_reader->bricksize();
    retValues.push_back(std::make_pair("Brick size", sizeToString(bricksize)));

    switch (m_reader->datatype())
    {
    default:
    case OpenZGY::SampleDataType::unknown:
        tmp = "Unknown";
        break;
    case OpenZGY::SampleDataType::int8:
        tmp = "Signed 8-bit";
        break;
    case OpenZGY::SampleDataType::int16:
        tmp = "Signed 16-bit";
        break;
    case OpenZGY::SampleDataType::float32:
        tmp = " 32-bit Floating point";
        break;
    }
    retValues.push_back(std::make_pair("Native data type", tmp));

    const auto& datarange = m_reader->datarange();
    retValues.push_back(std::make_pair("Data range", std::to_string(datarange[0]) + " to " + std::to_string(datarange[1])));

    const auto [zmin, zmax] = zRange();
    retValues.push_back(std::make_pair("Depth unit", m_reader->zunitname()));
    retValues.push_back(std::make_pair("Depth range", std::to_string(zmin) + " to " + std::to_string(zmax)));
    retValues.push_back(std::make_pair("Depth increment", std::to_string(m_reader->zinc())));

    retValues.push_back(std::make_pair("Horizontal unit", m_reader->hunitname()));

    const auto& annotstart = m_reader->annotstart();
    retValues.push_back(std::make_pair("First inline", std::to_string(annotstart[0])));
    retValues.push_back(std::make_pair("First crossline", std::to_string(annotstart[1])));

    const auto& annotinc = m_reader->annotinc();
    retValues.push_back(std::make_pair("Inline increment", std::to_string(annotinc[0])));
    retValues.push_back(std::make_pair("Crossline increment", std::to_string(annotinc[1])));

    const auto& annotsize = m_reader->size();

    retValues.push_back(std::make_pair("Inline size", std::to_string(annotsize[0])));
    retValues.push_back(std::make_pair("Crossline size", std::to_string(annotsize[1])));

    tmp = "";
    for (auto& c : m_reader->corners())
    {
        tmp += cornerToString(c) + " ";
    }
    retValues.push_back(std::make_pair("World coord. corners", tmp));

    tmp = "";
    for (auto& c : m_reader->indexcorners())
    {
        tmp += cornerToString(c) + " ";
    }
    retValues.push_back(std::make_pair("Brick index corners", tmp));

    tmp = "";
    for (auto& c : m_reader->annotcorners())
    {
        tmp += cornerToString(c) + " ";
    }
    retValues.push_back(std::make_pair("Inline/crossline corners", tmp));

    return retValues;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::string ZGYReader::cornerToString(std::array<double, 2> corner)
{
    std::string retval;
    
    retval = "(" + std::to_string(corner[0]) + ", ";
    retval += std::to_string(corner[1]);
    retval += ")";

    return retval;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::string ZGYReader::sizeToString(std::array<std::int64_t, 3> size)
{
    std::string retval;

    retval += std::to_string(size[0]) + "x";
    retval += std::to_string(size[1]) + "x";
    retval += std::to_string(size[2]);

    return retval;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::pair<int, int> ZGYReader::inlineRange() const
{
    if (m_reader == nullptr) return { 0, 0 };

    const auto& annotstart = m_reader->annotstart();
    const auto& annotsize = m_reader->size();
    const auto& annotinc = m_reader->annotinc();

    double stop = annotstart[0] + annotsize[0] * annotinc[0];

    return std::make_pair((int)annotstart[0], (int)stop);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int ZGYReader::inlineStep() const
{
    if (m_reader == nullptr) return 0;

    const auto& annotinc = m_reader->annotinc();
    return (int) annotinc[0];
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int ZGYReader::inlineSize() const
{
    if (m_reader == nullptr) return 0;

    const auto& annotsize = m_reader->size();
    return (int)annotsize[0];
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::pair<int, int> ZGYReader::xlineRange() const
{
    if (m_reader == nullptr) return { 0, 0 };

    const auto& annotstart = m_reader->annotstart();
    const auto& annotsize = m_reader->size();
    const auto& annotinc = m_reader->annotinc();

    double stop = annotstart[1] + annotsize[1] * annotinc[1];

    return std::make_pair((int)annotstart[1], (int)stop);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int ZGYReader::xlineStep() const
{
    if (m_reader == nullptr) return 0;

    const auto& annotinc = m_reader->annotinc();
    return (int)annotinc[1];
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int ZGYReader::xlineSize() const
{
    if (m_reader == nullptr) return 0;

    const auto& annotsize = m_reader->size();
    return (int)annotsize[1];
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::pair<double, double> ZGYReader::zRange() const
{
    if (m_reader == nullptr) return { 0.0, 0.0 };

    double zmin = m_reader->zstart();
    double zmax = zmin;

    zmax += m_reader->zinc() * zSize();
    
    return std::make_pair(zmin, zmax);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
double ZGYReader::zStep() const
{
    if (m_reader == nullptr) return 0.0;

    return 1.0 * m_reader->zinc();
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
int ZGYReader::zSize() const
{
    if (m_reader == nullptr) return 0;

    const auto bricksize = m_reader->bricksize()[2];
    const auto brickcount = m_reader->brickcount()[0][2];
    return bricksize * brickcount;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::pair<double, double> ZGYReader::dataRange() const
{
    if (m_reader == nullptr) return { 0.0, 0.0 };

    auto datarange = m_reader->datarange();

    return std::make_pair(datarange[0], datarange[1]);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
HistogramData* ZGYReader::histogram()
{
    m_histogram.reset();
    if (m_reader == nullptr) return &m_histogram;

    const auto hist = m_reader->histogram();

    int nVals = (int)hist.bins.size();
    if (nVals > 0)
    {
        double binsize = (hist.maxvalue - hist.minvalue) / nVals;
        double offset = binsize / 2;

        for (int i = 0; i < nVals; i++)
        {
            m_histogram.Xvalues.push_back(hist.minvalue + binsize * i + offset);
            m_histogram.Yvalues.push_back(1.0 * hist.bins[i]);
        }
    }
    return &m_histogram;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
Outline ZGYReader::seismicWorldOutline()
{
    Outline retval;

    if (m_reader == nullptr) return retval;

    auto& corners = m_reader->corners();

    for (auto& c : corners)
    {
        retval.addPoint(c[0], c[1]);
    }

    return retval;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::pair<double, double> ZGYReader::toWorldCoordinate(int inLine, int crossLine) const
{
    if (m_reader == nullptr) return { 0, 0 };

    std::array<double, 2> annotCoord{ 1.0 * inLine, 1.0 * crossLine };

    auto worldCoord = m_reader->annotToWorld(annotCoord);

    return std::make_pair(worldCoord[0], worldCoord[1]);
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::shared_ptr<SeismicSliceData> ZGYReader::seismicSlice(std::array<double, 3> worldStart, std::array<double, 3> worldStop)
{
    std::shared_ptr<SeismicSliceData> retval = std::make_shared<SeismicSliceData>(0, 0);

    if (m_reader == nullptr) return retval;



    return retval;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::shared_ptr<SeismicSliceData> ZGYReader::inlineSlice(int inlineIndex)
{
    if (m_reader == nullptr) return std::make_shared<SeismicSliceData>(0, 0);

    int width = xlineSize();
    int depth = zSize();

    std::shared_ptr<SeismicSliceData> retData = std::make_shared<SeismicSliceData>(width, depth);

    OpenZGY::IZgyMeta::size3i_t sliceStart;
    OpenZGY::IZgyMeta::size3i_t sliceSize;

    sliceStart[0] = inlineIndex;
    sliceSize[0] = 1;

    sliceStart[1] = 0;
    sliceSize[1] = width;

    sliceStart[2] = 0;
    sliceSize[2] = depth;

    try
    {
        m_reader->read(sliceStart, sliceSize, retData->values(), 0);
    }
    catch (OpenZGY::Errors::ZgyError &err)
    {
        retData->reset();
    }

    return retData;
}

//--------------------------------------------------------------------------------------------------
///
//--------------------------------------------------------------------------------------------------
std::shared_ptr<SeismicSliceData> ZGYReader::xlineSlice(int xlineIndex)
{
    if (m_reader == nullptr) return std::make_shared<SeismicSliceData>(0, 0);

    int width = inlineSize();
    int depth = zSize();

    std::shared_ptr<SeismicSliceData> retData = std::make_shared<SeismicSliceData>(width, depth);

    OpenZGY::IZgyMeta::size3i_t sliceStart;
    OpenZGY::IZgyMeta::size3i_t sliceSize;

    sliceStart[0] = 0;
    sliceSize[0] = width;

    sliceStart[1] = xlineIndex;
    sliceSize[1] = 1;

    sliceStart[2] = 0;
    sliceSize[2] = depth;

    try
    {
        m_reader->read(sliceStart, sliceSize, retData->values(), 0);
    }
    catch (OpenZGY::Errors::ZgyError &err)
    {
        retData->reset();
    }

    return retData;
}

}