

#include "zgyreader.h"
#include "openzgy.h"

namespace ZGYAccess
{

ZGYReader::ZGYReader()
{

}

ZGYReader::~ZGYReader()
{

}

bool ZGYReader::Open(std::string filename)
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


void ZGYReader::Close()
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

std::vector<std::pair<std::string, std::string>> ZGYReader::MetaData()
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

    retValues.push_back(std::make_pair("Depth unit", m_reader->zunitname()));
    retValues.push_back(std::make_pair("Depth offset", std::to_string(m_reader->zstart())));
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

std::string ZGYReader::cornerToString(std::array<double, 2> corner)
{
    std::string retval;
    
    retval = "(" + std::to_string(corner[0]) + ", ";
    retval += std::to_string(corner[1]);
    retval += ")";

    return retval;
}

std::string ZGYReader::sizeToString(std::array<std::int64_t, 3> size)
{
    std::string retval;

    retval += std::to_string(size[0]) + "x";
    retval += std::to_string(size[1]) + "x";
    retval += std::to_string(size[2]);

    return retval;
}

std::array<int, 2> ZGYReader::Origin()
{
    if (m_reader == nullptr) return { 0, 0 };

    const auto& annotstart = m_reader->annotstart();
    return { int(annotstart[0]), int(annotstart[1]) };
}

std::array<int, 2> ZGYReader::Size()
{
    if (m_reader == nullptr) return { 0, 0 };

    const auto& annotsize = m_reader->size();
    return { int(annotsize[0]), int(annotsize[1]) };
}

std::array<int, 2> ZGYReader::Step()
{
    if (m_reader == nullptr) return { 0, 0 };

    const auto& annotinc = m_reader->annotinc();
    return { int(annotinc[0]), int(annotinc[1]) };
}





}