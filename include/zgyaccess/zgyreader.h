
#include <string>
#include <vector>
#include <array>
#include <utility>
#include <memory>

#include "seismicslice.h"
#include "zgy_outline.h"

namespace OpenZGY
{
    class IZgyReader;
}

namespace ZGYAccess
{

    class HistogramData
    {
    public: 
        HistogramData() {};

        void reset()
        {
            Xvalues.clear();
            Yvalues.clear();
        };

        std::vector<double> Xvalues;
        std::vector<double> Yvalues;
    };

    class ZGYReader
    {
    public:
        ZGYReader();
        ~ZGYReader();

        bool open(std::string filename);
        void close();

        std::vector<std::pair<std::string, std::string>> metaData();

        std::pair<double, double> zRange() const;
        double zStep() const;
        int zSize() const;

        std::pair<int, int> inlineRange() const;
        int inlineStep() const;
        int inlineSize() const;

        std::pair<int, int> xlineRange() const;
        int xlineStep() const;
        int xlineSize() const;

        std::pair<double, double> dataRange() const;

        std::pair<double, double> toWorldCoordinate(int inLine, int crossLine) const;

        std::shared_ptr<SeismicSliceData> seismicSlice(std::array<double, 3> worldStart, std::array<double, 3> worldStop);

        std::shared_ptr<SeismicSliceData> inlineSlice(int inlineIndex);
        std::shared_ptr<SeismicSliceData> xlineSlice(int xlineIndex);

        HistogramData* histogram();

        Outline seismicWorldOutline();


    private:
        std::string cornerToString(std::array<double, 2> corner);
        std::string sizeToString(std::array<std::int64_t, 3> size);

    private:
        std::string                          m_filename;
        std::shared_ptr<OpenZGY::IZgyReader> m_reader;

        HistogramData m_histogram;
    };

}
