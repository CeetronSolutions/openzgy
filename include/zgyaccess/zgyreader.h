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

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <cmath>

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
        std::pair<int, int> toInlineXline(double worldX, double worldY) const;

        std::shared_ptr<SeismicSliceData> inlineSlice(int inlineIndex);
        std::shared_ptr<SeismicSliceData> inlineSlice(int inlineIndex, int zStartIndex, int zSize);

        std::shared_ptr<SeismicSliceData> xlineSlice(int xlineIndex);
        std::shared_ptr<SeismicSliceData> xlineSlice(int xlineIndex, int zStartIndex, int zSize);

        std::shared_ptr<SeismicSliceData> zSlice(int zIndex);

        std::shared_ptr<SeismicSliceData> zTrace(int inlineIndex, int xlineIndex);
        std::shared_ptr<SeismicSliceData> zTrace(int inlineIndex, int xlineIndex, int zStartIndex, int zSize);

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
