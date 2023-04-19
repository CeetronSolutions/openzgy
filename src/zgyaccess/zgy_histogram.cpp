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

#include "zgy_histogram.h"
#include "exception.h"
#include "api.h"
#include "histogramdata.h"
#include "histogrambuilder.h"

namespace ZGYAccess
{

std::unique_ptr<HistogramData> HistogramGenerator::getHistogram(std::vector<float> values, int nBins, float minVal, float maxVal)
{
    auto retData = std::make_unique<HistogramData>();

    InternalZGY::HistogramBuilder builder(nBins, minVal, maxVal);
    builder.add(values.begin(), values.end());

    auto& tmpHist = builder.gethisto();

    int nVals = (int)tmpHist.getsize();
    if (nVals > 0)
    {
        auto bins = tmpHist.getbins();

        double binsize = (maxVal - minVal) / nVals;
        double offset = binsize / 2;

        for (int i = 0; i < nVals; i++)
        {
            retData->Xvalues.push_back(minVal + binsize * i + offset);
            retData->Yvalues.push_back(1.0 * bins[i]);
        }
    }

    return retData;
}


}