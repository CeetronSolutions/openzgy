
#pragma once

#include <vector>


namespace ZGYAccess
{


class SeismicSliceData
{

public:
    SeismicSliceData(int width, int depth);
    ~SeismicSliceData();

    void addSubSliceData(int fromX, int fromY, int fromWidth, const std::vector<float> values, int toX, int toY);

    const std::vector<float> values();

private:
    int m_width;
    int m_depth;
    std::vector<float> m_values;
};




}