
#pragma once

#include <vector>
#include <memory>

namespace ZGYAccess
{


class SeismicSliceData
{

public:
    SeismicSliceData(int width, int depth);
    ~SeismicSliceData();

    //void addSubSliceData(int fromX, int fromY, int fromWidth, const std::vector<float> values, int toX, int toY);

    float* values();
    int size() const;

    void reset();

private:
    int m_width;
    int m_depth;
    std::unique_ptr<float> m_values;
};




}