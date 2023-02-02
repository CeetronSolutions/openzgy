

#include <vector>


namespace ZGYAccess
{


class SeismicSliceData
{

public:
    SeismicSliceData();
    ~SeismicSliceData();


private:
    int m_nTraces;
    int m_nSamplesInDepth;


    std::vector<float> m_values;

};




}