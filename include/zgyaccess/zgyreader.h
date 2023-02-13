

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <memory>

namespace OpenZGY
{
    class IZgyReader;
}

namespace ZGYAccess
{

    class ZGYReader
    {
    public:
        ZGYReader();
        ~ZGYReader();

        bool Open(std::string filename);
        void Close();

        std::vector<std::pair<std::string, std::string>> MetaData();

        std::array<int, 2> Origin();
        std::array<int, 2> Size();
        std::array<int, 2> Step();

        std::vector<std::pair<double, double>> WorldCorners() const;
        std::pair<double, double> ZRange() const;


    private:
        std::string cornerToString(std::array<double, 2> corner);
        std::string sizeToString(std::array<std::int64_t, 3> size);

    private:
        std::string                          m_filename;
        std::shared_ptr<OpenZGY::IZgyReader> m_reader;
    };

}
