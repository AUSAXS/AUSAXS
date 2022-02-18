#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <fstream>
#include <iostream>
#include <memory>

using std::string, std::cout, std::endl;

namespace CCP4 {
    /**
     * @brief The 1024-byte header of CCP4 files.
     *        It assumes the endian of the data is the same as the system. 
     */
    struct Header {
        // members CANNOT be reordered!
        int nx;
        int ny;
        int nz;
        int mode;
        int nxstart;
        int nystart;
        int nzstart;
        int mx;
        int my;
        int mz;
        float cella_x;
        float cella_y;
        float cella_z;
        float cellb_alpha;
        float cellb_beta;
        float cellb_gamma;
        int mapc;
        int mapr;
        int maps;
        float dmin;
        float dmax;
        float dmean;
        int ispg;
        int nsymbt;
        std::array<char, 8> extra1;
        std::array<char, 4> exttyp;
        int nversion;
        std::array<char, 84> extra2;
        float origin_x;
        float origin_y;
        float origin_z;
        std::array<char, 4> map;
        unsigned int machst;
        float rms;
        int nlabl;
        std::array<char, 800> label;
    };
    static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");

    
}

using std::vector;

enum class Type{int8, int16, uint16, float16, float32, complex64};

template<typename T>
vector<T> create_vector(Type type) {
    switch type {
        case Type::int8:
            return vector<short int>();
        case Type::int16:
            return vector<int>();
        case Type::uint16:
            return vector<unsigned int>();
        case Type::float16:
            return vector<short float>();
        case Type::float32:
            return vector<float>();
        default:
            //! throw error
    }
}

template<typename T>
class Image {
    public:
        Image(const CCP4::Header& header, std::ifstream& istream) : header(header) {
            read(istream);
        }

    private:
        const CCP4::Header header;
        vector<T> data;

        read(std::ifstream& istream) {
            data = vector<vector<vector<T>>>(header.nx, vector<vector<T>>(header.ny, vector<T>(header.nz)));

        }
};

int main(int argc, char const *argv[]) {
    std::ifstream input("data/A2M_map.ccp4", std::ios::binary);
    CCP4::Header header;
    input.read(reinterpret_cast<char*>(&header), sizeof(header));

    Type type;
    switch(header.mode) {
        case 0: { // int8 --> short int
            Image<short int> image(header, input);
        }
        case 1: { // int16 --> int
            Image<int> image(header, input);
        }
        case 2: { // float32 --> float
            Image<float> image(header, input);
        }
        case 6: { // uint16 --> unsigned int
            Image<unsigned int> image(header, input);
        }
        case 12: { // float16 --> short float
            // Image<short float> image(header, input);
            //! short float doesn't exist in pure c++
        }
        case 3:
        case 4:
            type = Type::complex64;
            break;
        default:
            //! throw error
    }
    vector<vector<vector<double>>> data(header.nx, vector<vector<double>>(header.ny, vector<double>(header.nz, 0)));

    std::cout << header.nx << ", " << header.ny << ", " << header.nz << ", " << header.mode << std::endl;
    std::cout << header.cella_x << ", " << header.cella_y << ", " << header.cella_z << std::endl;
    std::cout << header.nlabl << std::endl;
    std::cout << string(header.label.data()).substr(0, 80) << std::endl;
    return 0;
}