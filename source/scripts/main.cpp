#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <iostream>

using std::string, std::cout, std::endl;

namespace CCP4 {
    struct Header {
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

int main(int argc, char const *argv[]) {
    std::cout << "Struct has length: " << sizeof(CCP4::Header) << std::endl;
    return 0;
}