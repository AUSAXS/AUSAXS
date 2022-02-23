#pragma once

#include <array>

using std::array;

namespace em {
    namespace ccp4 {
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
            array<char, 8> extra1;
            array<char, 4> exttyp;
            int nversion;
            array<char, 84> extra2;
            float origin_x;
            float origin_y;
            float origin_z;
            array<char, 4> map;
            unsigned int machst;
            float rms;
            int nlabl;
            array<char, 800> label;
        };
        static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");
    }
}