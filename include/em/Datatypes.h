#pragma once

#include <array>

namespace em {
    namespace ccp4 {
        /**
         * @brief The 1024-byte header of CCP4 files.
         *        It assumes the endian of the data is the same as the system. 
         */
        struct Header {
            Header() {}

            // members CANNOT be reordered!
            int nx; // # of cols
            int ny; // # of rows
            int nz; // # of sections
            int mode; // bit format of data section
            int nxstart; // start position of image
            int nystart;
            int nzstart;
            int mx; // grid size
            int my;
            int mz;
            float cella_x; // cell size
            float cella_y;
            float cella_z;
            float cellb_alpha; // cell angles
            float cellb_beta;
            float cellb_gamma;
            int mapc;
            int mapr;
            int maps;
            float dmin; // minimum pixel value
            float dmax; // maximum pixel value
            float dmean; // mean pixel value
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

            std::string to_string() const;

            friend std::ostream& operator<<(std::ostream& os, const Header& h) {os << h.to_string(); return os;}

            size_t get_byte_size() const;
        };
        static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");
    }
}