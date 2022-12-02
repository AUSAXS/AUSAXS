#pragma once

#include <array>
#include <string>
#include <ostream>
#include <cstddef>

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

            /**
             * @brief Create a string representation of this object.
             */
            std::string to_string() const noexcept;

            /**
             * @brief Output this object to a stream.
             */
            friend std::ostream& operator<<(std::ostream& os, const Header& h) {os << h.to_string(); return os;}

            /**
             * @brief Get the byte size of each voxel.
             */
            unsigned int get_byte_size() const;

            /**
             * @brief Rotate the map contents. This does not affect the operation of this program.
             *        The arguments must be some permutation of {1, 2, 3}.
             */
            void rotate(int x, int y, int z) noexcept;
        };
        static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");
    }
}