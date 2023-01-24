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

            //! members CANNOT be reordered!
            int nx;             // Number of points along x-axis.
            int ny;             // Number of points along y-axis.
            int nz;             // Number of points along z-axis.
            int mode;           // Bit format of the data.
            int nxstart;        // Location of first column in the unit cell.
            int nystart;        // Location of first row in the unit cell.
            int nzstart;        // Location of first section in the unit cell. 
            int mx;             // Sampling rate along the x-axis in the unit cell.
            int my;             // Sampling rate along the y-axis in the unit cell.
            int mz;             // Sampling rate along the z-axis in the unit cell.
            float cella_x;      // Map size in Ångström along x-axis. 
            float cella_y;      // Map size in Ångström along y-axis.
            float cella_z;      // Map size in Ångström along z-axis.
            float cellb_alpha;  // Cell angle alpha in degrees.
            float cellb_beta;   // Cell angle beta in degrees.
            float cellb_gamma;  // Cell angle gamma in degrees.
            int mapc;           // Axis corresponding to columns (1,2,3 for X,Y,Z).
            int mapr;           // Axis corresponding to rows (1,2,3 for X,Y,Z).
            int maps;           // Axis corresponding to sections (1,2,3 for X,Y,Z).
            float dmin;         // Minimum pixel value.
            float dmax;         // Maximum pixel value.
            float dmean;        // Mean pixel value.
            int ispg;           // Space group number.
            int nsymbt;         // Size of extended header in bytes.
            std::array<char, 8> extra1;
            std::array<char, 4> exttyp;
            int nversion;       // Version number of the format.
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
        static_assert(sizeof(Header) == 1024, "em::ccp4::Header: Size of Header is wrong.");
    }
}