#pragma once

#include <array>
#include <Exceptions.h>

using std::array;

namespace em {
    namespace ccp4 {
        /**
         * @brief The 1024-byte header of CCP4 files.
         *        It assumes the endian of the data is the same as the system. 
         */
        struct Header {
            Header() {}

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

            size_t get_byte_size() const {
                switch(mode) {
                    case 0: // int8 --> short int
                        return sizeof(short int);

                    case 1: // int16 --> int
                        return sizeof(int);

                    case 2: // float32 --> float
                        return sizeof(float);

                    case 6: // uint16 --> unsigned int
                        return sizeof(unsigned int);

                    case 12: // float16 --> short float
                        throw except::parse_error("Short float data format is not currently supported.");

                    case 3: // complex32 (2x 16bit int)
                        throw except::parse_error("Complex data format is not currently supported.");

                    case 4: // complex64 (2x 32bit float)
                        throw except::parse_error("Complex data format is not currently supported.");

                    default:
                        throw except::parse_error("Unrecognized data format (mode " + std::to_string(mode) + ").");
                }
            }
        };
        static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");
    }
}