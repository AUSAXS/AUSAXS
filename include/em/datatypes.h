#pragma once

#include <array>
#include <Exceptions.h>
#include <iostream>
#include <sstream>

using std::array, std::endl;

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

            std::string to_string() const {
                std::stringstream s;
                s << "HEADER CONTENTS: " << endl;
                s << "Dimensions: (" << nx << ", " << ny << ", " << nz << ")" << endl;
                s << "Offsets: (" << nxstart << ", " << nystart << ", " << nzstart << ")" << endl;
                s << "Grid size: (" << mx << ", " << my << ", " << mz << ")" << endl;
                s << "Cell size in Å: (" << cella_x << ", " << cella_y << ", " << cella_z << ")" << endl; 
                s << "Cell angles in degrees: (" << cellb_alpha << ", " << cellb_beta << ", " << cellb_gamma << ")" << endl;
                s << "Order of dimensions: " << mapc << ", " << mapr << ", " << maps << endl;
                s << "Minimum pixel: " << dmin << ", maximum pixel: " << dmax << ", mean: " << dmean << endl;
                s << "Space group number: " << ispg << endl;
                s << "Extended header size: " << nsymbt << endl;
                return s.str();
            }

            friend std::ostream& operator<<(std::ostream& os, const Header& h) {os << h.to_string(); return os;}

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
                        throw except::parse_error("Error in Header::get_byte_size: Short float data format is not currently supported.");

                    case 3: // complex32 (2x 16bit int)
                        throw except::parse_error("Error in Header::get_byte_size: Complex data format is not currently supported.");

                    case 4: // complex64 (2x 32bit float)
                        throw except::parse_error("Error in Header::get_byte_size: Complex data format is not currently supported.");

                    default:
                        throw except::parse_error("Error in Header::get_byte_size: Unrecognized data format (mode " + std::to_string(mode) + ").");
                }
            }
        };
        static_assert(sizeof(Header) == 1024, "Size of CCP4 is wrong.");
    }
}