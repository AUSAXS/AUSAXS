#include <em/Datatypes.h>
#include <utility/Exceptions.h>

#include <iostream>
#include <sstream>

std::string em::ccp4::Header::to_string() const {
    std::stringstream s;
    s << "HEADER CONTENTS: \n";
    s << "Dimensions: (" << nx << ", " << ny << ", " << nz << ")\n";
    s << "Offsets: (" << nxstart << ", " << nystart << ", " << nzstart << ")\n";
    s << "Grid size: (" << mx << ", " << my << ", " << mz << ")\n";
    s << "Cell size in Å: (" << cella_x << ", " << cella_y << ", " << cella_z << ")\n";
    s << "Cell angles in degrees: (" << cellb_alpha << ", " << cellb_beta << ", " << cellb_gamma << ")\n";
    s << "Order of dimensions: " << mapc << ", " << mapr << ", " << maps << "\n";
    s << "Minimum pixel: " << dmin << ", maximum pixel: " << dmax << ", mean: " << dmean << "\n";
    s << "Space group number: " << ispg << "\n";
    s << "Extended header size: " << nsymbt << "\n";
    return s.str();
}

size_t em::ccp4::Header::get_byte_size() const {
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