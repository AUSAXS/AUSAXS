/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/header/RECHeader.h>
#include <utility/Exceptions.h>
#include <utility/Axis3D.h>
#include <io/ExistingFile.h>

#include <iostream>
#include <sstream>

using namespace em::detail::header;

RECHeader::RECHeader() : MapHeader(std::make_unique<RECData>()) {}
RECHeader::~RECHeader() = default;

bool RECHeader::is_rec(const io::ExistingFile& file) {
    const auto& extension = file.extension();
    if (extension == ".rec") {return true;}
    return false;
}

std::string RECHeader::to_string() const {
    auto& p = cast_data();

    std::stringstream s;
    s << "HEADER CONTENTS: \n";
    s << "Type: IMOD REC\n";
    s << "Dimensions: (" << p.nx << ", " << p.ny << ", " << p.nz << ")\n";
    s << "Offsets: (" << p.nxstart << ", " << p.nystart << ", " << p.nzstart << ")\n";
    s << "Grid size: (" << p.mx << ", " << p.my << ", " << p.mz << ")\n";
    s << "Cell size in Å: (" << p.cella_x << ", " << p.cella_y << ", " << p.cella_z << ")\n";
    s << "Cell angles in degrees: (" << p.cellb_alpha << ", " << p.cellb_beta << ", " << p.cellb_gamma << ")\n";
    s << "Order of dimensions: " << p.mapc << ", " << p.mapr << ", " << p.maps << "\n";
    s << "Minimum pixel: " << p.dmin << ", maximum pixel: " << p.dmax << ", mean: " << p.dmean << "\n";
    s << "Space group number: " << p.ispg << "\n";
    s << "Extended header size: " << p.nsymbt << "\n";
    s << "Byte size: " << get_byte_size() << " (mode: " << p.mode << ")\n";
    return s.str();
}

unsigned int RECHeader::get_header_size() const {
    return sizeof(RECData);
}

Axis3D RECHeader::get_axes() const noexcept {
    auto& p = cast_data();
    return Axis3D(
        Axis(0, p.cella_x, p.nx),
        Axis(0, p.cella_y, p.ny),
        Axis(0, p.cella_z, p.nz)
    );
}

std::tuple<unsigned int, unsigned int, unsigned int> RECHeader::get_axis_order() const noexcept {
    auto& p = cast_data();
    return std::make_tuple(p.mapc, p.mapr, p.maps);
}

em::detail::header::DataType RECHeader::get_data_type() const {
    auto& p = cast_data();

    switch(p.mode) {
        case 0:
            if (flags_enabled()) {
                if (flag_four_bit_vals()) {
                    throw except::parse_error("RECHeader::get_data_type: 4-bit values are not currently supported.");
                }
                if (flag_signed_bytes()) {
                    return em::detail::header::DataType::int8;
                } else {
                    return em::detail::header::DataType::uint8;
                }
            }
            return em::detail::header::DataType::uint8;

        case 1:
            return em::detail::header::DataType::int16;

        case 2:
            return em::detail::header::DataType::float32;

        case 6:
            return em::detail::header::DataType::uint16;

        case 12:
            return em::detail::header::DataType::float16;

        case 3:
            throw except::parse_error("RECHeader::get_data_type: Complex data format is not currently supported.");

        case 4:
            throw except::parse_error("RECHeader::get_data_type: Complex data format is not currently supported.");

        default:
            throw except::parse_error("RECHeader::get_data_type: Unrecognized data format (mode " + std::to_string(p.mode) + ").");
    }
}

RECData& RECHeader::cast_data() const noexcept {
    return *static_cast<RECData*>(MapHeader::get_data());
}

bool RECHeader::flags_enabled() const noexcept {
    auto& p = cast_data();
    return p.imodstamp == 1146047817;
}

bool RECHeader::flag_signed_bytes() const noexcept {
    auto& p = cast_data();
    return (p.imodflags & 1);
}

bool RECHeader::flag_four_bit_vals() const noexcept {
    auto& p = cast_data();
    return (p.imodflags & 16);
}
