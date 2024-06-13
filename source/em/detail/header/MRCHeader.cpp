/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/header/MRCHeader.h>
#include <utility/Exceptions.h>
#include <utility/Axis3D.h>
#include <io/ExistingFile.h>

#include <iostream>
#include <sstream>

using namespace em::detail::header;

MRCHeader::MRCHeader() : MapHeader(std::make_unique<MRCData>()) {}
MRCHeader::MRCHeader(MRCData&& data) : MapHeader(std::make_unique<MRCData>(std::move(data))) {}
MRCHeader::~MRCHeader() = default;

bool MRCHeader::is_mrc(const io::ExistingFile& file) {
    auto extension = file.extension();
    if (extension == ".mrc") {return true;}
    if (extension == ".ccp4") {return true;}
    if (extension == ".map") {return true;}
    return false;
}

std::string MRCHeader::to_string() const {
    auto& p = cast_data();

    std::stringstream s;
    s << "HEADER CONTENTS: \n";
    s << "Type: MRC\n";
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

unsigned int MRCHeader::get_header_size() const {
    return sizeof(MRCData);
}

Axis3D MRCHeader::get_axes() const noexcept {
    auto& p = cast_data();
    return Axis3D(
        Axis(0, p.cella_x, p.nx),
        Axis(0, p.cella_y, p.ny),
        Axis(0, p.cella_z, p.nz)
    );
}

std::tuple<unsigned int, unsigned int, unsigned int> MRCHeader::get_axis_order() const noexcept {
    auto& p = cast_data();
    return std::make_tuple(p.mapc, p.mapr, p.maps);
}

em::detail::header::DataType MRCHeader::get_data_type() const {
    auto& p = cast_data();

    switch(p.mode) {
        case 0:
            return em::detail::header::DataType::int8;

        case 1:
            return em::detail::header::DataType::int16;

        case 2:
            return em::detail::header::DataType::float32;

        case 6:
            return em::detail::header::DataType::uint16;

        case 12:
            return em::detail::header::DataType::float16;

        case 3:
            throw except::parse_error("MRCHeader::get_data_type: Complex data format is not currently supported.");

        case 4:
            throw except::parse_error("MRCHeader::get_data_type: Complex data format is not currently supported.");

        default:
            throw except::parse_error("MRCHeader::get_data_type: Unrecognized data format (mode " + std::to_string(p.mode) + ").");
    }
}

MRCData& MRCHeader::cast_data() const noexcept {
    return *static_cast<MRCData*>(MapHeader::get_data());
}

void MRCHeader::rotate(int x, int y, int z) noexcept {
    auto& p = cast_data();

    p.mapc = x;
    p.mapr = y;
    p.maps = z;
}
