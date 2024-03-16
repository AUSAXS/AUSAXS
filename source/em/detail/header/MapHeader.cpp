/*
This software is distributed under the GNU General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/header/MapHeader.h>
#include <utility/Exceptions.h>

#include <unordered_map>
#include <iostream>

using namespace em::detail::header;

MapHeader::MapHeader(std::unique_ptr<HeaderData> data) {
    this->data = std::move(data);
}

MapHeader::~MapHeader() = default;

unsigned int MapHeader::get_byte_size() const {
    if (byte_sizes.contains(get_data_type()) == false) {
        throw except::parse_error("MRCHeader::get_byte_size: Unknown data type.");
    };
    return byte_sizes.at(get_data_type());
}

observer_ptr<HeaderData> MapHeader::get_data() const noexcept {
    return data.get();}

void MapHeader::set_data(std::unique_ptr<HeaderData> data) {this->data = std::move(data);}

std::ostream& em::detail::header::operator<<(std::ostream& os, const MapHeader& h) {
    return os << h.to_string();
}