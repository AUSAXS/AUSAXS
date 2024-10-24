/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/header/MapHeader.h>
#include <em/detail/header/data/MRCData.h>
#include <em/detail/header/data/RECData.h>
#include <em/detail/header/data/DummyData.h>
#include <utility/Exceptions.h>

#include <unordered_map>
#include <iostream>

using namespace ausaxs::em::detail::header;

unsigned int IMapHeader::get_byte_size() const {
    if (byte_sizes.contains(get_data_type()) == false) {
        throw except::parse_error("MRCHeader::get_byte_size: Unknown data type.");
    };
    return byte_sizes.at(get_data_type());
}

std::ostream& IMapHeader::operator<<(std::ostream& os) {
    return os << this->to_string();
}

template<class T>
observer_ptr<T> MapHeader<T>::get_data() const noexcept {
    return data.get();}

template<class T>
void MapHeader<T>::set_data(std::unique_ptr<T> data) {this->data = std::move(data);}


template<class T>
MapHeader<T>::MapHeader(std::unique_ptr<T> data) {
    this->data = std::move(data);
}

template<class T>
MapHeader<T>::~MapHeader() = default;

template<class T>
char* MapHeader<T>::get_data_ptr() const {
    return reinterpret_cast<char*>(data.get());
}

template class em::detail::header::MapHeader<MRCData>;
template class em::detail::header::MapHeader<RECData>;
template class em::detail::header::MapHeader<DummyData>;