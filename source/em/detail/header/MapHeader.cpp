// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <em/detail/header/MapHeader.h>

#include <em/detail/header/data/DummyData.h>
#include <em/detail/header/data/MRCData.h>
#include <em/detail/header/data/RECData.h>
#include <utility/Axis3D.h>
#include <utility/Exceptions.h>

#include <array>
#include <iostream>
#include <unordered_map>

using namespace ausaxs;
using namespace ausaxs::em::detail::header;

namespace {
    // (mapc, mapr, maps) is only meaningful if it is a permutation of {1, 2, 3}; anything else is a corrupt header
    bool is_permutation(const std::array<int, 3>& map_order) noexcept {
        std::array<bool, 3> seen = {false, false, false};
        for (int axis : map_order) {
            if (axis < 1 || 3 < axis) {return false;}
            seen[axis-1] = true;
        }
        return seen[0] && seen[1] && seen[2];
    }
}

template<class T>
std::array<int, 3> em::detail::header::get_axis_order(const T& data) noexcept {
    std::array<int, 3> map_order = {data.mapc, data.mapr, data.maps};

    std::array<int, 3> order = {0, 1, 2};
    if (is_permutation(map_order)) {
        for (int i = 0; i < 3; ++i) {order[map_order[i]-1] = i;}
    }
    return order;
}

template std::array<int, 3> em::detail::header::get_axis_order<MRCData>(const MRCData&) noexcept;
template std::array<int, 3> em::detail::header::get_axis_order<RECData>(const RECData&) noexcept;

template<class T>
Axis3D em::detail::header::make_axes(const T& data) noexcept {
    std::array<int, 3> n = {data.nx, data.ny, data.nz};
    std::array<int, 3> m = {data.mx, data.my, data.mz};
    std::array<double, 3> cella = {data.cella_x, data.cella_y, data.cella_z};

    // the counts are in storage order, so the axis order tells us which of them belongs to each crystallographic axis
    std::array<int, 3> stored_axis = get_axis_order(data);

    std::array<Axis, 3> axes;
    for (int axis = 0; axis < 3; ++axis) {
        int bins = n[stored_axis[axis]];

        // if the sampling rate is unset, assume the stored region spans the full cell
        double width = 0 < m[axis] ? cella[axis]/m[axis] : (0 < bins ? cella[axis]/bins : 0);
        axes[axis] = Axis(0, bins*width, bins);
    }
    return {axes[0], axes[1], axes[2]};
}

template Axis3D em::detail::header::make_axes<MRCData>(const MRCData&) noexcept;
template Axis3D em::detail::header::make_axes<RECData>(const RECData&) noexcept;

int IMapHeader::get_byte_size() const {
    if (!byte_sizes.contains(get_data_type())) {
        throw except::parse_error("MRCHeader::get_byte_size: Unknown data type.");
    };
    return static_cast<int>(byte_sizes.at(get_data_type()));
}

std::ostream& IMapHeader::operator<<(std::ostream& os) const {
    return os << this->to_string();
}

template<class T>
observer_ptr<T> MapHeader<T>::get_data() const noexcept {
    return data.get();}

template<class T>
void MapHeader<T>::set_data(std::unique_ptr<T> data) {this->data = std::move(data);}


template<class T>
MapHeader<T>::MapHeader(std::unique_ptr<T> data) : data(std::move(data)) {}

template<class T>
MapHeader<T>::~MapHeader() = default;

template<class T>
char* MapHeader<T>::get_data_ptr() const {
    return reinterpret_cast<char*>(data.get());
}

template class em::detail::header::MapHeader<MRCData>;
template class em::detail::header::MapHeader<RECData>;
template class em::detail::header::MapHeader<DummyData>;