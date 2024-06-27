/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/detail/header/DummyHeader.h>
#include <em/detail/header/data/DummyData.h>
#include <utility/Axis3D.h>

#include <string>

using namespace em::detail::header;

DummyHeader::DummyHeader() : MapHeader(std::make_unique<DummyData>()) {}
DummyHeader::~DummyHeader() = default;

std::string DummyHeader::to_string() const {
    return "DummyHeader";
}

unsigned int DummyHeader::get_header_size() const {
    return sizeof(DummyData);
}

Axis3D DummyHeader::get_axes() const noexcept {
    return Axis3D(
        Axis(0, 0, 0),
        Axis(0, 0, 0),
        Axis(0, 0, 0)
    );
}

em::detail::header::DataType DummyHeader::get_data_type() const {
    return em::detail::header::DataType::NONE;
}