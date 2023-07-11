#include <em/detail/header/DummyHeader.h>
#include <em/detail/header/data/DummyData.h>
#include <utility/Axis3D.h>

using namespace em::detail::header;

DummyHeader::DummyHeader() : MapHeader(std::make_unique<DummyData>()) {}
DummyHeader::~DummyHeader() = default;

std::string DummyHeader::to_string() const noexcept {
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