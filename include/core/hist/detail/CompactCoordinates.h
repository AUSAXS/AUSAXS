#pragma once

#include <hist/detail/CompactCoordinatesTemplate.h>
#include <utility/Concepts.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_width>
    struct CompactCoordinates : public CompactCoordinatesTemplate<CoordinateTypeXYZW, variable_bin_width> {
        using CompactCoordinatesTemplate<CoordinateTypeXYZW, variable_bin_width>::CompactCoordinatesTemplate;

        float get_weight(unsigned int i) const {
            return this->get_non_coordinate_value(i);
        }
    };

    static_assert(supports_nothrow_move_v<CompactCoordinates<true>>,    "CompactCoordinates should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinates<false>>,   "CompactCoordinates should support nothrow move semantics.");
}