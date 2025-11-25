#pragma once

#include <hist/detail/CompactCoordinatesTemplate.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_width>
    struct CompactCoordinates : public CompactCoordinatesTemplate<CoordinateTypeXYZW, variable_bin_width> {
        using CompactCoordinatesTemplate<CoordinateTypeXYZW, variable_bin_width>::CompactCoordinatesTemplate;

        float get_weight(unsigned int i) const {
            return this->get_non_coordinate_value(i);
        }
    };

    template<bool variable_bin_width>
    struct CompactCoordinatesFF : public CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width> {
        using CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width>::CompactCoordinatesTemplate;

        int32_t get_ff_type(unsigned int i) const {
            return this->get_non_coordinate_value(i);
        }
    };

    static_assert(supports_nothrow_move_v<CompactCoordinates<true>>,    "CompactCoordinates should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinates<false>>,   "CompactCoordinates should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<true>>,  "CompactCoordinatesFF should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<false>>, "CompactCoordinatesFF should support nothrow move semantics.");
}