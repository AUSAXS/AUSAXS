// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinatesTemplate.h>
#include <utility/Concepts.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_width, bool explicit_ff = false>
    struct CompactCoordinatesFF : public CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width, explicit_ff> {
        CompactCoordinatesFF() = default;
        CompactCoordinatesFF(unsigned int size) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width, explicit_ff>(size) {}
        CompactCoordinatesFF(const std::vector<data::AtomFF>& atoms) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width, explicit_ff>(atoms) { validate_ff_types(); }
        CompactCoordinatesFF(const std::vector<data::Body>& bodies) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width, explicit_ff>(bodies) { validate_ff_types(); }
        CompactCoordinatesFF(const std::vector<data::Water>& atoms) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width, explicit_ff>(atoms) { validate_ff_types(); }

        int32_t get_ff_type(int i) const {
            return this->get_non_coordinate_value(i);
        }

        void validate_ff_types() const {
            for (int i = 0; i < static_cast<int>(this->data.size()); ++i) {
                if (get_ff_type(i) == static_cast<unsigned int>(ausaxs::form_factor::form_factor_t::UNKNOWN)) {
                    throw std::runtime_error(
                        "CompactCoordinatesFF: Attempted to use an atom with UNKNOWN form factor type.\n"
                        "Form factor information is required for the selected excluded volume model."
                    );
                }
            }
        }        
    };

    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<true, false>>,  "CompactCoordinatesFF should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<false, false>>, "CompactCoordinatesFF should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<true, true>>,  "CompactCoordinatesFF should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<false, true>>, "CompactCoordinatesFF should support nothrow move semantics.");
}