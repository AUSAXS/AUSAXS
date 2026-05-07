// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinatesTemplate.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/Concepts.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_width>
    struct CompactCoordinatesFF : public CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width> {
        CompactCoordinatesFF() = default;
        CompactCoordinatesFF(const std::vector<data::AtomFF>& atoms) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width>(atoms) {setup();}
        CompactCoordinatesFF(const std::vector<data::Body>& bodies) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width>(bodies) {setup();}
        CompactCoordinatesFF(const std::vector<data::Water>& atoms) : CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width>(atoms) {setup();}

        int32_t get_ff_type(int i) const {return this->get_non_coordinate_value(i);}
        int32_t& get_ff_type(int i) {return this->get_non_coordinate_value(i);}


        void setup() {
            auto map = form_factor::FormFactorManager::get_ff_mapping();
            for (unsigned int i = 0; i < this->data.size(); ++i) {
                get_ff_type(i) = map[get_ff_type(i)]; // remap to current form factor indices
                if (get_ff_type(i) == static_cast<int>(ausaxs::form_factor::form_factor_t::UNKNOWN)) {
                    throw std::runtime_error(
                        "CompactCoordinatesFF: Attempted to use an atom with UNKNOWN form factor type.\n"
                        "Form factor information is required for the selected excluded volume model."
                    );
                }
            }
        }
    };

    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<true>>,  "CompactCoordinatesFF should support nothrow move semantics.");
    static_assert(supports_nothrow_move_v<CompactCoordinatesFF<false>>, "CompactCoordinatesFF should support nothrow move semantics.");
}