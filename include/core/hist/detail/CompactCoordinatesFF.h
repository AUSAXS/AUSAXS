#pragma once

#include <hist/detail/CompactCoordinatesTemplate.h>
#include <utility/Concepts.h>

namespace ausaxs::hist::detail {
    template<bool variable_bin_width>
    struct CompactCoordinatesFF : public CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width> {
        using CompactCoordinatesTemplate<CoordinateTypeXYZFF, variable_bin_width>::CompactCoordinatesTemplate;

        int32_t get_ff_type(unsigned int i) const {
            return this->get_non_coordinate_value(i);
        }

        void init() const override {
            for (unsigned int i = 0; i < this->data.size(); ++i) {
                unsigned int ff_type = this->data.ff_types[i];
                if (ff_type == static_cast<unsigned int>(ausaxs::form_factor::form_factor_t::UNKNOWN)) {
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