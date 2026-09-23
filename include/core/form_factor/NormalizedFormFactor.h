// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactor.h>

#include <utility/Exceptions.h>

#include <array>
#include <utility>

namespace ausaxs::form_factor {
    struct NormalizedFormFactor : public FormFactor {
        constexpr NormalizedFormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : FormFactor(a, b, c) {set_normalization(1);}
        constexpr NormalizedFormFactor(const ExvFormFactor& ffx) : FormFactor(ffx) {set_normalization(1);}
        constexpr NormalizedFormFactor(const FormFactor& ff) : FormFactor(ff) {set_normalization(1);}
    };

    /**
     * The normalized vacuum form factors of all form factor types, as described by form_factor::detail::ff_info_table.
     */
    namespace lookup::atomic::normalized {
        namespace detail {
            constexpr auto table = [] <std::size_t... I> (std::index_sequence<I...>) {
                return std::array<NormalizedFormFactor, total_ff_count>{NormalizedFormFactor(raw::detail::table[I])...};
            }(std::make_index_sequence<total_ff_count>{});
        }

        constexpr const NormalizedFormFactor& get(form_factor_t type) {
            if (type == form_factor_t::UNKNOWN) {
                throw ausaxs::except::runtime_error(
                    "form_factor::lookup::atomic::normalized::get: Attempted to get the form factor of an UNKNOWN atom.\n"
                    "This typically occurs when performing species-dependent operations on data without form factor information."
                );
            }
            if (!form_factor::detail::is_tabulated(type)) {
                throw ausaxs::except::runtime_error("form_factor::lookup::atomic::normalized::get: Invalid form factor type (enum " + std::to_string(static_cast<int>(type)) + ")");
            }
            return detail::table[static_cast<int>(type)];
        }
    }
}
