// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactor.h>
#include <form_factor/NormalizedFormFactor.h>

namespace ausaxs::form_factor::xray::detail {
    struct RawFormFactorLookup {
        static constexpr const FormFactor& get(form_factor_t type) {
            return raw::get(type);
        }
    };

    struct NormalizedFormFactorLookup {
        static constexpr const NormalizedFormFactor& get(form_factor_t type) {
            return normalized::get(type);
        }
    };
}