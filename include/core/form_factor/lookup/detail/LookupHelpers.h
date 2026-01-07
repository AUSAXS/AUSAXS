// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactor.h>
#include <form_factor/NormalizedFormFactor.h>

namespace ausaxs::form_factor::lookup::detail {
    struct RawFormFactorLookup {
        static constexpr const FormFactor& get(form_factor_t type) {
            return lookup::atomic::raw::get(type);
        }
    };

    struct NormalizedFormFactorLookup {
        static constexpr const NormalizedFormFactor& get(form_factor_t type) {
            return lookup::atomic::normalized::get(type);
        }
    };
}