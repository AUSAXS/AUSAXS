// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/lookup/FormFactorProductBase.h>
#include <form_factor/NormalizedFormFactor.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    struct NormalizedFormFactorLookup {
        static constexpr const NormalizedFormFactor& get(form_factor_t type) {
            return lookup::atomic::normalized::get(type);
        }
    };

    auto ff_table = lookup::detail::generate_atomic_table<NormalizedFormFactorLookup>();
}

const NormalizedFormFactorProduct& form_factor::lookup::atomic::normalized::get_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const form_factor::lookup::atomic::table_t& form_factor::lookup::atomic::normalized::get_table() noexcept {
    return ff_table;
}