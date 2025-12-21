// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/lookup/FormFactorProductBase.h>
#include <form_factor/FormFactor.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    struct RawFormFactorLookup {
        static constexpr const FormFactor& get(form_factor_t type) {
            return lookup::atomic::raw::get(type);
        }
    };

    auto ff_table = lookup::detail::generate_atomic_table<RawFormFactorLookup>();
}

const FormFactorProduct& form_factor::lookup::atomic::raw::get_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const form_factor::lookup::atomic::table_t& form_factor::lookup::atomic::raw::get_table() noexcept {
    return ff_table;
}