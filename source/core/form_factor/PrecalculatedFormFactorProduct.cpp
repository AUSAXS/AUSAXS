// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <constants/Constants.h>

#if CONSTEXPR_TABLES
    #define CONST constexpr
#else
    #define CONST const
#endif

using namespace ausaxs;
using namespace ausaxs::form_factor;

CONST form_factor::storage::atomic::table_t generate_table() {
    form_factor::storage::atomic::table_t table;
    for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = PrecalculatedFormFactorProduct(
                storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::atomic::get_form_factor(static_cast<form_factor_t>(j))
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = PrecalculatedFormFactorProduct(
            storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
            storage::atomic::get_form_factor(static_cast<form_factor_t>(i))
        );
    }
    return table;
}

auto ff_table = generate_table();
const PrecalculatedFormFactorProduct& form_factor::storage::atomic::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const form_factor::storage::atomic::table_t& form_factor::storage::atomic::get_precalculated_form_factor_table() noexcept {
    return ff_table;
}