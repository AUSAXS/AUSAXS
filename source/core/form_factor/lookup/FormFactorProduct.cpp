// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/FormFactor.h>
#include <constants/Constants.h>

#if CONSTEXPR_TABLES
    #define CONST constexpr
#else
    #define CONST const
#endif

using namespace ausaxs;
using namespace ausaxs::form_factor;

namespace {
    CONST form_factor::lookup::atomic::table_t generate_table() {
        form_factor::lookup::atomic::table_t table;
        for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                table.index(i, j) = FormFactorProduct(
                    lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
                    lookup::atomic::raw::get(static_cast<form_factor_t>(j))
                );
                table.index(j, i) = table.index(i, j);
            }
            table.index(i, i) = FormFactorProduct(
                lookup::atomic::raw::get(static_cast<form_factor_t>(i)), 
                lookup::atomic::raw::get(static_cast<form_factor_t>(i))
            );
        }
        return table;
    }

    auto ff_table = generate_table();
}

const FormFactorProduct& form_factor::lookup::atomic::raw::get_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const form_factor::lookup::atomic::table_t& form_factor::lookup::atomic::raw::get_table() noexcept {
    return ff_table;
}