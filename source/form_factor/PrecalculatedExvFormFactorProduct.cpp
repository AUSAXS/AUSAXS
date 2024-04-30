/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <constants/Constants.h>

#include <cmath>

using namespace form_factor;

constexpr form_factor::storage::exv::table_t generate_exv_exv_table() {
    form_factor::storage::exv::table_t table;
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = PrecalculatedFormFactorProduct(
                storage::exv::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::exv::get_form_factor(static_cast<form_factor_t>(j))
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = PrecalculatedFormFactorProduct(
            storage::exv::get_form_factor(static_cast<form_factor_t>(i)), 
            storage::exv::get_form_factor(static_cast<form_factor_t>(i))
        );
    }
    return table;
}

auto ff_xx_table = generate_exv_exv_table();
const PrecalculatedFormFactorProduct& form_factor::storage::exv::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_xx_table.index(i, j);
}

const form_factor::storage::exv::table_t& form_factor::storage::exv::get_precalculated_form_factor_table() noexcept {
    return ff_xx_table;
}

constexpr form_factor::storage::cross::table_t generate_atom_exv_table() {
    form_factor::storage::cross::table_t table;
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
        for (unsigned int j = 0; j < form_factor::get_count_without_excluded_volume(); ++j) {
            table.index(i, j) = PrecalculatedFormFactorProduct(
                storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::exv::get_form_factor(static_cast<form_factor_t>(j))
            );
        }
    }
    return table;
}

auto ff_ax_table = generate_atom_exv_table();
const PrecalculatedFormFactorProduct& form_factor::storage::cross::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_ax_table.index(i, j);
}

const form_factor::storage::cross::table_t& form_factor::storage::cross::get_precalculated_form_factor_table() noexcept {
    return ff_ax_table;
}