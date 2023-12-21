#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <constants/Constants.h>

#include <cmath>

using namespace form_factor;

constexpr form_factor::storage::atomic::table_t generate_table() {
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

constexpr auto ff_table = generate_table();
const PrecalculatedFormFactorProduct& form_factor::storage::atomic::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const form_factor::storage::atomic::table_t& form_factor::storage::atomic::get_precalculated_form_factor_table() noexcept {
    return ff_table;
}