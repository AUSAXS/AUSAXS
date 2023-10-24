#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <constants/Constants.h>

#include <cmath>

using namespace form_factor;

constexpr container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()> generate_table() {
    container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()> table;
    for (unsigned int i = 0; i < form_factor::get_count(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = PrecalculatedFormFactorProduct(
                storage::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::get_form_factor(static_cast<form_factor_t>(j))
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = PrecalculatedFormFactorProduct(
            storage::get_form_factor(static_cast<form_factor_t>(i)), 
            storage::get_form_factor(static_cast<form_factor_t>(i))
        );
    }
    return table;
}

inline constexpr auto ff_table = generate_table();
const PrecalculatedFormFactorProduct& form_factor::storage::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count(), form_factor::get_count()>& form_factor::storage::get_precalculated_form_factor_table() noexcept {
    return ff_table;
}