#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <constants/Constants.h>

#include <cmath>

using namespace form_factor;

constexpr container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> generate_exv_exv_table() {
    container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> table;
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = std::move(
                PrecalculatedFormFactorProduct(
                    storage::exv::get_form_factor(static_cast<form_factor_t>(i)), 
                    storage::exv::get_form_factor(static_cast<form_factor_t>(j))
                )
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = std::move(
            PrecalculatedFormFactorProduct(
                storage::exv::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::exv::get_form_factor(static_cast<form_factor_t>(i))
            )
        );
    }
    return table;
}

inline constexpr auto ff_table = generate_exv_exv_table();
const PrecalculatedFormFactorProduct& form_factor::storage::exv::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>& form_factor::storage::exv::get_precalculated_form_factor_table() noexcept {
    return ff_table;
}

constexpr container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> generate_atom_exv_table() {
    container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()> table;
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = std::move(
                PrecalculatedFormFactorProduct(
                    storage::get_form_factor(static_cast<form_factor_t>(i)), 
                    storage::exv::get_form_factor(static_cast<form_factor_t>(j))
                )
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = std::move(
            PrecalculatedFormFactorProduct(
                storage::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::exv::get_form_factor(static_cast<form_factor_t>(i))
            )
        );
    }
    return table;
}

inline constexpr auto ff_table = generate_atom_exv_table();
const PrecalculatedFormFactorProduct& form_factor::storage::cross::get_precalculated_form_factor_product(unsigned int i, unsigned int j) noexcept {
    return ff_table.index(i, j);
}

const container::ArrayContainer2D<PrecalculatedFormFactorProduct, form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume()>& form_factor::storage::cross::get_precalculated_form_factor_table() noexcept {
    return ff_table;
}