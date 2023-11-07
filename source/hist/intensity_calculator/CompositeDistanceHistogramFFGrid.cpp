#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/GridSettings.h>

using namespace hist;
using namespace form_factor;

form_factor::storage::atomic::table_t CompositeDistanceHistogramFFGrid::generate_table() {
    form_factor::storage::atomic::table_t table;

    auto scaling_factor = settings::grid::exv_radius/settings::grid::width;
    FormFactor ffx = ExvFormFactor(std::pow(settings::grid::width*scaling_factor, 3));
    for (unsigned int i = 0; i < form_factor::get_count_without_excluded_volume(); ++i) {
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

        table.index(i, form_factor::exv_bin) = PrecalculatedFormFactorProduct(
            storage::atomic::get_form_factor(static_cast<form_factor_t>(i)), 
            ffx
        );
        table.index(form_factor::exv_bin, i) = table.index(i, form_factor::exv_bin);
        table.index(form_factor::exv_bin, form_factor::exv_bin) = PrecalculatedFormFactorProduct(
            ffx, 
            ffx
        );
    }
    return table;
}