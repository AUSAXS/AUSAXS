#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("PrecalculatedFormFactorProduct::evaluate") {
    for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
            const NormalizedFormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
            const NormalizedFormFactor& ff2_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff2));
            PrecalculatedFormFactorProduct ff(ff1_obj, ff2_obj);
            for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
            }
        }
    }
}

TEST_CASE("PrecalculatedFormFactorProduct::table") {
    const auto& table = storage::atomic::get_precalculated_form_factor_table();
    for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
            const NormalizedFormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
            const NormalizedFormFactor& ff2_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff2));
            const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
            for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
            }
        }
    }
}