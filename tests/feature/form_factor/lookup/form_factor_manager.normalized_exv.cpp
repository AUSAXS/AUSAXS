#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::evaluate") {
    SECTION("exv") {
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = exv_set.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (int ff1 = 0; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table") {
    SECTION("exv") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = exv_set.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        const auto& table = manager::get_active_product_tables()->normalized_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (int ff1 = 0; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}
