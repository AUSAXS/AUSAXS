#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <settings/ExvSettings.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvTableManager::set_custom_exv_table") {
    SECTION("set custom table") {
        auto original_setting = settings::exv::exv_set.value;

        constants::exv::detail::ExvSet custom_set = constants::exv::vdw;
        ExvTableManager::set_custom_exv_table(custom_set);

        const auto& table_exv = manager::get_active_product_tables()->raw_exv_table;
        const auto& table_cross = manager::get_active_product_tables()->raw_cross_table;

        REQUIRE(table_exv.index(form_factor::water_bin, form_factor::water_bin).evaluate(0) > 0);
        REQUIRE(table_cross.index(form_factor::water_bin, form_factor::water_bin).evaluate(0) > 0);

        settings::exv::exv_set = original_setting;
    }

    SECTION("products match direct calculation after custom table is set") {
        auto original_setting = settings::exv::exv_set.value;

        constants::exv::detail::ExvSet custom_set = constants::exv::Traube;
        ExvTableManager::set_custom_exv_table(custom_set);

        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(custom_set);
        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                for (int i = 0; i < 10; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    REQUIRE_THAT(table.index(ff1, ff2).evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        settings::exv::exv_set = original_setting;
    }
}

TEST_CASE("ExvSet switching") {
    SECTION("Traube") {
        auto original_setting = settings::exv::exv_set.value;
        settings::exv::exv_set = settings::exv::ExvSet::Traube;

        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::Traube);

        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }

        settings::exv::exv_set = original_setting;
    }

    SECTION("vdw") {
        auto original_setting = settings::exv::exv_set.value;
        settings::exv::exv_set = settings::exv::ExvSet::vdw;

        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::vdw);

        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }

        settings::exv::exv_set = original_setting;
    }
}
