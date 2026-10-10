#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <settings/ExvSettings.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvTableManager::set_custom_exv_table") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();

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

        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(custom_set);
        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                for (int i = 0; i < 10; ++i) {
                    double expected = ffset.get(t1).evaluate(constants::axes::q_vals[i]) * ffset.get(t2).evaluate(constants::axes::q_vals[i]);
                    REQUIRE_THAT(table.index(ff1, ff2).evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        settings::exv::exv_set = original_setting;
    }
}

TEST_CASE("ExvSet switching") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();

    SECTION("Traube") {
        auto original_setting = settings::exv::exv_set.value;
        settings::exv::exv_set = settings::exv::ExvSet::Traube;

        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::Traube);

        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const FormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ffset.get(t1).evaluate(constants::axes::q_vals[i]) * ffset.get(t2).evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }

        settings::exv::exv_set = original_setting;
    }

    SECTION("vdw") {
        auto original_setting = settings::exv::exv_set.value;
        settings::exv::exv_set = settings::exv::ExvSet::vdw;

        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::vdw);

        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const FormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ffset.get(t1).evaluate(constants::axes::q_vals[i]) * ffset.get(t2).evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }

        settings::exv::exv_set = original_setting;
    }
}
