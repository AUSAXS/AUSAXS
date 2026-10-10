#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/ExvSettings.h>
#include <support/form_factor_helper.h>

#include <utility>

using namespace ausaxs;
using namespace form_factor;

namespace {
    // activate a small set containing C and N, and return their active slots
    std::pair<int, int> use_carbon_and_nitrogen() {
        manager::detail::use_form_factors({
            static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
            static_cast<int>(form_factor_t::WATER),
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N),
            static_cast<int>(form_factor_t::OTHER)
        });
        auto mapping = manager::get_active_mapping();
        return {mapping[static_cast<int>(form_factor_t::C)], mapping[static_cast<int>(form_factor_t::N)]};
    }
}

TEST_CASE("ExvFormFactorProduct::comprehensive_exv_evaluation") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models

    SECTION("all exv form factor products match direct calculation") {
        test::form_factor::use_random_form_factors();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                ExvFormFactor exv1 = exv_set.get(static_cast<form_factor_t>(tables->ff_indices[ff1]));
                ExvFormFactor exv2 = exv_set.get(static_cast<form_factor_t>(tables->ff_indices[ff2]));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = exv1.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::comprehensive_cross_evaluation") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models

    SECTION("all cross form factor products match direct calculation") {
        test::form_factor::use_random_form_factors();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_cross_table;
        for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                const xray::FormFactor& ff1_obj = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[ff1]));
                ExvFormFactor exv2 = exv_set.get(static_cast<form_factor_t>(tables->ff_indices[ff2]));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::specific_exv_pairs") {
    auto [c, n] = use_carbon_and_nitrogen();

    SECTION("C exv form factor product") {
        const FormFactorProduct& ffp = manager::get_active_product_tables()->raw_exv_table.index(c, c);
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N exv cross product") {
        const FormFactorProduct& ffp = manager::get_active_product_tables()->raw_exv_table.index(c, n);
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}

TEST_CASE("ExvFormFactorProduct::specific_cross_pairs") {
    auto [c, n] = use_carbon_and_nitrogen();

    SECTION("C atomic-exv cross product") {
        const FormFactorProduct& ffp = manager::get_active_product_tables()->raw_cross_table.index(c, c);
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N atomic-exv cross product") {
        const FormFactorProduct& ffp = manager::get_active_product_tables()->raw_cross_table.index(c, n);
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}
