#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <settings/ExvSettings.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("manager::raw_exv_table") {
    const auto& table = manager::get_active_product_tables()->raw_exv_table;
    SECTION("single access") {
        const auto& ff = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );
        const auto& ff2 = table.index(
            static_cast<int>(form_factor_t::N),
            static_cast<int>(form_factor_t::C)
        );

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("manager::raw_exv_table: completeness") {
    SECTION("table access") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        const ExvFormFactor& C = exv_set.get(form_factor_t::C);
        const ExvFormFactor& N = exv_set.get(form_factor_t::N);
        const auto& ff = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (int ff1 = 1; ff1 < form_factor::total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < form_factor::total_ff_count; ++ff2) {
                const ExvFormFactor& ff1_obj = exv_set.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("manager::normalized_cross_table") {
    const auto& table = manager::get_active_product_tables()->normalized_cross_table;
    auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
    SECTION("single access") {
        const auto& ff = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

        const auto& ff = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }
}

TEST_CASE("manager::normalized_cross_table: completeness") {
    SECTION("table access") {
        const auto& table = manager::get_active_product_tables()->normalized_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

        const auto& ff = table.index(
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::N)
        );

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        const auto& table = manager::get_active_product_tables()->normalized_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (int ff1 = 0; ff1 < form_factor::total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < form_factor::total_ff_count; ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}
