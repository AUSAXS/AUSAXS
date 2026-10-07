#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <support/form_factor_helper.h>

#include <utility>

using namespace ausaxs;
using namespace form_factor;

namespace {
    // activate a small selection containing C and N, and return their slots
    std::pair<int, int> use_C_and_N() {
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

    // types without a volume in the current set have an empty exv profile
    double evaluate_exv(const form_factor::detail::ExvFormFactorSet& exv_set, form_factor_t type, double q) {
        return exv_set.contains(type) ? exv_set.get(type).evaluate(q) : 0.0;
    }
}

TEST_CASE("manager::raw_exv_table") {
    auto [c_slot, n_slot] = use_C_and_N();
    const auto& table = manager::get_active_product_tables()->raw_exv_table;
    SECTION("single access") {
        const auto& ff = table.index(c_slot, n_slot);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = table.index(c_slot, n_slot);
        const auto& ff2 = table.index(n_slot, c_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("manager::raw_exv_table: completeness") {
    SECTION("table access") {
        auto [c_slot, n_slot] = use_C_and_N();
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        const ExvFormFactor& C = exv_set.get(form_factor_t::C);
        const ExvFormFactor& N = exv_set.get(form_factor_t::N);
        const auto& ff = table.index(c_slot, n_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        test::form_factor::use_random_form_factors();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = evaluate_exv(exv_set, t1, constants::axes::q_vals[i]) * evaluate_exv(exv_set, t2, constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10) || Catch::Matchers::WithinAbs(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("manager::normalized_cross_table") {
    auto [c_slot, n_slot] = use_C_and_N();
    const auto& table = manager::get_active_product_tables()->normalized_cross_table;
    auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
    SECTION("single access") {
        const auto& ff = table.index(c_slot, n_slot);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

        const auto& ff = table.index(c_slot, n_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }
}

TEST_CASE("manager::normalized_cross_table: completeness") {
    SECTION("table access") {
        auto [c_slot, n_slot] = use_C_and_N();
        const auto& table = manager::get_active_product_tables()->normalized_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

        const auto& ff = table.index(c_slot, n_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        test::form_factor::use_random_form_factors();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->normalized_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(tables->ff_indices[ff1]));
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * evaluate_exv(exv_set, t2, constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10) || Catch::Matchers::WithinAbs(expected, 1e-10));
                }
            }
        }
    }
}
