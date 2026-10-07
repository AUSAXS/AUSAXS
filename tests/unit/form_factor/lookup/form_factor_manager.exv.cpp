#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/FormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <support/form_factor_helper.h>

#include <utility>
#include <vector>

using namespace ausaxs;
using namespace form_factor;

namespace {
    // the active slots from index 1 onwards, paired with their form factor types. 
    // with_volume restricts them to the types which have an excluded volume in the current set.
    std::vector<std::pair<int, form_factor_t>> active_slots(bool with_volume) {
        const auto* tables = manager::get_active_product_tables();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        std::vector<std::pair<int, form_factor_t>> slots;
        for (int i = 1; i < tables->active_count; ++i) {
            auto type = static_cast<form_factor_t>(tables->ff_indices[i]);
            if (with_volume && !exv_set.contains(type)) {continue;}
            slots.emplace_back(i, type);
        }
        return slots;
    }
}

TEST_CASE("ExvFormFactorProduct::raw_exv_table") {
    test::form_factor::use_random_form_factors();

    SECTION("product entries match direct calculation") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (auto [i, ti] : active_slots(true)) {
            for (auto [j, tj] : active_slots(true)) {
                const FormFactorProduct& product = table.index(i, j);

                ExvFormFactor exv1 = exv_set.get(ti);
                ExvFormFactor exv2 = exv_set.get(tj);
                
                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = exv1.evaluate(constants::axes::q_vals[k]) * exv2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("all table entries match direct calculation") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (auto [i, ti] : active_slots(true)) {
            for (auto [j, tj] : active_slots(true)) {
                const FormFactorProduct& product = table.index(i, j);
                
                ExvFormFactor exv1 = exv_set.get(ti);
                ExvFormFactor exv2 = exv_set.get(tj);
                
                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = exv1.evaluate(constants::axes::q_vals[k]) * exv2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::table symmetry") {
    test::form_factor::use_random_form_factors();

    SECTION("exv table is symmetric") {
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        for (auto [i, ti] : active_slots(false)) {
            for (auto [j, tj] : active_slots(false)) {
                const FormFactorProduct& product1 = table.index(i, j);
                const FormFactorProduct& product2 = table.index(j, i);
                
                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    CHECK_THAT(product1.evaluate(k), Catch::Matchers::WithinRel(product2.evaluate(k), 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::cross products") {
    test::form_factor::use_random_form_factors();

    SECTION("cross products return correct values") {
        const auto& table = manager::get_active_product_tables()->raw_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (auto [i, ti] : active_slots(false)) {
            for (auto [j, tj] : active_slots(true)) {
                const FormFactorProduct& product = table.index(i, j);
                const FormFactor& ff_atomic = lookup::atomic::raw::get(ti);
                ExvFormFactor exv = exv_set.get(tj);

                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff_atomic.evaluate(constants::axes::q_vals[k]) * exv.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("cross table returns correct values") {
        const auto& table = manager::get_active_product_tables()->raw_cross_table;
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (auto [i, ti] : active_slots(false)) {
            for (auto [j, tj] : active_slots(true)) {
                const FormFactorProduct& product = table.index(i, j);
                const FormFactor& ff_atomic = lookup::atomic::raw::get(ti);
                ExvFormFactor exv = exv_set.get(tj);

                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff_atomic.evaluate(constants::axes::q_vals[k]) * exv.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::product decreases with q") {
    test::form_factor::use_random_form_factors();

    SECTION("exv products decrease") {
        // slot 0 (EXCLUDED_VOLUME) has no explicit exv entries, so use WATER which is always present
        const auto& table = manager::get_active_product_tables()->raw_exv_table;
        const FormFactorProduct& product = table.index(water_bin, water_bin);

        double val1 = product.evaluate(0);
        double val2 = product.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = product.evaluate(constants::axes::q_axis.bins - 1);

        CHECK(val1 >= val2);
        CHECK(val2 >= val3);
    }

    SECTION("cross products decrease") {
        const auto& table = manager::get_active_product_tables()->raw_cross_table;
        const FormFactorProduct& product = table.index(water_bin, water_bin);

        double val1 = product.evaluate(0);
        double val2 = product.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = product.evaluate(constants::axes::q_axis.bins - 1);

        CHECK(val1 >= val2);
        CHECK(val2 >= val3);
    }
}
