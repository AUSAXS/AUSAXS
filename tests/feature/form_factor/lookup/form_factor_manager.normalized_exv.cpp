#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <settings/All.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::evaluate") {
    SECTION("exv") {
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        for (int ff1 = 1; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 1; ff2 < total_ff_count; ++ff2) {
                // not every type has a volume in the current set
                if (!exv_set.contains(static_cast<form_factor_t>(ff1)) || !exv_set.contains(static_cast<form_factor_t>(ff2))) {continue;}
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
                if (!exv_set.contains(static_cast<form_factor_t>(ff2))) {continue;}
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
    test::form_factor::use_random_form_factors();
    const auto* tables = manager::get_active_product_tables();

    // types without a volume in the current set have an empty exv profile
    auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
    auto exv = [&exv_set] (form_factor_t type, double q) {
        return exv_set.contains(type) ? exv_set.get(type).evaluate(q) : 0.0;
    };

    SECTION("exv") {
        const auto& table = tables->raw_exv_table;
        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = exv(t1, constants::axes::q_vals[i])*exv(t2, constants::axes::q_vals[i]);
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected) || Catch::Matchers::WithinAbs(expected, 1e-12));
                }
            }
        }
    }

    SECTION("cross") {
        const auto& table = tables->normalized_cross_table;
        for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(t1);
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i])*exv(t2, constants::axes::q_vals[i]);
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected) || Catch::Matchers::WithinAbs(expected, 1e-12));
                }
            }
        }
    }
}
