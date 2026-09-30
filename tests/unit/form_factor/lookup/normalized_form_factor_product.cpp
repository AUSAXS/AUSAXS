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
    // activate a small selection containing C and H, and return their slots
    std::pair<int, int> use_C_and_H() {
        manager::detail::use_form_factors({
            static_cast<int>(form_factor_t::EXCLUDED_VOLUME),
            static_cast<int>(form_factor_t::WATER),
            static_cast<int>(form_factor_t::H),
            static_cast<int>(form_factor_t::C),
            static_cast<int>(form_factor_t::OTHER)
        });
        auto mapping = manager::get_active_mapping();
        return {mapping[static_cast<int>(form_factor_t::C)], mapping[static_cast<int>(form_factor_t::H)]};
    }
}

TEST_CASE("NormalizedFormFactorProduct::constructor") {
    SECTION("from two NormalizedFormFactors") {
        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        const xray::NormalizedFormFactor& H = xray::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff(C, H);
        CHECK(ff.evaluate(0) > 0);
    }

    auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
    SECTION("from NormalizedFormFactor and ExvFormFactor") {
        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        const ExvFormFactor& exv = exv_set.get(form_factor_t::C);
        
        NormalizedFormFactorProduct ff(C, exv);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("from two ExvFormFactors") {
        const ExvFormFactor& exv1 = exv_set.get(form_factor_t::C);
        const ExvFormFactor& exv2 = exv_set.get(form_factor_t::N);
        
        NormalizedFormFactorProduct ff(exv1, exv2);
        CHECK(ff.evaluate(0) > 0);
    }
}

TEST_CASE("NormalizedFormFactorProduct::evaluate") {
    SECTION("matches manual calculation") {
        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        const xray::NormalizedFormFactor& H = xray::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff(C, H);
        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * H.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("symmetric") {
        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        const xray::NormalizedFormFactor& H = xray::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff1(C, H);
        NormalizedFormFactorProduct ff2(H, C);
        
        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }

    SECTION("same form factor squared") {
        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        
        NormalizedFormFactorProduct ff(C, C);
        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double c_val = C.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(c_val * c_val, 1e-10));
        }
    }
}

TEST_CASE("NormalizedFormFactorProduct::all_pairs") {
    SECTION("all atomic form factor pairs") {
        for (int ff1 = 0; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 0; ff2 < total_ff_count; ++ff2) {
                const xray::NormalizedFormFactor& ff1_obj = xray::normalized::get(static_cast<form_factor_t>(ff1));
                const xray::NormalizedFormFactor& ff2_obj = xray::normalized::get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("manager::raw_atomic_table") {
    auto [c_slot, h_slot] = use_C_and_H();
    const auto& table = manager::get_active_product_tables()->raw_atomic_table;
    SECTION("single access") {
        const auto& ff = table.index(c_slot, h_slot);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = table.index(c_slot, h_slot);
        const auto& ff2 = table.index(h_slot, c_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("manager::normalized_atomic_table") {
    SECTION("table access") {
        auto [c_slot, h_slot] = use_C_and_H();
        const auto& table = manager::get_active_product_tables()->normalized_atomic_table;

        const xray::NormalizedFormFactor& C = xray::normalized::get(form_factor_t::C);
        const xray::NormalizedFormFactor& H = xray::normalized::get(form_factor_t::H);

        const auto& ff = table.index(c_slot, h_slot);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * H.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        test::form_factor::use_random_form_factors();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->normalized_atomic_table;
        for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
            for (int ff2 = 0; ff2 < tables->active_count; ++ff2) {
                const xray::NormalizedFormFactor& ff1_obj = xray::normalized::get(static_cast<form_factor_t>(tables->ff_indices[ff1]));
                const xray::NormalizedFormFactor& ff2_obj = xray::normalized::get(static_cast<form_factor_t>(tables->ff_indices[ff2]));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}
