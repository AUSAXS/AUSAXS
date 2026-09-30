#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::comprehensive_evaluation") {
    SECTION("all form factor products match direct calculation") {
        for (int ff1 = 0; ff1 < total_ff_count; ++ff1) {
            for (int ff2 = 0; ff2 < total_ff_count; ++ff2) {
                const xray::FormFactor& ff1_obj = xray::raw::get(static_cast<form_factor_t>(ff1));
                const xray::FormFactor& ff2_obj = xray::raw::get(static_cast<form_factor_t>(ff2));
                FormFactorProduct ff(ff1_obj, ff2_obj);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table_comprehensive") {
    SECTION("all table entries match direct calculation") {
        test::form_factor::use_random_form_factors();
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_atomic_table;
        for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
            for (int ff2 = 0; ff2 < tables->active_count; ++ff2) {
                const xray::FormFactor& ff1_obj = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[ff1]));
                const xray::FormFactor& ff2_obj = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[ff2]));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::specific_pairs") {
    SECTION("H-H product") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::H);
        FormFactorProduct ffp(ff, ff);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N product") {
        const xray::FormFactor& ff_c = xray::raw::get(form_factor_t::C);
        const xray::FormFactor& ff_n = xray::raw::get(form_factor_t::N);
        FormFactorProduct ffp(ff_c, ff_n);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("O-S product") {
        const xray::FormFactor& ff_o = xray::raw::get(form_factor_t::O);
        const xray::FormFactor& ff_s = xray::raw::get(form_factor_t::S);
        FormFactorProduct ffp(ff_o, ff_s);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}
