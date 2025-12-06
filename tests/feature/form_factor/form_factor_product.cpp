#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/lookup/FormFactorProduct.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::comprehensive_evaluation") {
    SECTION("all form factor products match direct calculation") {
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const FormFactor& ff1_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff1));
                const FormFactor& ff2_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff2));
                FormFactorProduct ff(ff1_obj, ff2_obj);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table_comprehensive") {
    SECTION("all table entries match direct calculation") {
        const auto& table = lookup::atomic::raw::get_table();
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const FormFactor& ff1_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff1));
                const FormFactor& ff2_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::specific_pairs") {
    SECTION("H-H product") {
        const FormFactor& ff = lookup::atomic::raw::H;
        FormFactorProduct ffp(ff, ff);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N product") {
        const FormFactor& ff_c = lookup::atomic::raw::C;
        const FormFactor& ff_n = lookup::atomic::raw::N;
        FormFactorProduct ffp(ff_c, ff_n);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("O-S product") {
        const FormFactor& ff_o = lookup::atomic::raw::O;
        const FormFactor& ff_s = lookup::atomic::raw::S;
        FormFactorProduct ffp(ff_o, ff_s);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) > ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}
