#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::constructor") {
    SECTION("construct from two FormFactors") {
        const xray::FormFactor& ff1 = xray::raw::get(form_factor_t::H);
        const xray::FormFactor& ff2 = xray::raw::get(form_factor_t::C);
        FormFactorProduct ffp(ff1, ff2);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = ff1.evaluate(constants::axes::q_vals[i]) * ff2.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ffp.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("construct from same FormFactor") {
        const xray::FormFactor& ff = xray::raw::get(form_factor_t::C);
        FormFactorProduct ffp(ff, ff);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double ff_val = ff.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ffp.evaluate(i), Catch::Matchers::WithinRel(ff_val * ff_val, 1e-10));
        }
    }
}

TEST_CASE("FormFactorProduct::evaluate") {
    SECTION("product decreases with q") {
        const xray::FormFactor& ff1 = xray::raw::get(form_factor_t::C);
        const xray::FormFactor& ff2 = xray::raw::get(form_factor_t::N);
        FormFactorProduct ffp(ff1, ff2);

        double val1 = ffp.evaluate(0);
        double val2 = ffp.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = ffp.evaluate(constants::axes::q_axis.bins - 1);

        CHECK(val1 > val2);
        CHECK(val2 > val3);
    }

    SECTION("product is positive") {
        const xray::FormFactor& ff1 = xray::raw::get(form_factor_t::O);
        const xray::FormFactor& ff2 = xray::raw::get(form_factor_t::S);
        FormFactorProduct ffp(ff1, ff2);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK(ffp.evaluate(i) > 0);
        }
    }
}

TEST_CASE("FormFactorProduct::symmetry") {
    SECTION("product is symmetric") {
        const xray::FormFactor& ff1 = xray::raw::get(form_factor_t::H);
        const xray::FormFactor& ff2 = xray::raw::get(form_factor_t::O);
        FormFactorProduct ffp1(ff1, ff2);
        FormFactorProduct ffp2(ff2, ff1);

        for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ffp1.evaluate(i), Catch::Matchers::WithinRel(ffp2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("FormFactorProduct::raw_atomic_table") {
    test::form_factor::use_random_form_factors();

    SECTION("product entries match direct calculation") {
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_atomic_table;
        for (int i = 0; i < tables->active_count; ++i) {
            for (int j = 0; j < tables->active_count; ++j) {
                const xray::FormFactor& ff1 = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[i]));
                const xray::FormFactor& ff2 = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[j]));
                const FormFactorProduct& product = table.index(i, j);

                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff1.evaluate(constants::axes::q_vals[k]) * ff2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("all table entries match direct calculation") {
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_atomic_table;
        for (int i = 0; i < tables->active_count; ++i) {
            for (int j = 0; j < tables->active_count; ++j) {
                const xray::FormFactor& ff1 = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[i]));
                const xray::FormFactor& ff2 = xray::raw::get(static_cast<form_factor_t>(tables->ff_indices[j]));
                const FormFactorProduct& product = table.index(i, j);

                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff1.evaluate(constants::axes::q_vals[k]) * ff2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table symmetry") {
    test::form_factor::use_random_form_factors();

    SECTION("table is symmetric") {
        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_atomic_table;
        for (int i = 0; i < tables->active_count; ++i) {
            for (int j = 0; j < tables->active_count; ++j) {
                const FormFactorProduct& product1 = table.index(i, j);
                const FormFactorProduct& product2 = table.index(j, i);

                for (int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    CHECK_THAT(product1.evaluate(k), Catch::Matchers::WithinRel(product2.evaluate(k), 1e-10));
                }
            }
        }
    }
}
