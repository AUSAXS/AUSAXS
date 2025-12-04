#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::constructor") {
    SECTION("construct from two FormFactors") {
        const FormFactor& ff1 = lookup::atomic::raw::H;
        const FormFactor& ff2 = lookup::atomic::raw::C;
        FormFactorProduct ffp(ff1, ff2);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = ff1.evaluate(constants::axes::q_vals[i]) * ff2.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ffp.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("construct from same FormFactor") {
        const FormFactor& ff = lookup::atomic::raw::C;
        FormFactorProduct ffp(ff, ff);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double ff_val = ff.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ffp.evaluate(i), Catch::Matchers::WithinRel(ff_val * ff_val, 1e-10));
        }
    }
}

TEST_CASE("FormFactorProduct::evaluate") {
    SECTION("product decreases with q") {
        const FormFactor& ff1 = lookup::atomic::raw::C;
        const FormFactor& ff2 = lookup::atomic::raw::N;
        FormFactorProduct ffp(ff1, ff2);
        
        double val1 = ffp.evaluate(0);
        double val2 = ffp.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = ffp.evaluate(constants::axes::q_axis.bins - 1);
        
        CHECK(val1 > val2);
        CHECK(val2 > val3);
    }

    SECTION("product is positive") {
        const FormFactor& ff1 = lookup::atomic::raw::O;
        const FormFactor& ff2 = lookup::atomic::raw::S;
        FormFactorProduct ffp(ff1, ff2);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK(ffp.evaluate(i) > 0);
        }
    }
}

TEST_CASE("FormFactorProduct::symmetry") {
    SECTION("product is symmetric") {
        const FormFactor& ff1 = lookup::atomic::raw::H;
        const FormFactor& ff2 = lookup::atomic::raw::O;
        FormFactorProduct ffp1(ff1, ff2);
        FormFactorProduct ffp2(ff2, ff1);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ffp1.evaluate(i), Catch::Matchers::WithinRel(ffp2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("FormFactorProduct::lookup::atomic::raw") {
    SECTION("get_product returns correct product") {
        for (unsigned int i = 0; i < get_count(); ++i) {
            for (unsigned int j = 0; j < get_count(); ++j) {
                const FormFactor& ff1 = lookup::atomic::raw::get(static_cast<form_factor_t>(i));
                const FormFactor& ff2 = lookup::atomic::raw::get(static_cast<form_factor_t>(j));
                const FormFactorProduct& product = lookup::atomic::raw::get_product(i, j);
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff1.evaluate(constants::axes::q_vals[k]) * ff2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("get_table returns correct table") {
        const auto& table = lookup::atomic::raw::get_table();
        
        for (unsigned int i = 0; i < get_count(); ++i) {
            for (unsigned int j = 0; j < get_count(); ++j) {
                const FormFactor& ff1 = lookup::atomic::raw::get(static_cast<form_factor_t>(i));
                const FormFactor& ff2 = lookup::atomic::raw::get(static_cast<form_factor_t>(j));
                const FormFactorProduct& product = table.index(i, j);
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff1.evaluate(constants::axes::q_vals[k]) * ff2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table symmetry") {
    SECTION("table is symmetric") {
        const auto& table = lookup::atomic::raw::get_table();
        
        for (unsigned int i = 0; i < get_count(); ++i) {
            for (unsigned int j = 0; j < get_count(); ++j) {
                const FormFactorProduct& product1 = table.index(i, j);
                const FormFactorProduct& product2 = table.index(j, i);
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    CHECK_THAT(product1.evaluate(k), Catch::Matchers::WithinRel(product2.evaluate(k), 1e-10));
                }
            }
        }
    }
}
