#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/ExvTable.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvFormFactorProduct::lookup::exv::raw") {
    SECTION("get_product returns correct product") {
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < get_count_without_excluded_volume(); ++j) {
                const FormFactorProduct& product = lookup::exv::raw::get_product(i, j);
                
                ExvFormFactor exv1 = lookup::exv::standard.get(static_cast<form_factor_t>(i));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(j));
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = exv1.evaluate(constants::axes::q_vals[k]) * exv2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("get_table returns correct table") {
        const auto& table = lookup::exv::raw::get_table();
        
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < get_count_without_excluded_volume(); ++j) {
                const FormFactorProduct& product = table.index(i, j);
                
                ExvFormFactor exv1 = lookup::exv::standard.get(static_cast<form_factor_t>(i));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(j));
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = exv1.evaluate(constants::axes::q_vals[k]) * exv2.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::table symmetry") {
    SECTION("exv table is symmetric") {
        const auto& table = lookup::exv::raw::get_table();
        
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < get_count_without_excluded_volume(); ++j) {
                const FormFactorProduct& product1 = table.index(i, j);
                const FormFactorProduct& product2 = table.index(j, i);
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    CHECK_THAT(product1.evaluate(k), Catch::Matchers::WithinRel(product2.evaluate(k), 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::cross products") {
    SECTION("cross products return correct values") {
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < get_count_without_excluded_volume(); ++j) {
                const FormFactorProduct& product = lookup::cross::raw::get_product(i, j);
                
                const FormFactor& ff_atomic = lookup::atomic::raw::get(static_cast<form_factor_t>(i));
                ExvFormFactor exv = lookup::exv::standard.get(static_cast<form_factor_t>(j));
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff_atomic.evaluate(constants::axes::q_vals[k]) * exv.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }

    SECTION("cross table returns correct values") {
        const auto& table = lookup::cross::raw::get_table();
        
        for (unsigned int i = 0; i < get_count_without_excluded_volume(); ++i) {
            for (unsigned int j = 0; j < get_count_without_excluded_volume(); ++j) {
                const FormFactorProduct& product = table.index(i, j);
                
                const FormFactor& ff_atomic = lookup::atomic::raw::get(static_cast<form_factor_t>(i));
                ExvFormFactor exv = lookup::exv::standard.get(static_cast<form_factor_t>(j));
                
                for (unsigned int k = 0; k < constants::axes::q_axis.bins; ++k) {
                    double expected = ff_atomic.evaluate(constants::axes::q_vals[k]) * exv.evaluate(constants::axes::q_vals[k]);
                    CHECK_THAT(product.evaluate(k), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::product decreases with q") {
    SECTION("exv products decrease") {
        const FormFactorProduct& product = lookup::exv::raw::get_product(0, 0);
        
        double val1 = product.evaluate(0);
        double val2 = product.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = product.evaluate(constants::axes::q_axis.bins - 1);
        
        CHECK(val1 >= val2);
        CHECK(val2 >= val3);
    }

    SECTION("cross products decrease") {
        const FormFactorProduct& product = lookup::cross::raw::get_product(0, 0);
        
        double val1 = product.evaluate(0);
        double val2 = product.evaluate(constants::axes::q_axis.bins / 2);
        double val3 = product.evaluate(constants::axes::q_axis.bins - 1);
        
        CHECK(val1 >= val2);
        CHECK(val2 >= val3);
    }
}
