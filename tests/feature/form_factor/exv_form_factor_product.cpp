#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvFormFactorProduct.h>
#include <form_factor/ExvTable.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvFormFactorProduct::comprehensive_exv_evaluation") {
    SECTION("all exv form factor products match direct calculation") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                ExvFormFactor exv1 = lookup::exv::standard.get(static_cast<form_factor_t>(ff1));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = lookup::exv::raw::get_product(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = exv1.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::comprehensive_cross_evaluation") {
    SECTION("all cross form factor products match direct calculation") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff1));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = lookup::cross::raw::get_product(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::exv_table_comprehensive") {
    SECTION("all exv table entries match direct calculation") {
        const auto& table = lookup::exv::raw::get_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                ExvFormFactor exv1 = lookup::exv::standard.get(static_cast<form_factor_t>(ff1));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = exv1.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::cross_table_comprehensive") {
    SECTION("all cross table entries match direct calculation") {
        const auto& table = lookup::cross::raw::get_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = lookup::atomic::raw::get(static_cast<form_factor_t>(ff1));
                ExvFormFactor exv2 = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * exv2.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactorProduct::specific_exv_pairs") {
    SECTION("C exv form factor product") {
        const FormFactorProduct& ffp = lookup::exv::raw::get_product(1, 1);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N exv cross product") {
        const FormFactorProduct& ffp = lookup::exv::raw::get_product(1, 2);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}

TEST_CASE("ExvFormFactorProduct::specific_cross_pairs") {
    SECTION("C atomic-exv cross product") {
        const FormFactorProduct& ffp = lookup::cross::raw::get_product(1, 1);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }

    SECTION("C-N atomic-exv cross product") {
        const FormFactorProduct& ffp = lookup::cross::raw::get_product(1, 2);
        
        CHECK(ffp.evaluate(0) > 0);
        CHECK(ffp.evaluate(constants::axes::q_axis.bins - 1) > 0);
        CHECK(ffp.evaluate(0) >= ffp.evaluate(constants::axes::q_axis.bins - 1));
    }
}
