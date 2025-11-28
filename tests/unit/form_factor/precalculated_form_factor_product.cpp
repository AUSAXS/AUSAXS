#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("PrecalculatedFormFactorProduct::constructor") {
    SECTION("from two NormalizedFormFactors") {
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const NormalizedFormFactor& H = storage::atomic::get_form_factor(form_factor_t::H);
        
        PrecalculatedFormFactorProduct ff(C, H);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("from NormalizedFormFactor and ExvFormFactor") {
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const ExvFormFactor& exv = storage::exv::standard.get_form_factor(form_factor_t::C);
        
        PrecalculatedFormFactorProduct ff(C, exv);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("from two ExvFormFactors") {
        const ExvFormFactor& exv1 = storage::exv::standard.get_form_factor(form_factor_t::C);
        const ExvFormFactor& exv2 = storage::exv::standard.get_form_factor(form_factor_t::N);
        
        PrecalculatedFormFactorProduct ff(exv1, exv2);
        CHECK(ff.evaluate(0) > 0);
    }
}

TEST_CASE("PrecalculatedFormFactorProduct::evaluate") {
    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const NormalizedFormFactor& H = storage::atomic::get_form_factor(form_factor_t::H);
        
        PrecalculatedFormFactorProduct ff(C, H);
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * H.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("symmetric") {
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const NormalizedFormFactor& H = storage::atomic::get_form_factor(form_factor_t::H);
        
        PrecalculatedFormFactorProduct ff1(C, H);
        PrecalculatedFormFactorProduct ff2(H, C);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }

    SECTION("same form factor squared") {
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        
        PrecalculatedFormFactorProduct ff(C, C);
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double c_val = C.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(c_val * c_val, 1e-10));
        }
    }
}

TEST_CASE("PrecalculatedFormFactorProduct::all_pairs") {
    SECTION("all atomic form factor pairs") {
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const NormalizedFormFactor& ff2_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff2));
                PrecalculatedFormFactorProduct ff(ff1_obj, ff2_obj);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("storage::atomic::get_precalculated_form_factor_product") {
    SECTION("single access") {
        const auto& ff = storage::atomic::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::H)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = storage::atomic::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::H)
        );
        const auto& ff2 = storage::atomic::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::H),
            static_cast<unsigned int>(form_factor_t::C)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("storage::atomic::get_precalculated_form_factor_table") {
    SECTION("table access") {
        const auto& table = storage::atomic::get_precalculated_form_factor_table();
        
        const NormalizedFormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const NormalizedFormFactor& H = storage::atomic::get_form_factor(form_factor_t::H);
        
        const auto& ff = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::H)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * H.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        const auto& table = storage::atomic::get_precalculated_form_factor_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const NormalizedFormFactor& ff2_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}
