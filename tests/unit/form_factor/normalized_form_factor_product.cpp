#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("NormalizedFormFactorProduct::constructor") {
    SECTION("from two NormalizedFormFactors") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const NormalizedFormFactor& H = lookup::atomic::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff(C, H);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("from NormalizedFormFactor and ExvFormFactor") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& exv = lookup::exv::standard.get(form_factor_t::C);
        
        NormalizedFormFactorProduct ff(C, exv);
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("from two ExvFormFactors") {
        const ExvFormFactor& exv1 = lookup::exv::standard.get(form_factor_t::C);
        const ExvFormFactor& exv2 = lookup::exv::standard.get(form_factor_t::N);
        
        NormalizedFormFactorProduct ff(exv1, exv2);
        CHECK(ff.evaluate(0) > 0);
    }
}

TEST_CASE("NormalizedFormFactorProduct::evaluate") {
    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const NormalizedFormFactor& H = lookup::atomic::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff(C, H);
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * H.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("symmetric") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const NormalizedFormFactor& H = lookup::atomic::normalized::get(form_factor_t::H);
        
        NormalizedFormFactorProduct ff1(C, H);
        NormalizedFormFactorProduct ff2(H, C);
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }

    SECTION("same form factor squared") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        
        NormalizedFormFactorProduct ff(C, C);
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double c_val = C.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(c_val * c_val, 1e-10));
        }
    }
}

TEST_CASE("NormalizedFormFactorProduct::all_pairs") {
    SECTION("all atomic form factor pairs") {
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const NormalizedFormFactor& ff2_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("lookup::atomic::normalized::get_product") {
    SECTION("single access") {
        const auto& ff = lookup::atomic::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::H)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = lookup::atomic::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::H)
        );
        const auto& ff2 = lookup::atomic::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::H),
            static_cast<unsigned int>(form_factor_t::C)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("lookup::atomic::normalized::get_table") {
    SECTION("table access") {
        const auto& table = lookup::atomic::normalized::get_table();
        
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const NormalizedFormFactor& H = lookup::atomic::normalized::get(form_factor_t::H);
        
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
        const auto& table = lookup::atomic::normalized::get_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const NormalizedFormFactor& ff2_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}
