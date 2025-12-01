#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/NormalizedExvFormFactorProduct.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("lookup::exv::normalized::get_product") {
    SECTION("single access") {
        const auto& ff = lookup::exv::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = lookup::exv::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        const auto& ff2 = lookup::exv::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::N),
            static_cast<unsigned int>(form_factor_t::C)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("lookup::exv::normalized::get_table") {
    SECTION("table access") {
        const auto& table = lookup::exv::normalized::get_table();
        
        const ExvFormFactor& C = lookup::exv::standard.get(form_factor_t::C);
        const ExvFormFactor& N = lookup::exv::standard.get(form_factor_t::N);
        
        const auto& ff = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        const auto& table = lookup::exv::normalized::get_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("lookup::cross::normalized::get_product") {
    SECTION("single access") {
        const auto& ff = lookup::cross::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = lookup::exv::standard.get(form_factor_t::N);
        
        const auto& ff = lookup::cross::normalized::get_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }
}

TEST_CASE("lookup::cross::normalized::get_table") {
    SECTION("table access") {
        const auto& table = lookup::cross::normalized::get_table();
        
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = lookup::exv::standard.get(form_factor_t::N);
        
        const auto& ff = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }

    SECTION("table completeness") {
        const auto& table = lookup::cross::normalized::get_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("lookup::detail::set_custom_exv_table") {
    SECTION("set custom table") {
        auto original_setting = settings::molecule::exv_set;
        settings::molecule::exv_set = settings::molecule::ExvSet::Custom;
        
        constants::exv::detail::ExvSet custom_set = constants::exv::vdw;
        lookup::exv::normalized::detail::set_custom_table(custom_set);
        
        const auto& table_exv = lookup::exv::normalized::get_table();
        const auto& table_cross = lookup::cross::normalized::get_table();
        
        CHECK(table_exv.index(0, 0).evaluate(0) > 0);
        CHECK(table_cross.index(0, 0).evaluate(0) > 0);
        
        settings::molecule::exv_set = original_setting;
    }
}

TEST_CASE("ExvSet switching") {
    SECTION("Traube") {
        auto original_setting = settings::molecule::exv_set;
        settings::molecule::exv_set = settings::molecule::ExvSet::Traube;
        
        const auto& table = lookup::exv::normalized::get_table();
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::Traube);
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        
        settings::molecule::exv_set = original_setting;
    }

    SECTION("vdw") {
        auto original_setting = settings::molecule::exv_set;
        settings::molecule::exv_set = settings::molecule::ExvSet::vdw;
        
        const auto& table = lookup::cross::normalized::get_table();
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::vdw);
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        
        settings::molecule::exv_set = original_setting;
    }
}
