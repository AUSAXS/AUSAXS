#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <form_factor/FormFactor.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/MoleculeSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("storage::exv::get_precalculated_form_factor_product") {
    SECTION("single access") {
        const auto& ff = storage::exv::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = storage::exv::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        const auto& ff2 = storage::exv::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::N),
            static_cast<unsigned int>(form_factor_t::C)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            CHECK_THAT(ff1.evaluate(i), Catch::Matchers::WithinRel(ff2.evaluate(i), 1e-10));
        }
    }
}

TEST_CASE("storage::exv::get_precalculated_form_factor_table") {
    SECTION("table access") {
        const auto& table = storage::exv::get_precalculated_form_factor_table();
        
        const ExvFormFactor& C = storage::exv::standard.get_form_factor(form_factor_t::C);
        const ExvFormFactor& N = storage::exv::standard.get_form_factor(form_factor_t::N);
        
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
        const auto& table = storage::exv::get_precalculated_form_factor_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("storage::cross::get_precalculated_form_factor_product") {
    SECTION("single access") {
        const auto& ff = storage::cross::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("matches manual calculation") {
        const FormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const ExvFormFactor& N_exv = storage::exv::standard.get_form_factor(form_factor_t::N);
        
        const auto& ff = storage::cross::get_precalculated_form_factor_product(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        
        for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
            double expected = C.evaluate(constants::axes::q_vals[i]) * N_exv.evaluate(constants::axes::q_vals[i]);
            CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
        }
    }
}

TEST_CASE("storage::cross::get_precalculated_form_factor_table") {
    SECTION("table access") {
        const auto& table = storage::cross::get_precalculated_form_factor_table();
        
        const FormFactor& C = storage::atomic::get_form_factor(form_factor_t::C);
        const ExvFormFactor& N_exv = storage::exv::standard.get_form_factor(form_factor_t::N);
        
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
        const auto& table = storage::cross::get_precalculated_form_factor_table();
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("storage::detail::set_custom_exv_table") {
    SECTION("set custom table") {
        auto original_setting = settings::molecule::exv_set;
        settings::molecule::exv_set = settings::molecule::ExvSet::Custom;
        
        constants::exv::detail::ExvSet custom_set = constants::exv::vdw;
        storage::detail::set_custom_exv_table(custom_set);
        
        const auto& table_exv = storage::exv::get_precalculated_form_factor_table();
        const auto& table_cross = storage::cross::get_precalculated_form_factor_table();
        
        CHECK(table_exv.index(0, 0).evaluate(0) > 0);
        CHECK(table_cross.index(0, 0).evaluate(0) > 0);
        
        settings::molecule::exv_set = original_setting;
    }
}

TEST_CASE("ExvSet switching") {
    SECTION("Traube") {
        auto original_setting = settings::molecule::exv_set;
        settings::molecule::exv_set = settings::molecule::ExvSet::Traube;
        
        const auto& table = storage::exv::get_precalculated_form_factor_table();
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::Traube);
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                
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
        
        const auto& table = storage::cross::get_precalculated_form_factor_table();
        auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::vdw);
        
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        
        settings::molecule::exv_set = original_setting;
    }
}
