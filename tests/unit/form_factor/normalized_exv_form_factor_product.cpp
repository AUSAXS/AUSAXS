#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/ExvFormFactor.h>
#include <settings/ExvSettings.h>
#include <constants/Constants.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("lookup::exv::normalized::get_product") {
    auto& table = FormFactorManager::normalized_exv_table();
    SECTION("single access") {
        const auto& ff = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("symmetric access") {
        const auto& ff1 = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        const auto& ff2 = table.index(
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
        auto& table = FormFactorManager::normalized_exv_table();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
        const ExvFormFactor& C = exv_set.get(form_factor_t::C);
        const ExvFormFactor& N = exv_set.get(form_factor_t::N);
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
        const auto& table = FormFactorManager::normalized_exv_table();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = exv_set.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
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
    auto& table = FormFactorManager::normalized_cross_table();
    auto exv_set = ExvTableManager::get_current_exv_form_factor_set();
    SECTION("single access") {
        const auto& ff = table.index(
            static_cast<unsigned int>(form_factor_t::C),
            static_cast<unsigned int>(form_factor_t::N)
        );
        CHECK(ff.evaluate(0) > 0);
    }

    SECTION("matches manual calculation") {
        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

        const auto& ff = table.index(
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
        const auto& table = FormFactorManager::normalized_cross_table();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        const NormalizedFormFactor& C = lookup::atomic::normalized::get(form_factor_t::C);
        const ExvFormFactor& N_exv = exv_set.get(form_factor_t::N);

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
        auto& table = FormFactorManager::normalized_cross_table();
        auto exv_set = ExvTableManager::get_current_exv_form_factor_set();

        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = exv_set.get(static_cast<form_factor_t>(ff2));
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
        auto original_setting = settings::exv::exv_set;

        constants::exv::detail::ExvSet custom_set = constants::exv::vdw;
        ExvTableManager::set_custom_exv_table(custom_set);

        const auto& table_exv = FormFactorManager::normalized_exv_table();
        const auto& table_cross = FormFactorManager::raw_cross_table();

        REQUIRE(table_exv.index(0, 0).evaluate(0) > 0);
        REQUIRE(table_cross.index(0, 0).evaluate(0) > 0);

        settings::exv::exv_set = original_setting;
    }

    SECTION("raw and normalized tables stay in sync") {
        auto original_setting = settings::exv::exv_set;

        // Set custom table using normalized interface
        constants::exv::detail::ExvSet custom_set = constants::exv::Traube;
        ExvTableManager::set_custom_exv_table(custom_set);

        // Verify both raw and normalized tables are updated
        const auto& norm_table = FormFactorManager::normalized_exv_table();
        const auto& raw_table = FormFactorManager::raw_exv_table();

        // The excluded volume form factors are the same for raw and normalized
        // (normalization only affects atomic form factors)
        auto ffset = form_factor::detail::ExvFormFactorSet(custom_set);
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                
                for (unsigned int i = 0; i < 10; ++i) { // Check first 10 q-values
                    double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
                    REQUIRE_THAT(norm_table.index(ff1, ff2).evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                    REQUIRE_THAT(raw_table.index(ff1, ff2).evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
        settings::exv::exv_set = original_setting;
    }
}

TEST_CASE("ExvSet switching") {
    SECTION("Traube") {
        auto original_setting = settings::exv::exv_set;
        settings::exv::exv_set = settings::exv::ExvSet::Traube;

        const auto& table = FormFactorManager::raw_exv_table();
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

        settings::exv::exv_set = original_setting;
    }

    // SECTION("vdw") {
    //     auto original_setting = settings::exv::exv_set;
    //     settings::exv::exv_set = settings::exv::ExvSet::vdw;

    //     const auto& table = FormFactorManager::raw_cross_table();
    //     auto ffset = form_factor::detail::ExvFormFactorSet(constants::exv::vdw);

    //     for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
    //         for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
    //             const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
    //             const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
    //             const NormalizedFormFactorProduct& ff = table.index(ff1, ff2);

    //             for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
    //                 double expected = ff1_obj.evaluate(constants::axes::q_vals[i]) * ff2_obj.evaluate(constants::axes::q_vals[i]);
    //                 CHECK_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
    //             }
    //         }
    //     }

    //     settings::exv::exv_set = original_setting;
    // }
}
