#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("PrecalculatedFormFactorProduct::evaluate") {
    SECTION("exv") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                PrecalculatedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                PrecalculatedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}

TEST_CASE("PrecalculatedFormFactorProduct::table") {
    SECTION("exv") {
        const auto& table = storage::exv::get_precalculated_form_factor_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        const auto& table = storage::cross::get_precalculated_form_factor_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const FormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = storage::exv::standard.get_form_factor(static_cast<form_factor_t>(ff2));
                const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactor: switch volumes") {
    auto test = [] (const constants::displaced_volume::detail::DisplacedVolumeSet& vols) {
        SECTION("exv") {
            const auto& table = storage::exv::get_precalculated_form_factor_table();
            auto ffset = form_factor::detail::ExvFormFactorSet(vols);
            for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                    const ExvFormFactor& ff1_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff1));
                    const ExvFormFactor& ff2_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff2));
                    const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                    for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                    }
                }
            }
        }

        SECTION("cross") {
            const auto& table = storage::cross::get_precalculated_form_factor_table();
            auto ffset = form_factor::detail::ExvFormFactorSet(vols);
            for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                    const FormFactor& ff1_obj = storage::atomic::get_form_factor(static_cast<form_factor_t>(ff1));
                    const ExvFormFactor& ff2_obj = ffset.get_form_factor(static_cast<form_factor_t>(ff2));
                    const PrecalculatedFormFactorProduct& ff = table.index(ff1, ff2);
                    for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                    }
                }
            }
        }
    };

    SECTION("Traube") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Traube;
        test(constants::displaced_volume::Traube);
    }

    SECTION("Voronoi_explicit_H") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Voronoi_explicit_H;
        test(constants::displaced_volume::Voronoi_explicit_H);
    }

    SECTION("Voronoi_implicit_H") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::Voronoi_implicit_H;
        test(constants::displaced_volume::Voronoi_implicit_H);
    }

    SECTION("MinimumFluctutation_explicit_H") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::MinimumFluctutation_explicit_H;
        test(constants::displaced_volume::MinimumFluctuation_explicit_H);
    }

    SECTION("MinimumFluctutation_implicit_H") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::MinimumFluctutation_implicit_H;
        test(constants::displaced_volume::MinimumFluctuation_implicit_H);
    }

    SECTION("vdw") {
        settings::molecule::displaced_volume_set = settings::molecule::DisplacedVolumeSet::vdw;
        test(constants::displaced_volume::vdw);
    }
}