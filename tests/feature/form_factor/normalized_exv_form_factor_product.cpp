#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/NormalizedFormFactorProduct.h>
#include <form_factor/lookup/NormalizedExvFormFactorProduct.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("FormFactorProduct::evaluate") {
    SECTION("exv") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                NormalizedFormFactorProduct ff(ff1_obj, ff2_obj);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}

TEST_CASE("FormFactorProduct::table") {
    SECTION("exv") {
        const auto& table = lookup::exv::normalized::get_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const ExvFormFactor& ff1_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }

    SECTION("cross") {
        const auto& table = lookup::cross::normalized::get_table();
        for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                const ExvFormFactor& ff2_obj = lookup::exv::standard.get(static_cast<form_factor_t>(ff2));
                const FormFactorProduct& ff = table.index(ff1, ff2);
                for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                    REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                }
            }
        }
    }
}

TEST_CASE("ExvFormFactor: switch volumes") {
    auto test = [] (const constants::exv::detail::ExvSet& vols) {
        SECTION("exv") {
            const auto& table = lookup::exv::normalized::get_table();
            auto ffset = form_factor::detail::ExvFormFactorSet(vols);
            for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                    const ExvFormFactor& ff1_obj = ffset.get(static_cast<form_factor_t>(ff1));
                    const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                    const FormFactorProduct& ff = table.index(ff1, ff2);
                    for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                    }
                }
            }
        }

        SECTION("cross") {
            const auto& table = lookup::cross::normalized::get_table();
            auto ffset = form_factor::detail::ExvFormFactorSet(vols);
            for (unsigned int ff1 = 0; ff1 < get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < get_count_without_excluded_volume(); ++ff2) {
                    const NormalizedFormFactor& ff1_obj = lookup::atomic::normalized::get(static_cast<form_factor_t>(ff1));
                    const ExvFormFactor& ff2_obj = ffset.get(static_cast<form_factor_t>(ff2));
                    const FormFactorProduct& ff = table.index(ff1, ff2);
                    for (unsigned int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(ff1_obj.evaluate(constants::axes::q_vals[i])*ff2_obj.evaluate(constants::axes::q_vals[i])));
                    }
                }
            }
        }
    };

    SECTION("Traube") {
        settings::molecule::exv_set = settings::molecule::ExvSet::Traube;
        test(constants::exv::Traube);
    }

    SECTION("Voronoi_explicit_H") {
        settings::molecule::exv_set = settings::molecule::ExvSet::Voronoi_explicit_H;
        test(constants::exv::Voronoi_explicit_H);
    }

    SECTION("Voronoi_implicit_H") {
        settings::molecule::exv_set = settings::molecule::ExvSet::Voronoi_implicit_H;
        test(constants::exv::Voronoi_implicit_H);
    }

    SECTION("MinimumFluctutation_explicit_H") {
        settings::molecule::exv_set = settings::molecule::ExvSet::MinimumFluctutation_explicit_H;
        test(constants::exv::MinimumFluctuation_explicit_H);
    }

    SECTION("MinimumFluctutation_implicit_H") {
        settings::molecule::exv_set = settings::molecule::ExvSet::MinimumFluctutation_implicit_H;
        test(constants::exv::MinimumFluctuation_implicit_H);
    }

    SECTION("vdw") {
        settings::molecule::exv_set = settings::molecule::ExvSet::vdw;
        test(constants::exv::vdw);
    }
}