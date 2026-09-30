#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/NormalizedFormFactor.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <settings/All.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvFormFactor: switch volumes") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();

    auto test = [&] (const constants::exv::detail::ExvSet& vols) {
        auto ffset = form_factor::detail::ExvFormFactorSet(vols);
        auto exv = [&ffset] (form_factor_t type, double q) {return ffset.get(type).evaluate(q);};

        SECTION("exv") {
            const auto* tables = manager::get_active_product_tables();
            const auto& table = tables->raw_exv_table;
            for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
                for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                    auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                    auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                    const FormFactorProduct& ff = table.index(ff1, ff2);
                    for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        double expected = exv(t1, constants::axes::q_vals[i])*exv(t2, constants::axes::q_vals[i]);
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected));
                    }
                }
            }
        }

        SECTION("cross") {
            const auto* tables = manager::get_active_product_tables();
            const auto& table = tables->normalized_cross_table;
            for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
                for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                    auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                    auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                    const xray::NormalizedFormFactor& ff1_obj = xray::normalized::get(t1);
                    const FormFactorProduct& ff = table.index(ff1, ff2);
                    for (int i = 0; i < constants::axes::q_axis.bins; ++i) {
                        double expected = ff1_obj.evaluate(constants::axes::q_vals[i])*exv(t2, constants::axes::q_vals[i]);
                        REQUIRE_THAT(ff.evaluate(i), Catch::Matchers::WithinRel(expected));
                    }
                }
            }
        }
    };

    SECTION("Traube") {
        settings::exv::exv_set = settings::exv::ExvSet::Traube;
        test(constants::exv::Traube);
    }

    SECTION("Voronoi_explicit_H") {
        settings::exv::exv_set = settings::exv::ExvSet::Voronoi_explicit_H;
        test(constants::exv::Voronoi_explicit_H);
    }

    SECTION("Voronoi_implicit_H") {
        settings::exv::exv_set = settings::exv::ExvSet::Voronoi_implicit_H;
        test(constants::exv::Voronoi_implicit_H);
    }

    SECTION("MinimumFluctutation_explicit_H") {
        settings::exv::exv_set = settings::exv::ExvSet::MinimumFluctutation_explicit_H;
        test(constants::exv::MinimumFluctuation_explicit_H);
    }

    SECTION("MinimumFluctutation_implicit_H") {
        settings::exv::exv_set = settings::exv::ExvSet::MinimumFluctutation_implicit_H;
        test(constants::exv::MinimumFluctuation_implicit_H);
    }

    SECTION("vdw") {
        settings::exv::exv_set = settings::exv::ExvSet::vdw;
        test(constants::exv::vdw);
    }

    settings::exv::exv_set = settings::exv::ExvSet::Default;
}
