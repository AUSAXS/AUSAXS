#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <form_factor/ExvFormFactor.h>
#include <form_factor/lookup/ExvTableManager.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/lookup/FormFactorProduct.h>
#include <form_factor/xray/XrayFormFactor.h>
#include <settings/ExvSettings.h>
#include <settings/GeneralSettings.h>
#include <support/form_factor_helper.h>

using namespace ausaxs;
using namespace form_factor;

TEST_CASE("ExvTableManager::set_custom_exv_table") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();

    SECTION("set custom table") {
        constants::exv::detail::ExvSet custom_set = constants::exv::vdw;
        ExvTableManager::set_custom_exv_table(custom_set);

        const auto& table_exv = manager::get_active_product_tables()->raw_exv_table;
        const auto& table_cross = manager::get_active_product_tables()->raw_cross_table;

        const int other = manager::get_active_product_tables()->active_count - 1; // OTHER is always the last active slot
        REQUIRE(table_exv.index(other, other).evaluate(0) > 0);
        REQUIRE(table_cross.index(other, other).evaluate(0) > 0);
    }

    SECTION("products match direct calculation after custom table is set") {
        constants::exv::detail::ExvSet custom_set = constants::exv::Traube;
        ExvTableManager::set_custom_exv_table(custom_set);

        const auto* tables = manager::get_active_product_tables();
        const auto& table = tables->raw_exv_table;
        auto ffset = form_factor::detail::ExvFormFactorSet(custom_set, constants::charge::density::water);
        for (int ff1 = start_index_for_explicit_exv(); ff1 < tables->active_count; ++ff1) {
            for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                for (int i = 0; i < 10; ++i) {
                    double expected = test::form_factor::evaluate_exv(ffset, t1, constants::axes::q_vals[i]) * test::form_factor::evaluate_exv(ffset, t2, constants::axes::q_vals[i]);
                    REQUIRE_THAT(table.index(ff1, ff2).evaluate(i), Catch::Matchers::WithinRel(expected, 1e-10));
                }
            }
        }
    }
}

TEST_CASE("ExvSet switching") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();

    auto test = [&] (const constants::exv::detail::ExvSet& vols) {
        auto ffset = form_factor::detail::ExvFormFactorSet(vols, constants::charge::density::water);
        auto exv = [&ffset] (form_factor_t type, double q) {return test::form_factor::evaluate_exv(ffset, type, q);};

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
            const auto& table = tables->raw_cross_table;
            for (int ff1 = 0; ff1 < tables->active_count; ++ff1) {
                for (int ff2 = start_index_for_explicit_exv(); ff2 < tables->active_count; ++ff2) {
                    auto t1 = static_cast<form_factor_t>(tables->ff_indices[ff1]);
                    auto t2 = static_cast<form_factor_t>(tables->ff_indices[ff2]);
                    const xray::FormFactor& ff1_obj = xray::raw::get(t1);
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
}

TEST_CASE("ExvTableManager: kept explicit hydrogens use the explicit-H tables") {
    settings::exv::exv_method = settings::exv::ExvMethod::Fraser; // the explicit exv tables are only used by the Fraser-based models
    test::form_factor::use_random_form_factors();
    using Set = settings::exv::ExvSet;

    SECTION("the implicit-H sets switch to their explicit-H variants") {
        settings::exv::exv_set = Set::MinimumFluctutation_implicit_H;
        settings::general::keep_hydrogens = false;
        CHECK(*ExvTableManager::get_current_exv_table() == constants::exv::MinimumFluctuation_implicit_H);
        settings::general::keep_hydrogens = true;
        CHECK(*ExvTableManager::get_current_exv_table() == constants::exv::MinimumFluctuation_explicit_H);

        settings::exv::exv_set = Set::Voronoi_implicit_H;
        CHECK(*ExvTableManager::get_current_exv_table() == constants::exv::Voronoi_explicit_H);
    }

    SECTION("sets without hydrogen variants are unchanged") {
        settings::general::keep_hydrogens = true;
        settings::exv::exv_set = Set::Traube;
        CHECK(*ExvTableManager::get_current_exv_table() == constants::exv::Traube);
        settings::exv::exv_set = Set::vdw;
        CHECK(*ExvTableManager::get_current_exv_table() == constants::exv::vdw);
    }

    SECTION("toggling keep_hydrogens rebuilds the active tables") {
        settings::exv::exv_set = Set::MinimumFluctutation_implicit_H;
        settings::general::keep_hydrogens = false;
        // the active exv tables must follow the switch without a new form factor selection
        auto exv_at_zero = [] () {
            const auto* tables = manager::get_active_product_tables();
            double sum = 0;
            for (int ff = start_index_for_explicit_exv(); ff < tables->active_count; ++ff) {sum += tables->raw_exv_table.index(ff, ff).evaluate(0);}
            return sum;
        };
        double before = exv_at_zero();
        settings::general::keep_hydrogens = true;
        CHECK(exv_at_zero() != before);
    }
}
