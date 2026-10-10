#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <settings/All.h>

#include <support/hist_test_helper.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::form_factor;

static const std::vector<int>& identity() {
    static std::vector<int> identity;
    if (identity.empty()) {
        identity = std::vector<int>(total_ff_count);
        std::iota(identity.begin(), identity.end(), 0);
    }
    return identity;
}

// the form factor set selected for the molecule, with all but EXCLUDED_VOLUME, WATER, and OTHER in random order
static std::vector<int> shuffled(const data::Molecule& protein) {
    manager::use_form_factors(protein);
    const auto* tables = manager::get_active_product_tables();
    std::vector<int> selected(tables->ff_indices.begin(), tables->ff_indices.begin() + tables->active_count);

    auto shuffled = selected;
    std::mt19937 g(std::random_device{}());
    std::shuffle(shuffled.begin()+2, shuffled.end()-1, g);
    if (shuffled == selected) {std::reverse(shuffled.begin()+2, shuffled.end()-1);} // make sure the slots actually move
    return shuffled;
}

template<template<bool> class MANAGER>
static void run_comparison(const data::Molecule& protein) {
    manager::use_form_factors(protein);
    auto i1 = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2 = MANAGER<true>(&protein).calculate_all()->debye_transform();

    manager::detail::use_form_factors(shuffled(protein));
    auto i1s = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2s = MANAGER<true>(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i1, i1s));
    REQUIRE(compare_hist(i2, i2s));
}

template<typename MANAGER>
static void run_comparison(const data::Molecule& protein) {
    manager::use_form_factors(protein);
    auto i = MANAGER(&protein).calculate_all()->debye_transform();

    manager::detail::use_form_factors(shuffled(protein));
    auto is = MANAGER(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i, is));
}

TEST_CASE("manager ff set change scattering consistent across all managers") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    invoke_for_all_histogram_manager_variants(
        []<typename MANAGER>(const data::Molecule& protein) {
            run_comparison<MANAGER>(protein);
        },
        []<template<bool> class MANAGER>(const data::Molecule& protein) {
            run_comparison<MANAGER>(protein);
        },
        protein
    );
}

TEST_CASE("manager ff set change scattering consistent for special exv calculators") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    auto run_exv = [&](settings::exv::ExvMethod method) {
        settings::exv::exv_method = method;
        run_comparison<hist::HistogramManagerMTFFExplicit>(protein);
    };

    SECTION("FoXS") {
        run_exv(settings::exv::ExvMethod::FoXS);
    }

    SECTION("Pepsi") {
        run_exv(settings::exv::ExvMethod::Pepsi);
    }

    SECTION("CRYSOL") {
        run_exv(settings::exv::ExvMethod::CRYSOL);
    }

    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
}

// The molecule-derived form factor set is truncated to the types the molecule actually contains, which
// also shrinks every histogram dimension indexed by form factor type. The dropped slots are all-zero, so
// the scattering must be unchanged - this verifies the truncated allocations and the runtime packed-index
// stride agree with each other in every manager.
template<template<bool> class MANAGER>
static void run_truncation_comparison(data::Molecule& protein) {
    manager::detail::use_form_factors(identity());
    auto i1 = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2 = MANAGER<true>(&protein).calculate_all()->debye_transform();

    manager::use_form_factors(protein);
    REQUIRE(get_active_count() < form_factor::total_ff_count);
    auto i1t = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2t = MANAGER<true>(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i1, i1t));
    REQUIRE(compare_hist(i2, i2t));
}

template<typename MANAGER>
static void run_truncation_comparison(data::Molecule& protein) {
    manager::detail::use_form_factors(identity());
    auto i = MANAGER(&protein).calculate_all()->debye_transform();

    manager::use_form_factors(protein);
    REQUIRE(get_active_count() < form_factor::total_ff_count);
    auto it = MANAGER(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i, it));
}

TEST_CASE("form_factor_manager: truncated ff set scattering consistent across all managers") {
    // the full identity selection is the untruncated reference, so it needs every slot, and only the absent types may be dropped
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    settings::form_factor::max_types = total_ff_count;
    settings::form_factor::min_fraction = 0;
    settings::general::verbose = false;

    auto run = [] () {
        data::Molecule protein("tests/files/2epe.pdb");
        protein.generate_new_hydration();

        invoke_for_all_histogram_manager_variants(
            []<typename MANAGER>(data::Molecule& protein) {
                run_truncation_comparison<MANAGER>(protein);
            },
            []<template<bool> class MANAGER>(data::Molecule& protein) {
                run_truncation_comparison<MANAGER>(protein);
            },
            protein
        );
    };

    SECTION("implicit hydrogens") { // 2epe contains neither H nor any of the bare elements beyond C/N/O/S, so the set shrinks
        settings::molecule::implicit_hydrogens = true;
        run();
    }

    SECTION("explicit hydrogens") { // only C/N/O/S are present, so the set shrinks to 7
        settings::molecule::implicit_hydrogens = false;
        run();
    }

    settings::molecule::implicit_hydrogens = true;
    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}

TEST_CASE("form_factor_manager: truncated ff set scattering consistent for special exv calculators") {
    // the FoXS product tables are only filled over the active sub-block, so they need the same check as the histograms.
    // Pepsi and CRYSOL share the manager tables, but switch them to the Traube volumes - and all of these are only reachable through their exv models

    // the full identity selection is the untruncated reference, so it needs every slot, and only the absent types may be dropped
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    settings::form_factor::max_types = total_ff_count;
    settings::form_factor::min_fraction = 0;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    auto run_exv = [&](settings::exv::ExvMethod method) {
        settings::exv::exv_method = method;
        run_truncation_comparison<hist::HistogramManagerMTFFExplicit>(protein);
    };

    SECTION("FoXS")   {run_exv(settings::exv::ExvMethod::FoXS);}
    SECTION("Pepsi")  {run_exv(settings::exv::ExvMethod::Pepsi);}
    SECTION("CRYSOL") {run_exv(settings::exv::ExvMethod::CRYSOL);}

    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}

TEST_CASE("form_factor_manager: use_form_factors(Molecule) reproduces identity scattering") {
    // the full identity selection is the untruncated reference, so it needs every slot, and only the absent types may be dropped
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    settings::form_factor::max_types = total_ff_count;
    settings::form_factor::min_fraction = 0;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    // baseline: the full identity form factor ordering
    manager::detail::use_form_factors(identity());
    auto I = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();

    // the form factor set selected from the molecular composition drops the types the molecule does
    // not contain, and those slots are all-zero, so the resulting scattering must be unchanged
    manager::use_form_factors(protein);
    auto I2 = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(I, I2));
    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}

// Everything reaching a histogram through a Molecule - the API, pyAUSAXS, the rigidbody optimizer, the EM fitter and the CLI alike - builds its
// manager through the factory, so that is where the form factor set is selected. A manager constructed by hand leaves the caller's selection alone.
TEST_CASE("form_factor_manager: the factory selects the molecule's form factor set") {
    // the full identity selection is the untruncated reference, so it needs every slot, and only the absent types may be dropped
    auto original_max = settings::form_factor::max_types;
    auto original_fraction = settings::form_factor::min_fraction;
    settings::form_factor::max_types = total_ff_count;
    settings::form_factor::min_fraction = 0;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Average;
    manager::detail::use_form_factors(identity());

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    SECTION("a hand-built manager keeps the caller's set") {
        auto I = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();
        REQUIRE(get_active_count() == total_ff_count);
        CHECK(I.size() != 0);
    }

    SECTION("a manager from the factory truncates to the molecule") {
        auto I = hist::HistogramManagerMTFFAvg<false>(&protein).calculate_all()->debye_transform();
        REQUIRE(get_active_count() == total_ff_count);

        auto manager = hist::factory::construct_histogram_manager(&protein, settings::hist::HistogramManagerChoice::HistogramManagerMT, false);
        REQUIRE(get_active_count() < total_ff_count);

        // the dropped slots are empty, so the profile must be the one the full set produced
        REQUIRE(compare_hist(I, manager->calculate_all()->debye_transform()));
    }

    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    settings::molecule::implicit_hydrogens = true;
    settings::form_factor::max_types = original_max;
    settings::form_factor::min_fraction = original_fraction;
}
