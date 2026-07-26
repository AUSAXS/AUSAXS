#include <catch2/catch_test_macros.hpp>

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <settings/All.h>

#include "hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

template<template<bool> class MANAGER> constexpr settings::hist::HistogramManagerChoice choice_for();
template<template<bool, bool> class MANAGER> constexpr settings::hist::HistogramManagerChoice choice_for();

template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManager>() {return settings::hist::HistogramManagerChoice::HistogramManager;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMT>() {return settings::hist::HistogramManagerChoice::HistogramManagerMT;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMTFFAvg>() {return settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMTFFExplicit>() {return settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::SymmetryManagerMT>() {return settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::PartialHistogramManager>() {return settings::hist::HistogramManagerChoice::PartialHistogramManager;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::PartialHistogramManagerMT>() {return settings::hist::HistogramManagerChoice::PartialHistogramManagerMT;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::PartialSymmetryManagerMT>() {return settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMTFFGrid>() {return settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMTFFGridSurface>() {return settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface;}
template<> constexpr settings::hist::HistogramManagerChoice choice_for<hist::HistogramManagerMTFFGridScalableExv>() {return settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv;}

TEST_CASE("HistogramManagerFactory: resolves partial and symmetry preferences") {
    auto exv = settings::exv::exv_method;
    auto threads = settings::general::threads;
    auto prefer_partial = settings::flags::prefer_partial_manager;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    settings::general::threads = 4;
    settings::hist::weighted_bins = true;
    settings::flags::custom_bin_width = false;

    Molecule plain({Body{SimpleCube::get_atoms()}});

    Molecule symmetric({Body{SimpleCube::get_atoms()}});
    symmetric.get_body(0).symmetry().add(symmetry::type::c2);
    REQUIRE(symmetric.symmetry().has_symmetries());

    SECTION("without a partial preference") {
        settings::flags::prefer_partial_manager = false;

        // symmetry-awareness is derived from the molecule, never from a setting
        CHECK(dynamic_cast<hist::HistogramManagerMT<true, false>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::SymmetryManagerMT<true, false>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);
    }

    SECTION("with a partial preference") {
        settings::flags::prefer_partial_manager = true;

        CHECK(dynamic_cast<hist::PartialHistogramManagerMT<true, false>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::PartialSymmetryManagerMT<true, false>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);

        // the single-threaded partial manager has no symmetry-aware counterpart, so it upgrades to the MT one
        settings::general::threads = 1;
        CHECK(dynamic_cast<hist::PartialHistogramManager<true, false>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::PartialSymmetryManagerMT<true, false>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);
    }

    SECTION("preference is dropped when the excluded volume method has no partial implementation") {
        settings::flags::prefer_partial_manager = true;
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;

        // the excluded volume model wins: it changes the result, whereas dropping the partial preference only costs time
        CHECK(dynamic_cast<hist::HistogramManagerMTFFExplicit<true, false>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
    }

    settings::exv::exv_method = exv;
    settings::general::threads = threads;
    settings::flags::prefer_partial_manager = prefer_partial;
}

TEST_CASE("HistogramManagerFactory: creates expected manager") {
    Molecule protein({Body{SimpleCube::get_atoms()}});

    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein) {
            settings::flags::custom_bin_width = false;
            auto hm_w = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<false>*>(hm_w.get()) != nullptr);

            settings::flags::custom_bin_width = true;
            auto hm = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<true>*>(hm.get()) != nullptr);
        },
        []<template<bool, bool> class MANAGER>(const Molecule& protein) {
            settings::flags::custom_bin_width = false;
            settings::hist::weighted_bins = false;
            auto hm = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<false, false>*>(hm.get()) != nullptr);

            settings::hist::weighted_bins = true;
            auto hm_w = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<true, false>*>(hm_w.get()) != nullptr);

            settings::flags::custom_bin_width = true;
            settings::hist::weighted_bins = false;
            auto hm_vbw = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<false, true>*>(hm_vbw.get()) != nullptr);

            settings::hist::weighted_bins = true;
            auto hm_w_vbw = hist::factory::construct_histogram_manager(&protein, choice_for<MANAGER>());
            REQUIRE(dynamic_cast<MANAGER<true, true>*>(hm_w_vbw.get()) != nullptr);
        },
        protein
    );
}