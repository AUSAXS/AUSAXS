#include <catch2/catch_test_macros.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/BodySymmetryFacade.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>

using namespace ausaxs;
using namespace ausaxs::data;

// the kind of manager and the excluded volume model which the factory turns into MANAGER
struct Request {
    settings::hist::HistogramManagerChoice kind;
    settings::exv::ExvMethod exv;
};
template<typename MANAGER> constexpr Request request_for();
template<template<bool> class MANAGER> constexpr Request request_for();

using Choice = settings::hist::HistogramManagerChoice;
using Exv = settings::exv::ExvMethod;
template<> constexpr Request request_for<hist::HistogramManager>() {return {Choice::HistogramManager, Exv::Simple};}
template<> constexpr Request request_for<hist::HistogramManagerMT>() {return {Choice::HistogramManagerMT, Exv::Simple};}
template<> constexpr Request request_for<hist::HistogramManagerMTFFAvg>() {return {Choice::HistogramManagerMT, Exv::Average};}
template<> constexpr Request request_for<hist::HistogramManagerMTFFExplicit>() {return {Choice::HistogramManagerMT, Exv::Fraser};}
template<> constexpr Request request_for<hist::SymmetryManagerMT>() {return {Choice::HistogramSymmetryManagerMT, Exv::Simple};}
template<> constexpr Request request_for<hist::PartialHistogramManager>() {return {Choice::PartialHistogramManager, Exv::Simple};}
template<> constexpr Request request_for<hist::PartialHistogramManagerMT>() {return {Choice::PartialHistogramManagerMT, Exv::Simple};}
template<> constexpr Request request_for<hist::PartialSymmetryManagerMT>() {return {Choice::PartialHistogramSymmetryManagerMT, Exv::Simple};}
template<> constexpr Request request_for<hist::HistogramManagerMTFFGrid>() {return {Choice::HistogramManagerMT, Exv::Grid};}
template<> constexpr Request request_for<hist::HistogramManagerMTFFGridSurface>() {return {Choice::HistogramManagerMT, Exv::GridSurface};}
template<> constexpr Request request_for<hist::HistogramManagerMTFFGridScalableExv>() {return {Choice::HistogramManagerMT, Exv::GridScalable};}

TEST_CASE("HistogramManagerFactory: resolves partial and symmetry preferences") {
    auto exv = settings::exv::exv_method.value;
    auto threads = settings::general::threads;
    auto prefer_partial = settings::internal_state::prefer_partial_manager;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    settings::general::threads = 4;
    settings::hist::weighted_bins = true;

    Molecule plain({Body{SimpleCube::get_atoms()}});

    Molecule symmetric({Body{SimpleCube::get_atoms()}});
    symmetric.get_body(0).symmetry().add(symmetry::type::c2);
    REQUIRE(symmetric.symmetry().has_symmetries());

    SECTION("without a partial preference") {
        settings::internal_state::prefer_partial_manager = false;

        // symmetry-awareness is derived from the molecule, never from a setting
        CHECK(dynamic_cast<hist::HistogramManagerMT<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::SymmetryManagerMT<true>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);
    }

    SECTION("with a partial preference") {
        settings::internal_state::prefer_partial_manager = true;

        CHECK(dynamic_cast<hist::PartialHistogramManagerMT<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::PartialSymmetryManagerMT<true>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);

        // the single-threaded partial manager has no symmetry-aware counterpart, so it upgrades to the MT one
        settings::general::threads = 1;
        CHECK(dynamic_cast<hist::PartialHistogramManager<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::PartialSymmetryManagerMT<true>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);
    }

    SECTION("form factor models get the form factor variant of each kind") {
        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        settings::internal_state::prefer_partial_manager = false;
        CHECK(dynamic_cast<hist::HistogramManagerMTFFExplicit<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::SymmetryManagerMTFF<true>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);

        settings::internal_state::prefer_partial_manager = true;
        CHECK(dynamic_cast<hist::PartialHistogramManagerMTFF<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
        CHECK(dynamic_cast<hist::PartialSymmetryManagerMTFF<true>*>(hist::factory::construct_histogram_manager(&symmetric).get()) != nullptr);

        // the single-threaded partial manager is a weighted reference implementation, so the form factor models use the MT one
        settings::general::threads = 1;
        CHECK(dynamic_cast<hist::PartialHistogramManagerMTFF<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);

        settings::internal_state::prefer_partial_manager = false;
        settings::exv::exv_method = settings::exv::ExvMethod::Average;
        CHECK(dynamic_cast<hist::HistogramManagerMTFFAvg<true>*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
    }

    SECTION("preference is dropped when the excluded volume model has no partial implementation") {
        settings::internal_state::prefer_partial_manager = true;
        settings::exv::exv_method = settings::exv::ExvMethod::Grid;

        // the excluded volume model wins: it changes the result, whereas dropping the partial preference only costs time
        CHECK(dynamic_cast<hist::HistogramManagerMTFFGrid*>(hist::factory::construct_histogram_manager(&plain).get()) != nullptr);
    }

    settings::exv::exv_method = exv;
    settings::general::threads = threads;
    settings::internal_state::prefer_partial_manager = prefer_partial;
}

TEST_CASE("HistogramManagerFactory: creates expected manager") {
    Molecule protein({Body{SimpleCube::get_atoms()}});

    // the surface and scalable grid managers are only used when the excluded volume is fitted
    auto fit_excluded_volume = settings::fit::fit_excluded_volume;
    settings::fit::fit_excluded_volume = true;

    invoke_for_all_histogram_manager_variants(
        []<typename MANAGER>(const Molecule& protein) {
            constexpr auto request = request_for<MANAGER>();
            auto hm = hist::factory::construct_histogram_manager(&protein, request.kind, true, request.exv);
            REQUIRE(dynamic_cast<MANAGER*>(hm.get()) != nullptr);
        },
        []<template<bool> class MANAGER>(const Molecule& protein) {
            constexpr auto request = request_for<MANAGER>();
            auto hm = hist::factory::construct_histogram_manager(&protein, request.kind, false, request.exv);
            REQUIRE(dynamic_cast<MANAGER<false>*>(hm.get()) != nullptr);

            auto hm_w = hist::factory::construct_histogram_manager(&protein, request.kind, true, request.exv);
            REQUIRE(dynamic_cast<MANAGER<true>*>(hm_w.get()) != nullptr);
        },
        protein
    );
    settings::fit::fit_excluded_volume = fit_excluded_volume;
}