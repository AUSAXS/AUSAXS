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