#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>

#include "hist/intensity_calculator/DistanceHistogram.h"
#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

template<template<bool, bool> class MANAGER>
void run_nongrid_test1(const Molecule& protein) {
    auto h1 = MANAGER<false, false>(&protein).calculate_all();
    REQUIRE(h1->get_d_axis().size() == settings::axes::bin_count);
    auto h2 = MANAGER<true, false>(&protein).calculate_all();
    REQUIRE(h2->get_d_axis().size() == settings::axes::bin_count);
    auto h3 = MANAGER<false, true>(&protein).calculate_all();
    REQUIRE(h3->get_d_axis().size() == settings::axes::bin_count);
    auto h4 = MANAGER<true, true>(&protein).calculate_all();
    REQUIRE(h4->get_d_axis().size() == settings::axes::bin_count);
}

template<template<bool> class MANAGER>
void run_grid_test1(const Molecule& protein) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE(h1->get_d_axis().size() >= 0.95*settings::axes::bin_count);
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(h2->get_d_axis().size() >= 0.95*settings::axes::bin_count);
}

TEST_CASE("Custom bin count: respected by managers") {
    settings::general::verbose = false;
    settings::axes::bin_count = GENERATE(4000, 5000, 6000);

    // for nongrid managers it should be possible to have distances up to max_dist
    double max_dist = settings::axes::bin_width*(settings::axes::bin_count-1);
    std::vector atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::C)
    };
    Molecule protein({Body{atoms}});
    invoke_for_all_nongrid_histogram_manager_variants(
        []<template<bool, bool> class MANAGER>(const Molecule& protein) {
            run_nongrid_test1<MANAGER>(protein);
        },
        protein
    );

    // have to be careful with the expanded exv in grid-based managers
    max_dist = 0.95*settings::axes::bin_width*(settings::axes::bin_count);
    atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::C)
    };
    invoke_for_all_grid_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein) {
            run_grid_test1<MANAGER>(protein);
        },
        protein
    );
}