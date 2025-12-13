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
#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <settings/All.h>

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
    REQUIRE(h1->get_d_axis().size() > 0.95*settings::axes::bin_count);
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(h2->get_d_axis().size() > 0.95*settings::axes::bin_count);
}
TEST_CASE("Custom bin count: respected by managers") {
    settings::general::verbose = false;
    settings::axes::bin_count = GENERATE(4000, 5000, 6000);

    // for nongrid managers it should be possible to have distances up to max_dist
    settings::grid::min_exv_radius = 0;
    double max_dist = settings::axes::bin_width*(settings::axes::bin_count-1);
    std::vector atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::H),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::H)
    };
    Molecule protein({Body{atoms}});
    invoke_for_all_nongrid_histogram_manager_variants(
        []<template<bool, bool> class MANAGER>(const Molecule& protein) {
            run_nongrid_test1<MANAGER>(protein);
        },
        protein
    );

    // have to be careful with the expanded exv in grid-based managers
    max_dist = 0.95*settings::axes::bin_width*(settings::axes::bin_count-1);
    atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::C)
    };
    protein = Molecule({Body{atoms}});
    invoke_for_all_grid_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein) {
            run_grid_test1<MANAGER>(protein);
        },
        protein
    );
}

template<template<bool, bool> class MANAGER>
void run_test2(const Molecule& protein) {
    auto h1 = MANAGER<false, false>(&protein).calculate_all();
    REQUIRE_THAT(h1->get_d_axis()[1] - h1->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
    auto h2 = MANAGER<true, false>(&protein).calculate_all();
    REQUIRE_THAT(h2->get_d_axis()[1] - h2->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
    auto h3 = MANAGER<false, true>(&protein).calculate_all();
    REQUIRE_THAT(h3->get_d_axis()[1] - h3->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
    auto h4 = MANAGER<true, true>(&protein).calculate_all();
    REQUIRE_THAT(h4->get_d_axis()[1] - h4->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
}

template<template<bool> class MANAGER>
void run_test2(const Molecule& protein) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE_THAT(h1->get_d_axis()[1] - h1->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE_THAT(h2->get_d_axis()[1] - h2->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
}
TEST_CASE("Custom bin width: respected by managers") {
    settings::general::verbose = false;
    settings::axes::bin_width = GENERATE(0.1, 0.05, 0.02);

    Molecule protein({Body{SimpleCube::get_atoms()}});
    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein) {
            run_test2<MANAGER>(protein);
        },
        []<template<bool, bool> class MANAGER>(const Molecule& protein) {
            run_test2<MANAGER>(protein);
        },
        protein
    );
}

template<template<bool> class MANAGER>
void run_test3(const Molecule& protein, const auto& target) {
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h2.get()), target));
}
template<template<bool, bool> class MANAGER>
void run_test3(const Molecule& protein, const auto& target) {
    auto h3 = MANAGER<false, true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h3.get()), target));
    auto h4 = MANAGER<true, true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h4.get()), target));
}
TEST_CASE("Custom bin width: varying widths agree with analytical result") {
    settings::general::verbose = false;

    static auto calc_exp = [] (double width) {
        std::vector<double> res(std::round(3.5/width)+1);
        res[0] = 8;
        res[std::round(2/width)] += 8*3;
        res[std::round(std::sqrt(8)/width)] += 8*3;
        res[std::round(std::sqrt(12)/width)] += 8*1;
        return res;
    };

    settings::axes::bin_width = GENERATE(0.1, 0.05, 0.02);
    Molecule protein({Body{SimpleCube::get_atoms()}});
    set_unity_charge(protein);

    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein, const auto& target) {
            run_test3<MANAGER>(protein, target);
        },
        []<template<bool, bool> class MANAGER>(const Molecule& protein, const auto& target) {
            run_test3<MANAGER>(protein, target);
        },
        protein, calc_exp(settings::axes::bin_width)
    );
}

template<template<bool> class MANAGER>
void run_test4(const Molecule& protein) {
    auto iq = MANAGER<false>(&protein).calculate_all()->debye_transform();
    settings::axes::bin_width = constants::axes::d_axis.width();
    settings::axes::bin_count = 8000/settings::axes::bin_width;
    auto iq2 = MANAGER<true>(&protein).calculate_all()->debye_transform();
    REQUIRE(compare_hist(iq, iq2, 1e-6, 0.005));
}
template<template<bool, bool> class MANAGER>
void run_test4(const Molecule& protein) {
    settings::axes::bin_width = constants::axes::d_axis.width();
    settings::axes::bin_count = 8000/settings::axes::bin_width;

    auto iq = MANAGER<false, false>(&protein).calculate_all()->debye_transform();
    auto iq2 = MANAGER<false, true>(&protein).calculate_all()->debye_transform();
    REQUIRE(compare_hist(iq, iq2, 1e-6, 0.005));

    iq = MANAGER<true, false>(&protein).calculate_all()->debye_transform();
    iq2 = MANAGER<true, true>(&protein).calculate_all()->debye_transform();
    REQUIRE(compare_hist(iq, iq2, 1e-6, 0.005));
}
TEST_CASE("Custom bin width: fixed and variable widths agree") {
    settings::general::verbose = false;
    Molecule protein("tests/files/2epe.pdb");
    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein) {
            run_test4<MANAGER>(protein);
        },
        []<template<bool, bool> class MANAGER>(const Molecule& protein) {
            run_test4<MANAGER>(protein);
        },
        protein
    );
}

auto avg_deviation = [] (const std::vector<double>& a, const std::vector<double>& b) {
    double total_dev = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        total_dev += std::abs(a[i]-b[i])/b[i];
    }
    return total_dev/a.size();
};
template<template<bool, bool> class MANAGER>
void run_test5(const Molecule& protein, const std::vector<double>& exact) {
    settings::axes::bin_width = 0.5;
    settings::axes::bin_count = 8000/settings::axes::bin_width;
    auto target_dev = avg_deviation(
        MANAGER<true, false>(&protein).calculate_all()->debye_transform().get_counts(),
        exact
    );
    for (auto width : {0.25, 0.15, 0.1}) {
        settings::axes::bin_width = width;
        settings::axes::bin_count = 8000/settings::axes::bin_width;
        auto iq = MANAGER<true, true>(&protein).calculate_all()->debye_transform().get_counts();
        REQUIRE(avg_deviation(iq, exact) <= target_dev*1.001); // allow numerical noise
    }
}
TEST_CASE("Custom bin width: smaller widths increase accuracy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    Molecule protein("tests/files/c60.pdb");
    auto exact = hist::exact_debye_transform(protein, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector());
    invoke_for_all_nongrid_histogram_manager_variants(
        []<template<bool, bool> class MANAGER>(const Molecule& protein, const std::vector<double>& exact) {
            run_test5<MANAGER>(protein, exact);
        },
        protein, exact
    );
}