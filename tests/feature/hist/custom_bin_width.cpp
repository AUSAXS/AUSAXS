#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::data;

template<template<bool> class MANAGER>
static void run_nongrid_test1(const Molecule& protein, std::size_t expected_bins) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE(h1->get_d_axis().size() == expected_bins);
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(h2->get_d_axis().size() == expected_bins);
}

template<typename MANAGER>
static void run_grid_test1(const Molecule& protein, std::size_t min_bins) {
    auto h = MANAGER(&protein).calculate_all();
    REQUIRE(h->get_d_axis().size() >= min_bins);
}
TEST_CASE("Deduced bin count: axis covers the structure") {
    settings::general::verbose = false;
    double max_dist = GENERATE(250., 500., 1000.);
    auto expected_bins = static_cast<std::size_t>(std::round(max_dist/settings::axes::bin_width) + 1);

    settings::grid::min_exv_radius = 0;
    std::vector atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::H),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::H)
    };
    Molecule protein({Body{atoms}});
    invoke_for_all_nongrid_histogram_manager_variants(
        [expected_bins]<template<bool> class MANAGER>(const Molecule& protein) {
            run_nongrid_test1<MANAGER>(protein, expected_bins);
        },
        protein
    );

    // the grid-based managers histogram an expanded excluded volume, so they can only be bounded below
    atoms = {
        data::AtomFF({0, 0, 0}, form_factor::form_factor_t::C),
        data::AtomFF({max_dist, 0, 0}, form_factor::form_factor_t::C)
    };
    protein = Molecule({Body{atoms}});
    invoke_for_all_grid_histogram_manager_variants(
        [expected_bins]<typename MANAGER>(const Molecule& protein) {
            run_grid_test1<MANAGER>(protein, expected_bins);
        },
        protein
    );
}

template<template<bool> class MANAGER>
static void run_test2(const Molecule& protein) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE_THAT(h1->get_d_axis()[1] - h1->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE_THAT(h2->get_d_axis()[1] - h2->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
}

template<typename MANAGER>
static void run_test2(const Molecule& protein) {
    auto h = MANAGER(&protein).calculate_all();
    REQUIRE_THAT(h->get_d_axis()[1] - h->get_d_axis()[0], Catch::Matchers::WithinAbs(settings::axes::bin_width, 1e-9));
}
TEST_CASE("Custom bin width: respected by managers") {
    settings::general::verbose = false;
    settings::axes::bin_width = GENERATE(0.1, 0.05, 0.02);

    Molecule protein({Body{SimpleCube::get_atoms()}});
    invoke_for_all_histogram_manager_variants(
        []<typename MANAGER>(const Molecule& protein) {
            run_test2<MANAGER>(protein);
        },
        []<template<bool> class MANAGER>(const Molecule& protein) {
            run_test2<MANAGER>(protein);
        },
        protein
    );
}

template<typename MANAGER>
static void run_test3(const Molecule& protein, const auto& target) {
    auto h = MANAGER(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h.get()), target));
}
template<template<bool> class MANAGER>
static void run_test3(const Molecule& protein, const auto& target) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h1.get()), target));
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h2.get()), target));
}
TEST_CASE("Custom bin width: varying widths agree with analytical result") {
    settings::general::verbose = false;

    static auto calc_exp = [] (double width) {
        std::vector<double> res(static_cast<std::size_t>(std::round(3.5/width)+1));
        res[0] = 8;
        res[static_cast<std::size_t>(std::round(2/width))] += 8*3;
        res[static_cast<std::size_t>(std::round(std::sqrt(8)/width))] += 8*3;
        res[static_cast<std::size_t>(std::round(std::sqrt(12)/width))] += 8*1;
        return res;
    };

    settings::axes::bin_width = GENERATE(0.1, 0.05, 0.02);
    Molecule protein({Body{SimpleCube::get_atoms()}});
    set_unity_charge(protein);

    invoke_for_all_histogram_manager_variants(
        []<typename MANAGER>(const Molecule& protein, const auto& target) {
            run_test3<MANAGER>(protein, target);
        },
        []<template<bool> class MANAGER>(const Molecule& protein, const auto& target) {
            run_test3<MANAGER>(protein, target);
        },
        protein, calc_exp(settings::axes::bin_width)
    );
}

static auto avg_deviation = [] (const std::vector<double>& a, const std::vector<double>& b) {
    double total_dev = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        total_dev += std::abs(a[i]-b[i])/b[i];
    }
    return total_dev/static_cast<double>(a.size());
};
template<template<bool> class MANAGER>
static void run_test5(const Molecule& protein, const std::vector<double>& exact) {
    settings::axes::bin_width = 0.5;
    auto target_dev = avg_deviation(
        MANAGER<true>(&protein).calculate_all()->debye_transform().get_counts(),
        exact
    );
    for (auto width : {0.25, 0.15, 0.1}) {
        settings::axes::bin_width = width;
        auto iq = MANAGER<true>(&protein).calculate_all()->debye_transform().get_counts();
        REQUIRE(avg_deviation(iq, exact) <= target_dev*1.001); // allow numerical noise
    }
}
TEST_CASE("Custom bin width: smaller widths increase accuracy") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    Molecule protein("tests/files/c60.pdb");
    auto exact = hist::exact_debye_transform(protein, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector());
    invoke_for_all_nongrid_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein, const std::vector<double>& exact) {
            run_test5<MANAGER>(protein, exact);
        },
        protein, exact
    );
}