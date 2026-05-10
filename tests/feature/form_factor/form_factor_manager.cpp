#include <catch2/catch_test_macros.hpp>

#include <form_factor/lookup/FormFactorManager.h>
#include <form_factor/FormFactorType.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

#include <hist/hist_test_helper.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>

#include <numeric>
#include <vector>
#include <algorithm>
#include <random>

using namespace ausaxs;
using namespace ausaxs::form_factor;

const std::vector<int>& identity() {
    static std::vector<int> identity;
    if (identity.empty()) {
        identity = std::vector<int>(get_total_ff_count());
        std::iota(identity.begin(), identity.end(), 0);
    }
    return identity;
}

const std::vector<int>& shuffled() {
    static std::vector<int> shuffled;
    if (shuffled.empty()) {
        shuffled = identity();
        std::mt19937 g(std::random_device{}());
        std::shuffle(shuffled.begin()+2, shuffled.end()-1, g);
    }
    return shuffled;
}

template<template<bool, bool> class MANAGER>
void run_comparison(const data::Molecule& protein) {
    manager::detail::use_form_factors(identity());
    auto i1 = MANAGER<false, false>(&protein).calculate_all()->debye_transform();
    auto i2 = MANAGER<true, false>(&protein).calculate_all()->debye_transform();
    auto i3 = MANAGER<false, true>(&protein).calculate_all()->debye_transform();
    auto i4 = MANAGER<true, true>(&protein).calculate_all()->debye_transform();

    manager::detail::use_form_factors(shuffled());
    auto i1s = MANAGER<false, false>(&protein).calculate_all()->debye_transform();
    auto i2s = MANAGER<true, false>(&protein).calculate_all()->debye_transform();
    auto i3s = MANAGER<false, true>(&protein).calculate_all()->debye_transform();
    auto i4s = MANAGER<true, true>(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i1, i1s));
    REQUIRE(compare_hist(i2, i2s));
    REQUIRE(compare_hist(i3, i3s));
    REQUIRE(compare_hist(i4, i4s));
}

template<template<bool> class MANAGER>
void run_comparison(const data::Molecule& protein) {
    manager::detail::use_form_factors(identity());
    auto i1 = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2 = MANAGER<true>(&protein).calculate_all()->debye_transform();

    manager::detail::use_form_factors(shuffled());
    auto i1s = MANAGER<false>(&protein).calculate_all()->debye_transform();
    auto i2s = MANAGER<true>(&protein).calculate_all()->debye_transform();

    REQUIRE(compare_hist(i1, i1s));
    REQUIRE(compare_hist(i2, i2s));
}

TEST_CASE("manager ff set change scattering consistent across all managers") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    data::Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();

    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const data::Molecule& protein) {
            run_comparison<MANAGER>(protein);
        },
        []<template<bool, bool> class MANAGER>(const data::Molecule& protein) {
            run_comparison<MANAGER>(protein);
        },
        protein
    );
    manager::detail::use_form_factors(identity());
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
    manager::detail::use_form_factors(identity());
}

TEST_CASE("form_factor_manager: single form factor") {
    data::Molecule molecule("tests/files/2epe.pdb");
    molecule.clear_hydration();
    for (auto& a : molecule.iterate_atoms()) {a.form_factor_type() = form_factor_t::C;}
    auto I = molecule.get_total_histogram()->debye_transform();

    data::Molecule molecule2("tests/files/2epe.pdb");
    molecule2.clear_hydration();
    manager::detail::use_form_factors({static_cast<int>(form_factor_t::C)});
    auto I2 = molecule2.get_total_histogram()->debye_transform();

    REQUIRE(compare_hist(I, I2));
    manager::detail::use_form_factors(identity());
}