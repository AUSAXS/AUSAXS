#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/state/StateManager.h>
#include <data/state/Signaller.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
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
#include <io/ExistingFile.h>
#include <settings/All.h>
#include <constants/Constants.h>
#include <utility/Utility.h>

#include <hist/hist_test_helper.h>

using namespace ausaxs;
using namespace ausaxs::data;

struct analytical_histogram {
    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    std::vector<double> calc_exp() {
        auto width = constants::axes::d_axis.width();
        std::vector<double> res(std::round(3.5/width)+1);
        res[0] = 8;
        res[std::round(2/width)] += 8*3;
        res[std::round(std::sqrt(8)/width)] += 8*3;
        res[std::round(std::sqrt(12)/width)] += 8*1;
        return res;
    }

    std::vector<double> p_exp = calc_exp();
};

template<template<bool> class MANAGER>
void run_test1(const Molecule& protein, const auto& target) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h1.get()), target));
    
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h2.get()), target));
}

template<template<bool, bool> class MANAGER>
void run_test1(const Molecule& protein, const auto& target) {
    auto h1 = MANAGER<false, false>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h1.get()), target));

    auto h2 = MANAGER<false, true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h2.get()), target));

    auto h3 = MANAGER<true, false>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h3.get()), target));

    auto h4 = MANAGER<true, true>(&protein).calculate_all();
    REQUIRE(compare_hist(get_raw_counts(h4.get()), target));
}

TEST_CASE_METHOD(analytical_histogram, "HistogramManager::calculate_all") {
    settings::molecule::implicit_hydrogens = false;
    settings::general::verbose = false;

    SECTION("analytical") {
        SECTION("atoms only") {
            // the following just describes the eight corners of a cube centered at origo
            std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
            std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
            std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
            std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
            std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
            Molecule protein(a);
            protein.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManager);
            set_unity_charge(protein);

            invoke_for_all_histogram_manager_variants(
                []<template<bool> class MANAGER>(const Molecule& protein, const auto& target) {
                    run_test1<MANAGER>(protein, target);
                },
                []<template<bool, bool> class MANAGER>(const Molecule& protein, const auto& target) {
                    run_test1<MANAGER>(protein, target);
                },
                protein, p_exp
            );
        }

        SECTION("waters only") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            std::vector<AtomFF> a = {};
            std::vector<Water> w = {
                Water({-1, -1, -1}), Water({-1, 1, -1}), 
                Water({ 1, -1, -1}), Water({ 1, 1, -1}), 
                Water({-1, -1,  1}), Water({-1, 1,  1}),
                Water({ 1, -1,  1}), Water({ 1, 1,  1})
            };
            Molecule protein({Body{a, w}});
            set_unity_charge(protein);

            invoke_for_all_nongrid_histogram_manager_variants(
                []<template<bool, bool> class MANAGER>(const Molecule& protein, const auto& target) {
                    run_test1<MANAGER>(protein, target);
                },
                protein, p_exp
            );
            // grid-based doesn't make sense for a water-only system and will throw an exception - excluded volume is based on atomic volumes
            // hm_mt_ff_grid 
            // hm_mt_ff_grid_surface
            // hm_mt_ff_grid_scalable_exv
        }

        SECTION("both waters and atoms") {
            // the following just describes the eight corners of a cube centered at origo, with an additional atom at the very middle
            std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
            std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
            std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
            std::vector<Water> w = {Water({1, -1,  1}), Water({1, 1, 1})};
            std::vector<Body> a = {Body(b1, w), Body(b2), Body(b3)};
            Molecule protein(a);
            set_unity_charge(protein);

            invoke_for_all_histogram_manager_variants(
                []<template<bool> class MANAGER>(const Molecule& protein, const auto& target) {
                    run_test1<MANAGER>(protein, target);
                },
                []<template<bool, bool> class MANAGER>(const Molecule& protein, const auto& target) {
                    run_test1<MANAGER>(protein, target);
                },
                protein, p_exp                
            );
        }
    }
}

template<template<bool> class MANAGER>
void run_test2(const Molecule& protein, const auto& target) {
    auto h1 = MANAGER<false>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h1->get_weighted_counts(), target));
    
    auto h2 = MANAGER<true>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h2->get_weighted_counts(), target));
}

template<template<bool, bool> class MANAGER>
void run_test2(const Molecule& protein, const auto& target) {
    auto h1 = MANAGER<false, false>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h1->get_weighted_counts(), target));

    auto h2 = MANAGER<false, true>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h2->get_weighted_counts(), target));

    auto h3 = MANAGER<true, false>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h3->get_weighted_counts(), target));

    auto h4 = MANAGER<true, true>(&protein).calculate_all();
    REQUIRE(compare_hist_approx(h4->get_weighted_counts(), target));
}

TEST_CASE("HistogramManager::calculate_all real data") {
    settings::molecule::implicit_hydrogens = false;
    settings::general::verbose = false;

    Molecule protein("tests/files/2epe.pdb");
    protein.generate_new_hydration();
    auto p_exp = protein.get_histogram();

    invoke_for_all_histogram_manager_variants(
        []<template<bool> class MANAGER>(const Molecule& protein, const auto& target) {
            run_test2<MANAGER>(protein, target);
        },
        []<template<bool, bool> class MANAGER>(const Molecule& protein, const auto& target) {
            run_test2<MANAGER>(protein, target);
        },
        protein, p_exp->get_weighted_counts()
    );
}

TEST_CASE("PartialHistogramManager::get_probe") {
    settings::general::verbose = false;
    Molecule protein("tests/files/2epe.pdb");
    auto phm = hist::PartialHistogramManager<false, false>(&protein);
    auto sm = phm.get_state_manager();

    // check the signalling object is correct
    CHECK(phm.get_probe(0) == sm->get_probe(0)); 

    // check that it links to the state manager
    sm->reset_to_false();
    phm.get_probe(0)->modified_external();
    CHECK(sm->is_externally_modified(0));
}

TEST_CASE("PartialHistogramManager::signal_modified_hydration_layer") {
    settings::general::verbose = false;
    Molecule protein("tests/files/2epe.pdb");
    auto phm = hist::PartialHistogramManager<false, false>(&protein);
    auto sm = phm.get_state_manager();
    sm->reset_to_false();
    phm.signal_modified_hydration_layer();
    CHECK(sm->is_modified_hydration());
}