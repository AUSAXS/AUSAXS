#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/Symmetry.h>
#include <data/state/Signaller.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"
#include "settings/GridSettings.h"

#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

auto test = [] (data::Molecule& protein) {
    // no changes
    auto p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    auto phm_res = protein.get_histogram()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // add symmetry
    // protein.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});
    // p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    // phm_res = phm(protein)->debye_transform();
    // REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry
    std::cout << "\nMODIFIED SYMMETRY" << std::endl;
    protein.get_body(0).symmetry().get(0).translate = {0, 1, 0};
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & hydration simultanously
    std::cout << "\nMODIFIED SYMMETRY & HYDRATION" << std::endl;
    protein.get_body(0).symmetry().get(0).translate = {0, -1, 0};
    protein.get_waters().clear();
    protein.signal_modified_hydration_layer();
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & external simultanously
    std::cout << "\nMODIFIED SYMMETRY & EXTERNAL" << std::endl;
    protein.get_body(0).symmetry().get(0).translate = {0, 1, 0};
    protein.get_body(0).translate({2, 0, 0});
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

auto test_random = [] (data::Molecule& protein) {
    auto p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    auto phm_res = protein.get_histogram()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    bool modify_symmetry = true;
    bool modify_external = GENERATE(false, true);
    bool modify_internal = GENERATE(false, true);
    bool modify_hydration = GENERATE(false, true);

    static std::random_device seed;
    static std::mt19937 gen(seed());
    static std::uniform_int_distribution<> ri(0, 100);
    static std::uniform_real_distribution<> rd(-10, 10);

    std::cout << "\nRANDOM" << std::endl;
    if (modify_symmetry) {std::cout << "\tMODIFY SYMMETRY" << std::endl;}
    if (modify_external) {std::cout << "\tMODIFY EXTERNAL" << std::endl;}
    if (modify_internal) {std::cout << "\tMODIFY INTERNAL" << std::endl;}
    if (modify_hydration) {std::cout << "\tMODIFY HYDRATION" << std::endl;}
    for (int i = 0; i < 10; ++i) {
        if (modify_symmetry) {
            int body_index = ri(gen) % protein.size_body();
            if (protein.get_body(body_index).size_symmetry() != 0) {
                int symmetry_index = ri(gen) % protein.get_body(body_index).size_symmetry();
                protein.get_body(body_index).symmetry().get(symmetry_index).translate = {rd(gen), rd(gen), rd(gen)};    
            }
        }

        if (modify_external) {
            int body_index = ri(gen) % protein.size_body();
            protein.get_body(body_index).translate({rd(gen), rd(gen), rd(gen)});
        }

        if (modify_internal) {
            int body_index = ri(gen) % protein.size_body();
            int atom_index = ri(gen) % protein.get_body(body_index).size_atom();
            for (int j = 0; j < ri(gen); ++j) {
                protein.get_body(body_index).get_atom(atom_index).weight() = rd(gen);
            }
            protein.get_body(body_index).get_signaller()->modified_internal();
        }

        if (modify_hydration) {
            protein.generate_new_hydration();
        }

        auto p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
        auto phm_res = protein.get_histogram()->get_total_counts();
        REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
    }
};

// Test that subsequent calculations are correct
TEST_CASE("PartialSymmetryManagerMT: subsequent calculations") {
    settings::general::threads = 1;
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT;

    SECTION("simple") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.get_body(0).symmetry().add({{-1, 0, 0}});

        // SECTION("deterministic") {
        //     test(protein);
        // }

        SECTION("random") {
            test_random(protein);
        }
    }

    // SECTION("simple with waters") {
    //     data::Molecule protein({
    //         Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}}, 
    //         Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}}
    //     });
    //     protein.get_body(0).symmetry().add({{-1, 0, 0}});

    //     SECTION("deterministic") {
    //         test(protein);
    //     }

    //     SECTION("random") {
    //         test_random(protein);
    //     }
    // }

    // SECTION("2epe") {
    //     data::Molecule protein({
    //         Body("tests/files/2epe.pdb"), 
    //         Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
    //     });
    //     protein.get_body(0).symmetry().add({Vector3<double>(-1, 0, 0)});
    //     protein.generate_new_hydration();

    //     SECTION("deterministic") {
    //         test(protein);
    //     }
        
    //     SECTION("random") {
    //         test_random(protein);
    //     }
    // }
}