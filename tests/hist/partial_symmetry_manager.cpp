#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/Symmetry.h>
#include <data/state/Signaller.h>
#include <rigidbody/BodySplitter.h>
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
    protein.get_body(0).symmetry().get(0).repeat_relation.translate = {0, 1, 0};
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & hydration simultanously
    protein.get_body(0).symmetry().get(0).repeat_relation.translate = {0, -1, 0};
    protein.get_waters().clear();
    protein.signal_modified_hydration_layer();
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & external simultanously
    protein.get_body(0).symmetry().get(0).repeat_relation.translate = {0, 1, 0};
    protein.get_body(0).translate({2, 0, 0});
    phm_res = protein.get_histogram()->get_total_counts();
    p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & internal simultanously
    protein.get_body(0).symmetry().get(0).repeat_relation.translate = {0, -1, 0};
    protein.get_body(0).get_atom(0).weight() = 2;
    protein.get_body(0).get_signaller()->modified_internal();
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

    for (int i = 0; i < 5; ++i) {
        if (modify_symmetry) {
            int body_index = ri(gen) % protein.size_body();
            if (protein.get_body(body_index).size_symmetry() != 0) {
                int symmetry_index = ri(gen) % protein.get_body(body_index).size_symmetry();
                protein.get_body(body_index).symmetry().get(symmetry_index).repeat_relation.translate = {rd(gen), rd(gen), rd(gen)};    
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
                protein.get_body(body_index).get_atom(atom_index).weight() = std::abs(rd(gen));
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
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT;

    // settings::general::threads = 1;

    SECTION("simple") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.get_body(0).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{-1, 0, 0}, {0, 0, 0}}, 1}));

        test(protein);
        test_random(protein);
    }

    SECTION("simple with waters") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}}
        });
        protein.get_body(0).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{-1, 0, 0}, {0, 0, 0}}, 1}));

        test(protein);
        test_random(protein);
    }

    SECTION("2epe") {
        data::Molecule protein({
            Body("tests/files/2epe.pdb"), 
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.get_body(0).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{-1, 0, 0}, {0, 0, 0}}, 1}));
        protein.generate_new_hydration();

        test(protein);
        test_random(protein);
    }

    SECTION("multiple symmetries") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}},
        });
        protein.get_body(0).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{-1, 0, 0}, {0, 0, 0}}, 1}));
        protein.get_body(0).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}}, 1}));
        protein.get_body(1).symmetry().add(symmetry::Symmetry({{{0, 0, 0}, {0, 0, 0}}, {{0,  1, 0}, {0, 0, 0}}, 1}));

        test(protein);
        test_random(protein);
    }

    SECTION("symmetry-heavy") {
        // split the 2epe file into 10 smaller bodies
        auto protein = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {10, 20, 30, 40, 50, 60, 70, 80, 90});
        protein.generate_new_hydration();

        static std::random_device seed;
        static std::mt19937 gen(seed());
        static std::uniform_int_distribution<> ri(1, 10);
        static std::uniform_real_distribution<> rd(-10, 10);
        for (unsigned int i = 0; i < protein.size_body(); ++i) {
            auto& body = protein.get_body(i);
            for (int j = 0; j < ri(gen); ++j) {
                // symmetry with up to 4 repeats
                symmetry::Symmetry sym({{0, 0, 0}, {0, 0, 0}}, {{rd(gen), rd(gen), rd(gen)}, {0, 0, 0}}, (ri(gen) % 4)+1);
                body.symmetry().add(std::move(sym));
            }
        }

        test(protein);
        test_random(protein);
    }
}