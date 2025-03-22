#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <data/state/Signaller.h>
#include <settings/All.h>

#include <random>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

// Test that the first calculation is correct
TEST_CASE("PartialHistogramManager: initial calculation") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;
    std::vector<std::string> files = {
        "2epe",
        "6lyz",
        "c60",
        "diamond",
        "LAR1-2"
    };
    for (auto f : files) {
        {   // no hydration
            data::Molecule protein("tests/files/" + f + ".pdb");
            protein.clear_hydration();
            auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
            {   // phm
                auto phm = hist::PartialHistogramManager<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, phm, 0, 1e-2));
            }
            {   // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, phm_mt, 0, 1e-2));
            }
            {   // pshm_mt
                auto pshm_mt = hist::PartialSymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, pshm_mt, 0, 1e-2));
            }
        }

        {   // with hydration
            data::Molecule protein("tests/files/" + f + ".pdb");
            protein.generate_new_hydration();
            auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
            {   // phm
                auto phm = hist::PartialHistogramManager<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, phm, 0, 1e-2));
            }
            {   // phm_mt
                auto phm_mt = hist::PartialHistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, phm_mt, 0, 1e-2));
            }
            {   // pshm_mt
                auto pshm_mt = hist::PartialSymmetryManagerMT<true>(&protein).calculate_all()->get_total_counts();
                REQUIRE(compare_hist(p_exp, pshm_mt, 0, 1e-2));
            }
        }
    }
}

auto test = [] (data::Molecule& protein, auto&& phm) {
    // no changes
    auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
    auto phm_res = phm(protein)->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // change hydration
    protein.generate_new_hydration();
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
    phm_res = phm(protein)->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // external change
    protein.get_body(1).translate({1, 1, 1});
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
    phm_res = phm(protein)->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // internal change
    protein.get_body(1).get_atom(0).weight() = 100;
    protein.get_body(1).get_signaller()->modified_internal();
    p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
    phm_res = phm(protein)->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

auto test_random = [] (data::Molecule& protein, auto&& phm) {
    auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
    auto phm_res = phm(protein)->get_total_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    bool modify_external = GENERATE(false, true);
    bool modify_internal = GENERATE(false, true);
    bool modify_hydration = GENERATE(false, true);

    static std::random_device seed;
    static std::mt19937 gen(seed());
    static std::uniform_int_distribution<> ri(0, 100);
    static std::uniform_real_distribution<> rd(-10, 10);

    for (int i = 0; i < 10; ++i) {
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

        auto p_exp = hist::HistogramManagerMT<true>(&protein).calculate_all()->get_total_counts();
        auto phm_res = phm(protein)->get_total_counts();
        REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
    }
};

// Test that subsequent calculations are correct
TEST_CASE("PartialHistogramManager: subsequent calculations") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    data::Molecule protein({
        Body("tests/files/2epe.pdb"), 
        Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
    });

    protein.generate_new_hydration();

    // deterministic
    test(protein, [] (const Molecule& protein) {return hist::PartialHistogramManager<true>(&protein).calculate_all();});
    test(protein, [] (const Molecule& protein) {return hist::PartialHistogramManagerMT<true>(&protein).calculate_all();});
    // test(protein, [] (const Molecule& protein) {return hist::PartialSymmetryManagerMT<true>(&protein).calculate_all();});    

    // random
    test_random(protein, [] (const Molecule& protein) {return hist::PartialHistogramManager<true>(&protein).calculate_all();});
    test_random(protein, [] (const Molecule& protein) {return hist::PartialHistogramManagerMT<true>(&protein).calculate_all();});
    // test_random(protein, [] (const Molecule& protein) {return hist::PartialSymmetryManagerMT<true>(&protein).calculate_all();});    
}