#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
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

auto make_unique_cyclic_sym = [] (
    const Vector3<double>& initial_relation, const Vector3<double>& per_step_translation, const Vector3<double>& axis, double angle, int repetitions = 1
) {
    return std::make_unique<symmetry::CyclicSymmetry>(initial_relation, per_step_translation, axis, angle, repetitions);
};

auto make_unique_point_sym = [] (const Vector3<double>& translation, const Vector3<double>& rotation) {
    return std::make_unique<symmetry::PointSymmetry>(translation, rotation);
};

auto test = [] (data::Molecule& protein) {
    auto cast = [&protein](int body_idx, int sym_idx) {
        return static_cast<symmetry::CyclicSymmetry*>(protein.get_body(body_idx).symmetry().get(sym_idx));
    };

    // no changes
    auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = protein.get_histogram()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // add symmetry
    // protein.get_body(0).symmetry().add({(1, 0, 0)});
    // p_exp = hist::SymmetryManagerMT<true>(&protein).calculate_all()->debye_transform();
    // phm_res = phm(protein)->debye_transform();
    // REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry
    cast(0, 0)->_repeat_relation.translation = {0, 1, 0};
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & hydration simultanously
    cast(0, 0)->_repeat_relation.translation = {0, -1, 0};
    protein.get_waters().clear();
    protein.signal_modified_hydration_layer();
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & external simultanously
    cast(0, 0)->_repeat_relation.translation = {0, 1, 0};
    protein.get_body(0).translate({2, 0, 0});
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify symmetry & internal simultanously
    cast(0, 0)->_repeat_relation.translation = {0, -1, 0};
    protein.get_body(0).get_atom(0).weight() = 2;
    protein.get_body(0).get_signaller()->modified_internal();
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

auto test_random = [] (data::Molecule& protein) {
    auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = protein.get_histogram()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    bool modify_symmetry = true;
    bool modify_external = GENERATE(false, true);
    bool modify_internal = GENERATE(false, true);
    bool modify_hydration = GENERATE(false, true);

    static std::random_device seed;
    static std::mt19937 gen(seed());
    static std::uniform_int_distribution<> ri(0, 100);
    static std::uniform_real_distribution<> rd(-10, 10);

    for (int i = 0; i < 2; ++i) {
        if (modify_symmetry) {
            int body_index = ri(gen) % protein.size_body();
            if (protein.get_body(body_index).size_symmetry() != 0) {
                int symmetry_index = ri(gen) % protein.get_body(body_index).size_symmetry();
                static_cast<symmetry::CyclicSymmetry*>(protein.get_body(body_index).symmetry().get(symmetry_index))
                    ->_repeat_relation.translation = {rd(gen), rd(gen), rd(gen)}
                ;
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

        auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
        auto phm_res = protein.get_histogram()->get_weighted_counts();
        REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
    }
};

// Test that subsequent calculations are correct
TEST_CASE("PartialSymmetryManagerMT: subsequent calculations") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;

    SECTION("simple") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        // initial_relation = {0,0,0}, per-step translation = {-1,0,0}, no rotation, 1 repeat
        protein.get_body(0).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, {-1, 0, 0}, {0, 0, 1}, 0, 1));

        test(protein);
        test_random(protein);
    }

    SECTION("simple with waters") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, {-1, 0, 0}, {0, 0, 1}, 0.0, 1));

        test(protein);
        test_random(protein);
    }

    SECTION("2epe") {
        data::Molecule protein({
            Body("tests/files/2epe.pdb"), 
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, {-1, 0, 0}, {0, 0, 1}, 0.0, 1));
        protein.generate_new_hydration();

        test(protein);
        test_random(protein);
    }

    SECTION("multiple symmetries") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}}, 
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}},
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, {-1, 0, 0}, {0, 0, 1}, 0.0, 1));
        protein.get_body(0).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, { 0,-1, 0}, {0, 0, 1}, 0.0, 1));
        protein.get_body(1).symmetry().add(make_unique_cyclic_sym({0, 0, 0}, { 0, 1, 0}, {0, 0, 1}, 0.0, 1));

        test(protein);
        test_random(protein);
    }

    SECTION("symmetry-heavy") {
        // split the 2epe file into 10 smaller bodies
        auto protein = rigidbody::BodySplitter::split("tests/files/2epe.pdb", {10, 20, 30, 40, 50, 60, 70, 80, 90});
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.generate_new_hydration();

        static std::random_device seed;
        static std::mt19937 gen(seed());
        static std::uniform_int_distribution<> ri(1, 10);
        static std::uniform_real_distribution<> rd(-10, 10);
        for (unsigned int i = 0; i < protein.size_body(); ++i) {
            auto& body = protein.get_body(i);
            for (int j = 0; j < ri(gen); ++j) {
                // symmetry with up to 4 repeats
                body.symmetry().add(make_unique_cyclic_sym({0, 0, 0}, {rd(gen), rd(gen), rd(gen)}, {0, 0, 1}, 0.0, (ri(gen) % 4)+1));
            }
        }

        test(protein);
        test_random(protein);
    }
}

auto test_point = [] (data::Molecule& protein) {
    auto cast = [&protein](int body_idx, int sym_idx) {
        return static_cast<symmetry::PointSymmetry*>(protein.get_body(body_idx).symmetry().get(sym_idx));
    };

    // no changes
    auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = protein.get_histogram()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify translation
    cast(0, 0)->translation = {0, 1, 0};
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify translation & hydration simultaneously
    cast(0, 0)->translation = {0, -1, 0};
    protein.get_waters().clear();
    protein.signal_modified_hydration_layer();
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify translation & external simultaneously
    cast(0, 0)->translation = {0, 1, 0};
    protein.get_body(0).translate({2, 0, 0});
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    // modify translation & internal simultaneously
    cast(0, 0)->translation = {0, -1, 0};
    protein.get_body(0).get_atom(0).weight() = 2;
    protein.get_body(0).get_signaller()->modified_internal();
    phm_res = protein.get_histogram()->get_weighted_counts();
    p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));
};

auto test_point_random = [] (data::Molecule& protein) {
    auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
    auto phm_res = protein.get_histogram()->get_weighted_counts();
    REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

    bool modify_symmetry = true;
    bool modify_external = GENERATE(false, true);
    bool modify_internal = GENERATE(false, true);
    bool modify_hydration = GENERATE(false, true);

    static std::random_device seed;
    static std::mt19937 gen(seed());
    static std::uniform_int_distribution<> ri(0, 100);
    static std::uniform_real_distribution<> rd(-10, 10);

    for (int i = 0; i < 2; ++i) {
        if (modify_symmetry) {
            int body_index = ri(gen) % protein.size_body();
            if (protein.get_body(body_index).size_symmetry() != 0) {
                int symmetry_index = ri(gen) % protein.get_body(body_index).size_symmetry();
                static_cast<symmetry::PointSymmetry*>(protein.get_body(body_index).symmetry().get(symmetry_index))
                    ->translation = {rd(gen), rd(gen), rd(gen)}
                ;
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

        auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
        auto phm_res = protein.get_histogram()->get_weighted_counts();
        REQUIRE(compare_hist_approx(p_exp, phm_res, 0, 1e-2));
    }
};

// Test that PartialSymmetryManagerMT correctly handles PointSymmetry updates
TEST_CASE("PartialSymmetryManagerMT: PointSymmetry subsequent calculations") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::grid::min_bins = 100;

    SECTION("simple") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}},
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_point_sym({-1, 0, 0}, {0, 0, 0}));

        test_point(protein);
        test_point_random(protein);
    }

    SECTION("simple with waters") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}},
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_point_sym({-1, 0, 0}, {0, 0, 0}));

        test_point(protein);
        test_point_random(protein);
    }

    SECTION("2epe") {
        data::Molecule protein({
            Body("tests/files/2epe.pdb"),
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}}
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_point_sym({-1, 0, 0}, {0, 1, 0}));
        protein.generate_new_hydration();

        test_point(protein);
        test_point_random(protein);
    }

    SECTION("multiple PointSymmetries") {
        data::Molecule protein({
            Body{std::vector{AtomFF({0, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({0, 0, 1})}},
            Body{std::vector{AtomFF({1, 0, 0}, form_factor::form_factor_t::C)}, std::vector{Water({1, 0, 1})}},
        });
        protein.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT);
        protein.get_body(0).symmetry().add(make_unique_point_sym({-1, 0, 0}, {0, 0, 1}));
        protein.get_body(0).symmetry().add(make_unique_point_sym({ 0,-1, 0}, {1, 0, 0}));
        protein.get_body(1).symmetry().add(make_unique_point_sym({ 0, 1, 0}, {0, 1, 0}));

        // test_point only works for single PointSymmetry on body 0; do the full random check instead
        auto p_exp = hist::SymmetryManagerMT<true, false>(&protein).calculate_all()->get_weighted_counts();
        auto phm_res = protein.get_histogram()->get_weighted_counts();
        REQUIRE(compare_hist(p_exp, phm_res, 0, 1e-2));

        test_point_random(protein);
    }
}