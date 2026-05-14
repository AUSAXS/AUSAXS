#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_test_macros.hpp>

#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <settings/All.h>
#include <plots/All.h>

using namespace ausaxs;

TEST_CASE("DistanceHistogram::is_highly_ordered") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    SECTION("false") {
        auto false_file = GENERATE(
            "2epe",
            "6lyz",
            "LAR1-2",
            "SASDJQ4"
        );

        SECTION(false_file) {
            auto hist = data::Molecule("tests/files/" + std::string(false_file) + ".pdb").get_histogram();
            REQUIRE_FALSE(hist->is_highly_ordered());
        }
    }

    SECTION("true") {
        auto true_file = GENERATE(
            "6lyz_exv",
            "c60",
            "diamond"        
        );

        SECTION(true_file) {
            auto hist = data::Molecule("tests/files/" + std::string(true_file) + ".pdb").get_histogram(); 
            REQUIRE(hist->is_highly_ordered());
        }
    }
}

TEST_CASE("DistanceHistogram: check relative errors") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = true;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    settings::axes::qmax = constants::axes::q_axis.max;

    // only use small-ish files here
    auto file = GENERATE(
        "2epe",
        "6lyz",
        "LAR1-2"
    );

    INFO("Testing file: " << file);
    auto protein = data::Molecule("tests/files/" + std::string(file) + ".pdb");
    protein.clear_hydration();
    auto hist = protein.get_histogram();
    auto I_exact = hist::exact_debye_transform(protein, constants::axes::q_axis.as_vector());

    {   // default q-range
        auto I_estimated = hist->debye_transform();
        REQUIRE(I_estimated.size() == I_exact.size());
        for (unsigned int i = 0; i < I_estimated.size(); ++i) {
            auto rel_error = std::abs(I_estimated[i] - I_exact[i]) / I_exact[i];
            REQUIRE(rel_error < 0.02);
        }
    }

    {   // custom q-range: debye_transform(q) must give the same result as debye_transform()
        auto q_vec = constants::axes::q_axis.as_vector();
        auto I_default = hist->debye_transform();
        auto I_custom  = hist->debye_transform(q_vec);
        REQUIRE(I_custom.size() == I_default.size());
        for (unsigned int i = 0; i < I_custom.size(); ++i) {
            REQUIRE_THAT(I_custom.y(i), Catch::Matchers::WithinRelMatcher(I_default[i], 1e-6));
        }
    }
}

TEST_CASE("DistanceHistogram: crystals with fine binning") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = true;
    settings::exv::exv_method = settings::exv::ExvMethod::None;
    settings::axes::qmax = constants::axes::q_axis.max;
    settings::axes::bin_width = 0.05;

    auto file = GENERATE(
        "c60",
        "diamond"
    );

    auto protein = data::Molecule("tests/files/" + std::string(file) + ".pdb");
    protein.clear_hydration();
    auto hist = protein.get_histogram();
    auto I_exact = hist::exact_debye_transform(protein, constants::axes::q_axis.as_vector());

    auto I_estimated = hist->debye_transform();
    REQUIRE(I_estimated.size() == I_exact.size());
    double I0 = I_exact[0]; // forward scattering as normalization
    for (unsigned int i = 0; i < I_estimated.size(); ++i) {
        auto abs_err = std::abs(I_estimated[i] - I_exact[i]) / I0;
        REQUIRE(abs_err < 0.01);
    }
}

TEST_CASE("DistanceHistogram: extended q-range") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = true;
    settings::exv::exv_method = settings::exv::ExvMethod::None;

    auto protein = data::Molecule("tests/files/2epe.pdb");
    protein.clear_hydration();
    auto hist = protein.get_histogram();

    std::vector<double> q(200);
    for (unsigned int i = 0; i < q.size(); ++i) {
        q[i] = (i+1)*(10./q.size()); // q values from 0.05 to 10 Å⁻¹
    }
    auto I_dynamic = hist->debye_transform(q);
    auto I_exact   = hist::exact_debye_transform(protein, q);

    REQUIRE(I_dynamic.size() == q.size());
    for (unsigned int i = 0; i < I_dynamic.size(); ++i) {
        auto rel_error = std::abs(I_dynamic.y(i) - I_exact[i]) / I_exact[i];
        REQUIRE(rel_error < 0.03);
    }
}