#include "plots/PlotDataset.h"
#include "plots/PlotOptions.h"
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_test_macros.hpp>

#include <data/Molecule.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/All.h>
#include <plots/All.h>

TEST_CASE("DistanceHistogram::is_highly_ordered") {
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;

    SECTION("false") {
        auto false_file = GENERATE(
            "2epe",
            "6lyz",
            "LAR1-2",
            "SASDJQ4"
        );

        SECTION(false_file) {
            auto hist = data::Molecule(settings::general::output + "test/files/" + std::string(false_file) + ".pdb").get_histogram();
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
            auto hist = data::Molecule(settings::general::output + "test/files/" + std::string(true_file) + ".pdb").get_histogram(); 
            REQUIRE(hist->is_highly_ordered());
        }
    }
}

TEST_CASE("DistanceHistogram: check relative errors") {
    settings::general::verbose = false;
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    settings::hist::weighted_bins = true;

    // only use small-ish files here
    auto file = GENERATE(
        "2epe",
        "6lyz",
        "LAR1-2",
        "diamond",
        "c60"
    );

    SECTION(file) {
        std::cout << file << std::endl;
        auto protein = data::Molecule("test/files/" + std::string(file) + ".pdb");
        protein.clear_hydration();
        auto hist = protein.get_histogram();
        auto I_estimated = hist->debye_transform();
        auto I_exact = protein.debye_transform();

        plots::PlotDataset plot;
        plot.plot(I_estimated.as_dataset(), plots::PlotOptions({{"color", style::color::black}, {"logx", true}, {"logy", true}}));
        plot.plot(SimpleDataset(std::vector<double>(constants::axes::q_vals.begin(), constants::axes::q_vals.end()), I_exact), plots::PlotOptions({{"color", style::color::red}}));
        plot.save(settings::general::output + "test/hist/distance_histogram/relative_errors.png");

        for (unsigned int i = 0; i < I_estimated.size(); ++i) {
            auto rel_error = std::abs(I_estimated[i] - I_exact[i]) / I_exact[i];
            REQUIRE(rel_error < 0.01);
        }
    }
}