#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <em/ImageStack.h>
#include <plots/All.h>
#include <fitter/FitReporter.h>
#include <utility/Utility.h>
#include <settings/All.h>
#include <hist/Histogram.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Water.h>
#include <em/detail/ExtendedLandscape.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <mini/detail/FittedParameter.h>

using namespace data;

TEST_CASE("consistency_check", "[slow],[manual]") {
    unsigned int repeats = 100;

    settings::molecule::use_effective_charge = false;
    settings::em::sample_frequency = 2;
    // setting::axes::qmax = 0.4;

    // prepare measured data
    Molecule protein("data/A2M/native.pdb");
    SimpleDataset data = protein.get_histogram()->debye_transform();
    data.reduce(settings::fit::N, true);
    data.limit_x(Limit(settings::axes::qmin, settings::axes::qmax));
    data.simulate_errors();

    // prepare fit data
    em::ImageStack image("sim/native_23.ccp4");

    std::vector<double> optimal_vals;
    for (unsigned int i = 0; i < repeats; i++) {
        auto hist = protein.get_histogram();
        auto fit = image.fit(std::move(hist));
        optimal_vals.push_back(fit->get_parameter("cutoff").value);
    }
    hist::Histogram h(optimal_vals);
    h.generate_axis();

    plots::PlotHistogram::quick_plot(h, {}, settings::general::output + "/test/fitter/consistency_check.png");
}

TEST_CASE("excluded_volume") {
    settings::molecule::use_effective_charge = true;
    settings::hist::histogram_manager = GENERATE(
        settings::hist::HistogramManagerChoice::HistogramManager, 
        settings::hist::HistogramManagerChoice::PartialHistogramManager,
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
        settings::hist::HistogramManagerChoice::HistogramManagerMT
    );

    //! broken because bodies are not registered to the signaller at time of call?
    SECTION("simple, " + std::to_string(static_cast<int>(settings::hist::histogram_manager))) {
        std::string mfile = "test/files/2epe.dat";
        Molecule protein("test/files/2epe.pdb");

        fitter::HydrationFitter fitter(mfile, protein.get_histogram());
        auto fit1 = fitter.fit();

        protein.update_effective_charge(1.2);
        fitter.set_scattering_hist(protein.get_histogram());
        auto fit2 = fitter.fit();
        REQUIRE(fit1->fval != fit2->fval);

        protein.update_effective_charge(1.0);
        fitter.set_scattering_hist(protein.get_histogram());
        fit2 = fitter.fit();
        REQUIRE_THAT(fit1->fval, Catch::Matchers::WithinAbs(fit2->fval, 1e-6));
    }

}