#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>

#include <em/ImageStack.h>
#include <plots/all.h>
#include <fitter/SimpleIntensityFitter.h>
#include <fitter/FitReporter.h>
#include <utility/Multiset.h>
#include <utility/Utility.h>
#include <utility/Settings.h>
#include <hist/Histogram.h>

TEST_CASE("consistency_check", "[fitter],[slow],[manual]") {
    unsigned int repeats = 100;

    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    // setting::axes::qmax = 0.4;

    // prepare measured data
    Protein protein("data/A2M/native.pdb");
    SimpleDataset data = protein.get_histogram().calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit_x(Limit(setting::axes::qmin, setting::axes::qmax));
    data.simulate_errors();

    // prepare fit data
    em::ImageStack image("sim/native_23.ccp4");
    auto hist = protein.get_histogram();

    hist::Histogram optimal_vals;
    for (unsigned int i = 0; i < repeats; i++) {
        auto fit = image.fit(hist);
        optimal_vals.p.push_back(fit->get_parameter("cutoff").value);
    }
    optimal_vals.generate_axis();

    plots::PlotHistogram plot(optimal_vals);
    plot.save("figures/test/fitter/consistency_check.pdf");
}