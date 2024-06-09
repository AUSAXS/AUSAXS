#include "settings/GeneralSettings.h"
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

// TODO: get proper test data for this - the typical EM maps are too large to be included in the repo
// TEST_CASE("fitter: EM_consistency", "[manual]") {
//     unsigned int repeats = 100;

//     settings::molecule::use_effective_charge = false;
//     settings::em::sample_frequency = 2;
//     // setting::axes::qmax = 0.4;

//     // prepare measured data
//     Molecule protein("test/files/2epe.pdb");
//     SimpleDataset data = protein.get_histogram()->debye_transform();
//     data.reduce(settings::fit::N, true);
//     data.limit_x(Limit(settings::axes::qmin, settings::axes::qmax));
//     data.simulate_errors();

//     // prepare fit data
//     em::ImageStack image("test/files/2epe.ccp4");

//     std::vector<double> optimal_vals;
//     for (unsigned int i = 0; i < repeats; i++) {
//         auto hist = protein.get_histogram();
//         auto fit = image.fit(std::move(hist));
//         optimal_vals.push_back(fit->get_parameter("cutoff").value);
//     }
//     hist::Histogram h(optimal_vals);
//     h.generate_axis();

//     plots::PlotHistogram::quick_plot(h, {}, settings::general::output + "test/fitter/consistency_check.png");
// }

TEST_CASE("fitter: consistent_charge_scaling") {
    settings::molecule::use_effective_charge = true;
    settings::general::verbose = false;
    settings::hist::histogram_manager = GENERATE(
        settings::hist::HistogramManagerChoice::HistogramManager, 
        settings::hist::HistogramManagerChoice::PartialHistogramManager,
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
        settings::hist::HistogramManagerChoice::HistogramManagerMT,
        settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg,
        settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg,
        settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit,
        settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid
    );

    auto string_mapper = [] (settings::hist::HistogramManagerChoice choice) -> std::string {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager: return "HistogramManager";
            case settings::hist::HistogramManagerChoice::PartialHistogramManager: return "PartialHistogramManager";
            case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT: return "PartialHistogramManagerMT";
            case settings::hist::HistogramManagerChoice::HistogramManagerMT: return "HistogramManagerMT";
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg: return "HistogramManagerMTFFAvg";
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit: return "HistogramManagerMTFFExplicit";
            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: return "HistogramManagerMTFFGrid";
            default: return "Unknown";
        }
    };

    SECTION("simple, " + string_mapper(settings::hist::histogram_manager)) {
        std::string mfile = "test/files/2epe.dat";
        Molecule protein("test/files/2epe.pdb");
        protein.clear_hydration();
        fitter::HydrationFitter fitter(mfile, protein.get_histogram());
        auto fit1 = fitter.fit_chi2_only();

        protein.update_effective_charge(1.2);
        fitter.set_scattering_hist(protein.get_histogram());
        auto fit2 = fitter.fit_chi2_only();
        REQUIRE(fit1 != fit2);

        protein.update_effective_charge(1.0);
        fitter.set_scattering_hist(protein.get_histogram());
        fit2 = fitter.fit_chi2_only();
        REQUIRE_THAT(fit1, Catch::Matchers::WithinAbs(fit2, 1e-6));
    }
}