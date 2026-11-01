#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Molecule.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/scattering_model/ScatteringModelSimple.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;

// the plain managers subtract the displaced charge from each atomic weight before binning, while the partial managers apply it afterwards
// to the histogram of each form factor type; both must describe the same model. Every atom of these structures is fully occupied, since
// the form factor histograms carry no per-atom weights
TEST_CASE("make_simple_model: matches the plain managers") {
    settings::general::verbose = false;
    settings::exv::exv_method = settings::exv::ExvMethod::Simple;
    bool weighted_bins = GENERATE(false, true);
    const auto* file = GENERATE("tests/files/2epe.pdb", "tests/files/SASDJG5.pdb");
    auto partial = GENERATE(
        settings::hist::HistogramManagerChoice::PartialHistogramManagerMT,
        settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT
    );

    Molecule protein(file);
    protein.generate_new_hydration();
    auto plain = hist::factory::construct_histogram_manager(&protein, settings::hist::HistogramManagerChoice::HistogramManagerMT, weighted_bins)->calculate_all();
    auto ff = hist::factory::construct_histogram_manager(&protein, partial, weighted_bins)->calculate_all();
    REQUIRE(dynamic_cast<hist::ScatteringModelSimple*>(plain.get()) != nullptr);
    REQUIRE(dynamic_cast<hist::ScatteringModelSimple*>(ff.get()) != nullptr);

    // with unweighted bins, the only difference is the summation order. A weighted bin sits at its mean distance, which the plain managers
    // weigh by the shifted pair weights and the form factor managers by the count of each type pair. The excluded volume largely cancels
    // the atomic term, so this moves the sum by up to 0.15% near the minima
    double tolerance = weighted_bins ? 3e-3 : 1e-4;
    for (double cw : {1., 2.}) {
        plain->apply_water_scaling_factor(cw);
        ff->apply_water_scaling_factor(cw);
        auto I_plain = plain->debye_transform().get_intensity();
        auto I_ff = ff->debye_transform().get_intensity();
        REQUIRE(I_plain.size() == I_ff.size());
        for (unsigned int i = 0; i < I_plain.size(); ++i) {
            INFO(file << ", weighted " << weighted_bins << ", cw " << cw << ", q index " << i);
            CHECK_THAT(I_ff[i], Catch::Matchers::WithinRel(I_plain[i], tolerance));
        }
    }
}
