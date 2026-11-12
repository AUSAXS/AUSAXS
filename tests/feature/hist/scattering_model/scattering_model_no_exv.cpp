#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <data/Molecule.h>
#include <hist/scattering_model/IScatteringModelExv.h>
#include <settings/All.h>

#include <support/hist_test_helper.h>

#include <string>

using namespace ausaxs;

// without a solvent, the excluded volume of any model carries no density, so it must reduce to the vacuum model
TEST_CASE("ScatteringModelNoExv: agrees with an exv model at zero solvent density") {
    settings::general::verbose = false;
    std::string file = GENERATE("2epe", "6lyz", "LAR1-2", "SASDJQ4");
    bool weighted_bins = GENERATE(true, false);
    settings::hist::weighted_bins = weighted_bins; // pinned, since auto resolves differently for the two models
    SECTION(file + (weighted_bins ? " (weighted bins)" : " (unweighted bins)")) {
        data::Molecule protein("tests/files/" + file + ".pdb");

        settings::exv::exv_method = settings::exv::ExvMethod::Vacuum;
        protein.reset_histogram_manager();
        auto I_vacuum = protein.get_histogram()->debye_transform();

        settings::exv::exv_method = settings::exv::ExvMethod::Fraser;
        protein.reset_histogram_manager();
        auto hist = protein.get_histogram();
        auto exv_model = dynamic_cast<hist::IScatteringModelExv*>(hist.get());
        REQUIRE(exv_model != nullptr);
        exv_model->apply_solvent_density_scaling_factor(0);
        auto I_exv = exv_model->debye_transform();

        REQUIRE(I_vacuum.size() == I_exv.size());
        REQUIRE(compare_hist(I_vacuum.get_intensity(), I_exv.get_intensity(), 0, 1e-9));
    }
}
