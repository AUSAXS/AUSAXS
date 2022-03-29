#include <em/ImageStack.h>
#include <plots/PlotImage.h>
#include <plots/PlotIntensity.h>
#include <fitter/SimpleIntensityFitter.h>

#include "catch2/catch.hpp"

TEST_CASE("extract_image", "[em],[files],[manual]") {
    em::ImageStack image("data/A2M_map.ccp4"); 

    plots::PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}

TEST_CASE("test_model", "[em],[files]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    em::ImageStack image("data/maptest.ccp4");
    Protein protein("data/maptest.pdb");

    image.fit(protein.get_histogram());
}

TEST_CASE("staining_and_limits", "[em],[files]") {
    SECTION("maptest.ccp4") {
        em::ImageStack image("data/maptest.ccp4");
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, setting::fit::q_high));
    }

    SECTION("A2M_map.ccp4") {
        em::ImageStack image("data/A2M_map.ccp4");
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, setting::fit::q_high));
    }

    SECTION("native10.ccp4") {
        em::ImageStack image("data/native10.ccp4", 10);
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, 2*M_PI/10));
    }

    SECTION("native25.ccp4") {
        em::ImageStack image("data/native25.ccp4", 25);
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, 2*M_PI/25));
    }
}

TEST_CASE("plot_pdb_as_points", "[em],[files]") {
    Protein protein("data/maptest.pdb");

    auto h = protein.get_histogram();
    SAXSDataset data = h.calc_debye_scattering_intensity();
    data.reduce(100);        // reduce to 100 datapoints
    data.simulate_errors();  // simulate y errors
    data.scale_errors(1000); // scale all errors so we can actually see them

    for (unsigned int i = 0; i < data.size(); i++) {
        std::cout << "VALS: " << data.x[i] << ", " << data.y[i] << std::endl;
        std::cout << "ERRS: " << data.xerr[i] << ", " << data.yerr[i] << std::endl;
    }

    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve
    plot.plot_intensity(data);                          // plot simulated data points
    plot.save("plot_pdb_as_points_test.pdf");
}