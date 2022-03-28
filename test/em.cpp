#include <em/ImageStack.h>
#include <plots/PlotImage.h>
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
        CHECK(!image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, setting::fit::q_high));
    }

    SECTION("native10.ccp4") {
        em::ImageStack image("data/native10.ccp4");
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, 2*M_PI/10));
    }

    SECTION("native25.ccp4") {
        em::ImageStack image("data/native25.ccp4");
        CHECK(image.is_positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, 2*M_PI/25));
    }
}