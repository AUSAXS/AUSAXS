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
    setting::protein::use_effective_charge = false;
    em::ImageStack image("data/maptest.ccp4");
    Protein protein("data/maptest.pdb");

    image.fit(protein.get_histogram());
}