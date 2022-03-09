#include <em/ImageStack.h>
#include <plots/PlotImage.h>

#include "catch2/catch.hpp"

TEST_CASE("extract_image", "[em],[files]") {
    em::ImageStack image("data/A2M_map.ccp4"); 

    PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}