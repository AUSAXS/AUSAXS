#include <em/ImageStack.h>
#include <plots/PlotImage.h>
#include <plots/PlotIntensity.h>
#include <plots/PlotIntensityFit.h>
#include <plots/PlotIntensityFitResiduals.h>
#include <fitter/SimpleIntensityFitter.h>

#include "catch2/catch.hpp"

TEST_CASE("extract_image", "[em],[files],[manual]") {
    em::ImageStack image("data/A2M_map.ccp4"); 

    plots::PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}

TEST_CASE("test_model", "[em],[files],[slow]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::max_atoms = 50000;
    em::ImageStack image("data/native10.ccp4");
    Protein protein("data/native.pdb");
    auto res = image.fit(protein.get_histogram());

    // set optimal cutoff
    std::cout << "Optimal cutoff is " << res->params.at("cutoff") << std::endl;

    // Fit intensity plot (debug, should be equal to the next one)
    plots::PlotIntensity plot_i(protein.get_histogram(), kBlack);
    plot_i.plot_intensity(res, kBlue);
    plot_i.save("em_intensity.pdf");

    // Fit plot
    plots::PlotIntensityFit plot_f(res);
    plot_f.save("em_intensity_fit." + setting::figures::format);

    // Residual plot
    plots::PlotIntensityFitResiduals plot_r(res);
    plot_r.save("em_residuals." + setting::figures::format);
}

TEST_CASE("check_bound_savings", "[em],[files],[slow]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    em::ImageStack image("data/native10.ccp4");

    ObjectBounds3D bounds = image.minimum_volume(1);
    std::cout << "Cutoff = 1: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(2);
    std::cout << "Cutoff = 2: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(3);
    std::cout << "Cutoff = 3: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(4);
    std::cout << "Cutoff = 4: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;
}

TEST_CASE("plot_pdb_as_points", "[em],[files]") {
    Protein protein("data/maptest.pdb");

    auto h = protein.get_histogram();
    SAXSDataset data = h.calc_debye_scattering_intensity();
    data.set_resolution(25); // set the resolution
    data.reduce(100, true);  // reduce to 100 datapoints
    data.simulate_errors();  // simulate y errors
    // data.scale_errors(1000); // scale all errors so we can actually see them

    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve
    plot.plot_intensity(data);                          // plot simulated data points
    plot.save("plot_pdb_as_points_test.pdf");
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

TEST_CASE("minimum_area", "[em]") {
    Matrix data = Matrix<float>{{0, -5, -5, -5, 0, 0}, {0, -5, 5, 5, 0, 0}, {0, 0, 5, 5, 0, 0}, {0, 5, 0, 0, 5, 0}, {0, 5, 5, 5, 0, 0}, {0, -5, 0, 5, -5, 0}};
    em::Image image(data);

    SECTION("correct_bounds") {
        ObjectBounds2D bounds = image.setup_bounds(1);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 0);
        CHECK(bounds[0].max == 5);
        CHECK(bounds[1].min == 2);
        CHECK(bounds[1].max == 3);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 3);
        CHECK(bounds[3].min == 1);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 3);
        CHECK(bounds[5].max == 3);

        bounds = image.setup_bounds(-1);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 3);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 1);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 5);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 5);
        CHECK(bounds[4].min == 0);
        CHECK(bounds[4].max == 5);
        CHECK(bounds[5].min == 1);
        CHECK(bounds[5].max == 4);
    }

    SECTION("complicated_bounds") {
        Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::Image image2(data2);

        // cutoff = 1
        ObjectBounds2D bounds = image2.setup_bounds(1);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 5);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 1);
        CHECK(bounds[2].max == 4);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 4);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 5);

        // cutoff = 2
        bounds = image2.setup_bounds(2);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 2);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 2);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 2);
        CHECK(bounds[4].max == 2);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 3);

        // cutoff = 3
        bounds = image2.setup_bounds(2);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 3);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 4);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 0);
        CHECK(bounds[4].max == 0);
        CHECK(bounds[5].min == 0);
        CHECK(bounds[5].max == 2);
    }

    SECTION("correct_bounded_area") {
        ObjectBounds2D bounds = image.setup_bounds(1);
        CHECK(bounds.total_area() == 6*6);
        CHECK(bounds.bounded_area() == 18);

        bounds = image.setup_bounds(-1);
        CHECK(bounds.bounded_area() == 26);
    }

    SECTION("correct_bounds_imagestack") {
        Matrix data2 = Matrix<float>{{0, -5, -5, -5, 0, 0}, {0, -5, 5, 5, 0, 0}, {0, 0, 5, 5, 0, 0}, {0, 5, 0, 0, 5, 0}, {0, 5, 5, 5, 0, 0}, {0, -5, 0, 5, -5, 0}};
        em::ImageStack images({data, data2});
        
        ObjectBounds3D bounds = images.minimum_volume(1);
        CHECK(bounds.total_volume() == 2*6*6);
        CHECK(bounds.bounded_volume() == 36);

        bounds = images.minimum_volume(-1);
        CHECK(bounds.bounded_volume() == 52);
    }
}