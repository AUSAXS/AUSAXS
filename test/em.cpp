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
    setting::em::max_atoms = 10000;
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

TEST_CASE("check_simulated_errors", "[em],[files]") {
    SAXSDataset data1("data/2epe.RSR");
    data1.scale_y(10000);
    SAXSDataset data2 = data1;
    data2.simulate_errors();

    data1.plot_options.color = kBlack;
    data1.plot_options.draw_line = false;
    data1.plot_options.draw_markers = true;
    data1.plot_options.linewidth = 2;

    data2.plot_options.color = kOrange+2;
    data2.plot_options.draw_line = false;
    data2.plot_options.draw_markers = true;
    data2.plot_options.linewidth = 2;

    plots::PlotIntensity plot(data1);
    plot.plot_intensity(data2);
    plot.save("temp/compare_errors.pdf");
}

TEST_CASE("dataset_can_read_rsr", "[em],[dataset],[files]") {
    SAXSDataset data("data/2epe.RSR");
    vector<double>& x = data.x;
    vector<double>& y = data.y;
    vector<double>& yerr = data.yerr;

    vector<double> validate_x = {9.81300045E-03, 1.06309997E-02, 1.14489999E-02, 1.22659998E-02, 1.30840000E-02, 1.39020002E-02, 1.47200003E-02, 1.55379996E-02, 1.63550004E-02, 1.71729997E-02};
    vector<double> validate_y = {6.67934353E-03, 7.27293547E-03, 8.74083303E-03, 9.22449585E-03, 9.13867634E-03, 9.21153929E-03, 9.37998667E-03, 8.67372658E-03, 9.23649967E-03, 9.22480784E-03};
    vector<double> validate_yerr = {1.33646582E-03, 1.01892441E-03, 8.62116576E-04, 7.71059655E-04, 6.87870081E-04, 6.30189374E-04, 4.98525158E-04, 4.69041377E-04, 4.46073769E-04, 4.26004088E-04};

    REQUIRE(x.size() == 104);
    REQUIRE(y.size() == 104);
    REQUIRE(yerr.size() == 104);
    for (unsigned int i = 0; i < validate_x.size(); i++) {
        CHECK(x[i] == Approx(validate_x[i]));
        CHECK(y[i] == Approx(validate_y[i]));
        CHECK(yerr[i] == Approx(validate_yerr[i]));
    }
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

#include <em/PartialHistogramManager.h>
TEST_CASE("partial_histogram_manager", "[em]") {
    setting::protein::use_effective_charge = false;

    // SECTION("basic functionality works") {
    //     std::shared_ptr<em::ccp4::Header> header = std::make_shared<em::ccp4::Header>();
    //     header->cella_x = 1, header->cella_y = 1, header->cella_z = 1, header->nz = 1;

    //     Matrix data = Matrix<float>{{1, 2, 3, 4, 5, 6}, {0.5, 1.5, 2.5, 3.5, 4.5, 5.5}};
    //     em::Image image(data, header, 0);
    //     em::ImageStack images({image});

    //     em::PartialHistogramManager manager(images);
    //     manager.set_cutoff_levels({2, 4, 6, 8});
    //     manager.update_protein(0);
    //     std::shared_ptr<Protein> protein = manager.get_protein();

    //     REQUIRE(protein->body_size() == 4);
    //     CHECK(protein->bodies[0].protein_atoms.size() == 3);
    //     CHECK(protein->bodies[1].protein_atoms.size() == 4);
    //     CHECK(protein->bodies[2].protein_atoms.size() == 4);
    //     CHECK(protein->bodies[3].protein_atoms.size() == 1);

    //     manager.update_protein(3);
    //     REQUIRE(protein->body_size() == 4);
    //     CHECK(protein->bodies[0].protein_atoms.size() == 0);
    //     CHECK(protein->bodies[1].protein_atoms.size() == 2);
    //     CHECK(protein->bodies[2].protein_atoms.size() == 4);
    //     CHECK(protein->bodies[3].protein_atoms.size() == 1);
    // }

    SECTION("simple comparison with standard approach") {
        std::shared_ptr<em::ccp4::Header> header = std::make_shared<em::ccp4::Header>();
        header->cella_x = 1, header->cella_y = 1, header->cella_z = 1, header->nz = 1;

        Matrix data = Matrix<float>{{1, 2, 3, 4, 5, 6}, {0.5, 1.5, 2.5, 3.5, 4.5, 5.5}};
        em::Image image(data, header, 0);
        em::ImageStack images({image});

        em::PartialHistogramManager manager(images);
        manager.set_cutoff_levels({2, 4, 6, 8});

        // try an arbitrary cutoff level
        std::cout << "CHECKPOINT" << std::endl;
        ScatteringHistogram h1 = manager.get_histogram_slow(3);
        std::cout << "CHECKPOINT" << std::endl;
        ScatteringHistogram h2 = manager.get_histogram(3);
        REQUIRE(h1.p.size() == h2.p.size());
        for (unsigned int i = 0; i < h1.p.size(); i++) {
            if (h1.p[i] != h2.p[i]) {
                cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
                REQUIRE(false);
            }
        }

        // try a lower cutoff level
        h1 = manager.get_histogram_slow(1);
        h2 = manager.get_histogram(1);
        REQUIRE(h1.p.size() == h2.p.size());
        for (unsigned int i = 0; i < h1.p.size(); i++) {
            if (h1.p[i] != h2.p[i]) {
                cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
                REQUIRE(false);
            }
        }

        // try a higher cutoff level
        h1 = manager.get_histogram_slow(4);
        h2 = manager.get_histogram(4);
        REQUIRE(h1.p.size() == h2.p.size());
        for (unsigned int i = 0; i < h1.p.size(); i++) {
            if (h1.p[i] != h2.p[i]) {
                cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
                REQUIRE(false);
            }
        }
    }

    // SECTION("comparison with standard approach") {
    //     setting::em::max_atoms = 10000;
    //     em::ImageStack images("data/A2M_map.ccp4");
    //     em::PartialHistogramManager manager(images);

    //     // try an arbitrary cutoff level
    //     ScatteringHistogram h1 = manager.get_histogram_slow(4);
    //     ScatteringHistogram h2 = manager.get_histogram(4);
    //     REQUIRE(h1.p.size() == h2.p.size());
    //     for (unsigned int i = 0; i < h1.p.size(); i++) {
    //         if (h1.p[i] != h2.p[i]) {
    //             cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
    //             REQUIRE(false);
    //         }
    //     }

    //     // try a lower cutoff level
    //     h1 = manager.get_histogram_slow(2);
    //     h2 = manager.get_histogram(2);
    //     REQUIRE(h1.p.size() == h2.p.size());
    //     for (unsigned int i = 0; i < h1.p.size(); i++) {
    //         if (h1.p[i] != h2.p[i]) {
    //             cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
    //             REQUIRE(false);
    //         }
    //     }

    //     // try a higher cutoff level
    //     h1 = manager.get_histogram_slow(5);
    //     h2 = manager.get_histogram(5);
    //     REQUIRE(h1.p.size() == h2.p.size());
    //     for (unsigned int i = 0; i < h1.p.size(); i++) {
    //         if (h1.p[i] != h2.p[i]) {
    //             cout << "Failed on index " << i << ". Values: " << h1.p[i] << ", " << h2.p[i] << endl;
    //             REQUIRE(false);
    //         }
    //     }
    // }
}

TEST_CASE("minimum_area", "[em]") {
    Matrix data = Matrix<float>{{0, -5, -5, -5, 0, 0}, {0, -5, 5, 5, 0, 0}, {0, 0, 5, 5, 0, 0}, {0, 5, 0, 0, 5, 0}, {0, 5, 5, 5, 0, 0}, {0, -5, 0, 5, -5, 0}};
    em::Image image(data);

    SECTION("correct_bounds") {
        ObjectBounds2D bounds = image.setup_bounds(1);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 0);
        CHECK(bounds[0].max == 0);
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
        CHECK(bounds[2].max == 0);
        CHECK(bounds[3].min == 0);
        CHECK(bounds[3].max == 0);
        CHECK(bounds[4].min == 0);
        CHECK(bounds[4].max == 0);
        CHECK(bounds[5].min == 1);
        CHECK(bounds[5].max == 4);
    }

    SECTION("complicated_bounds") {
        Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::Image image2(data2);

        // cutoff = 1
        ObjectBounds2D bounds = image2.setup_bounds(1);
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

        // cutoff = 2
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
        CHECK(bounds.bounded_area() == 13);

        bounds = image.setup_bounds(-1);
        CHECK(bounds.bounded_area() == 11);
    }

    SECTION("correct_bounds_imagestack") {
        Matrix data2 = Matrix<float>{{0, -5, -5, -5, 0, 0}, {0, -5, 5, 5, 0, 0}, {0, 0, 5, 5, 0, 0}, {0, 5, 0, 0, 5, 0}, {0, 5, 5, 5, 0, 0}, {0, -5, 0, 5, -5, 0}};
        em::ImageStack images({data, data2});
        
        ObjectBounds3D bounds = images.minimum_volume(1);
        CHECK(bounds.total_volume() == 2*6*6);
        CHECK(bounds.bounded_volume() == 26);

        bounds = images.minimum_volume(-1);
        CHECK(bounds.bounded_volume() == 22);
    }
}