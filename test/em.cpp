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

TEST_CASE("extract_image", "[em],[files],[manual]") {
    em::ImageStack image("data/A2M_ma.ccp4"); 

    plots::PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}

TEST_CASE("test_model", "[em],[files],[slow],[manual]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("data/native.pdb");
    auto res = image.fit(protein.get_histogram());

    // set optimal cutoff
    std::cout << "Optimal cutoff is " << res->get_parameter("cutoff").value << std::endl;

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

    FitReporter::report(res);
}

TEST_CASE("check_chi2", "[em],[files]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    Protein protein("data/2epe.pdb");
    SAXSDataset data = protein.simulate_dataset();

    SimpleIntensityFitter fitter(data, protein.get_histogram());
    auto res = fitter.fit();
    REQUIRE_THAT(res->chi2/res->dof, Catch::Matchers::WithinAbs(1., 0.5));
    plots::PlotIntensityFit plot1(res);
    plot1.save("figures/test/em/check_chi2_1.pdf");

    // check that reduced chi2 is ~1
    std::cout << "Reduced chi2 is " << res->chi2/res->dof << std::endl;

    data.scale_errors(2);
    fitter = SimpleIntensityFitter(data, protein.get_histogram());
    res = fitter.fit();
    REQUIRE_THAT(res->chi2/res->dof, Catch::Matchers::WithinAbs(1., 0.5));
    plots::PlotIntensityFit plot2(res);
    plot2.save("figures/test/em/check_chi2_2.pdf");

    data.scale_errors(1./4);
    fitter = SimpleIntensityFitter(data, protein.get_histogram());
    res = fitter.fit();
    REQUIRE_THAT(res->chi2/res->dof, Catch::Matchers::WithinAbs(1., 0.5));
    plots::PlotIntensityFit plot3(res);
    plot3.save("figures/test/em/check_chi2_3.pdf");
}

TEST_CASE("generate_contour", "[em],[files],[slow],[manual]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("data/native.pdb");
    hist::ScatteringHistogram hist(protein.get_histogram());

    auto data = image.cutoff_scan_fit({1000, 1, 2}, hist);

    Dataset& scan = data.contour;
    Dataset& fit = data.fit.evaluated_points;
    fit.plot_options.set("markers", {{"color", kOrange+2}});

    plots::PlotDataset plot(scan);
    plot.plot(fit);
    plot.save("figures/test/em/chi2_landscape.pdf");
}

TEST_CASE("check_bound_savings", "[em],[files],[slow]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    em::ImageStack image("sim/native_10.ccp4");

    ObjectBounds3D bounds = image.minimum_volume(1);
    std::cout << "Cutoff = 1: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(2);
    std::cout << "Cutoff = 2: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(3);
    std::cout << "Cutoff = 3: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(4);
    std::cout << "Cutoff = 4: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;
}

TEST_CASE("repeat_chi2_contour", "[em],[files]") {
    unsigned int repeats = 50;

    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::fit::q_high = 0.4;

    // prepare measured data
    Protein protein("data/native.pdb");
    SAXSDataset data = protein.get_histogram().calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit(Limit(setting::fit::q_low, setting::fit::q_high));
    data.simulate_errors();

    // prepare fit data
    em::ImageStack image("sim/native_10.ccp4");
    auto hist = protein.get_histogram();

    vector<std::pair<double, double>> optvals;
    Multiset contours;
    Multiset evaluations;

    // comparison function. check if two datasets are exactly equal
    auto compare_contours = [&contours] (const Dataset& data) {
        Dataset& base = contours[0];
        REQUIRE(base.size() == data.size());
        for (unsigned int i = 0; i < base.size(); i++) {
            if (base.x[i] != data.x[i]) {
                utility::print_warning("Error: x values are not equal.");
                REQUIRE(base.x[i] == data.x[i]);
            }

            if (base.y[i] != data.y[i]) {
                utility::print_warning("Error: y values are not equal.");
                REQUIRE(base.y[i] == data.y[i]);
            }
        }
    };

    // SECTION("check_equality_no_noise") {
    //     setting::em::sample_frequency = 2;
    //     setting::em::simulation::noise = false;
    //     for (unsigned int i = 0; i < repeats; i++) {
    //         Dataset contour = image.cutoff_scan({10, 0, 6}, hist);
    //         contours.push_back(contour);
    //         compare_contours(contour);
    //     }
    // }

    SECTION("with_noise") {
        setting::em::simulation::noise = true;
        for (unsigned int i = 0; i < repeats; i++) {
            auto data = image.cutoff_scan_fit({100, 1.5, 4.5}, hist);
            Dataset& fit = data.fit.evaluated_points;
            optvals.push_back({fit.x[fit.size()-1], fit.y[fit.size()-1]});

            // chi2 contour plot
            Dataset& scan = data.contour;
            fit.plot_options.set("markers", {{"color", kOrange+2}});

            plots::PlotDataset plot_c(scan);
            plot_c.plot(fit);
            plot_c.save("figures/test/em/repeat_chi2_contours/" + std::to_string(i) + ".png");

            scan.plot_options.set("line", {{"color", kBlack}, {"alpha", 1.}});
            contours.push_back(scan);
            evaluations.push_back(fit);
        }

        Dataset fit_mins;
        fit_mins.set_plot_options(plots::PlotOptions("markers", {{"color", kOrange+2}, {"ms", 8}, {"s", 0.8}}));
        for (const auto& val : optvals) {
            fit_mins.push_back({val.first, val.second});
            std::cout << "(x, y): " << "(" << val.first << ", " << val.second << ")" << std::endl;
        }

        Dataset scan_mins;
        scan_mins.set_plot_options(plots::PlotOptions("markers", {{"color", kBlue+2}, {"ms", 8}, {"s", 0.8}}));
        for (const Dataset& contour : contours) {
            scan_mins.push_back(contour.find_minimum());
        }

        plots::PlotDataset plot_c(contours);
        plot_c.plot(fit_mins);
        plot_c.plot(scan_mins);
        plot_c.save("figures/test/em/repeat_chi2_contours.pdf");


        // Difference plot
        REQUIRE(scan_mins.size() == fit_mins.size());

        Dataset diff;
        diff.set_plot_options(plots::PlotOptions("markers", {{"color", kOrange+2}, {"ms", 8}, {"s", 0.8}, {"xlabel", "\\Delta cutoff"}, {"ylabel", "\\Delta \\chi^{2}"}}));
        for (unsigned int i = 0; i < scan_mins.size(); i++) {
            double delta_x = fit_mins.x[i] - scan_mins.x[i];
            double delta_y = fit_mins.y[i] - scan_mins.y[i];
            diff.push_back({delta_x, delta_y});
        }
        plots::PlotDataset::quick_plot(diff, "figures/test/em/diff.pdf");        
    }
}

TEST_CASE("plot_images", "[em],[files],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::fit::q_high = 0.4;

    string file = "data/A2M_ma.map";

    em::ImageStack image("sim/native_23.ccp4");
    for (unsigned int i = 0; i < image.size(); i++) {
        plots::PlotImage plot(image.image(i));
        // plot.plot_atoms(-1);
        plot.save("figures/test/em/images/" + utility::stem(file) + "/" + std::to_string(i) + ".png");
    }
}

TEST_CASE("voxelplot", "[em],[files],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::fit::q_high = 0.4;
    em::ImageStack image("sim/native_23.ccp4");

    SECTION("check voxel count") {
        CHECK(image.image(25).count_voxels(2) == 72);
        CHECK(image.image(30).count_voxels(10) == 38);
        CHECK(image.image(35).count_voxels(10) == 74);

        // The following generates the actual plots where the number of voxels can be counted.
        // plots::PlotImage p25(image.image(25));
        // p25.plot_atoms(2);
        // p25.save("figures/test/em/voxels/25.png");
        // std::cout << "25.png: " << image.image(25).count_voxels(2) << std::endl;

        // plots::PlotImage p30(image.image(30));
        // p30.plot_atoms(10);
        // p30.save("figures/test/em/voxels/30.png");
        // std::cout << "30.png: " << image.image(30).count_voxels(10) << std::endl;

        // plots::PlotImage p35(image.image(35));
        // p35.plot_atoms(10);
        // p35.save("figures/test/em/voxels/35.png");
        // std::cout << "35.png: " << image.image(35).count_voxels(10) << std::endl;
    }

    SECTION("make voxelcount plot") {
        Dataset data;
        Axis range(1000, 0, 6);
        for (const double& val : range.as_vector()) {
            data.push_back({val, double(image.count_voxels(val))});
        }

        data.add_plot_options("markers", {{"xlabel", "cutoff"}, {"ylabel", "number of voxels"}});
        plots::PlotDataset::quick_plot(data, "figures/test/em/voxel_count.pdf"); 
    }
}

TEST_CASE("plot_pdb_as_points", "[em],[files],[manual]") {
    Protein protein("data/maptest.pdb");

    auto h = protein.get_histogram();
    SAXSDataset data = h.calc_debye_scattering_intensity();
    data.set_resolution(25); // set the resolution
    data.reduce(100, true);  // reduce to 100 datapoints
    data.simulate_errors();  // simulate y errors
    // data.scale_errors(1000); // scale all errors so we can actually see them

    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve
    plot.plot_intensity(data);                          // plot simulated data points
    plot.save("figures/test/em/plot_pdb_as_points.pdf");
}

TEST_CASE("check_simulated_errors", "[em],[files],[manual]") {
    setting::fit::q_high = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    em::ImageStack image("data/native10.ccp4");
    auto hist = image.get_histogram(2);
    auto data = hist.calc_debye_scattering_intensity();
    data.normalize(1.1);
    data.simulate_errors();
    data.save("temp/em/simulated_errors.txt");

    data.add_plot_options("errors", {{"logx", true}, {"logy", true}});
    plots::PlotDataset plot(data);
    plot.save("temp/em/check_errors.pdf");
}

TEST_CASE("staining_and_limits", "[em],[files]") {
    SECTION("maptest.ccp4") {
        em::ImageStack image("data/maptest.ccp4");
        CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, setting::fit::q_high));
    }

    SECTION("A2M_ma.ccp4") {
        em::ImageStack image("data/A2M_ma.ccp4");
        CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, setting::fit::q_high));
    }

    SECTION("native10.ccp4") {
        em::ImageStack image("data/native10.ccp4", 10);
        CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::fit::q_low, 2*M_PI/10));
    }

    SECTION("native25.ccp4") {
        em::ImageStack image("data/native25.ccp4", 25);
        CHECK(image.positively_stained());
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