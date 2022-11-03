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
#include <filesystem>

using std::vector;

TEST_CASE("extract_image", "[em],[files],[manual]") {
    em::ImageStack image("data/A2M_ma/A2M_ma.ccp4"); 

    plots::PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}

TEST_CASE("test_model", "[em],[files],[slow],[manual]") {
    setting::axes::qmax = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("data/A2M_native/native.pdb");
    auto res = image.fit(protein.get_histogram());

    // set optimal cutoff
    std::cout << "Optimal cutoff is " << res->get_parameter("cutoff").value << std::endl;

    // Fit intensity plot (debug, should be equal to the next one)
    plots::PlotIntensity plot_i(protein.get_histogram(), color::black);
    plot_i.plot(res);
    plot_i.save("em_intensity.pdf");

    // Fit plot
    plots::PlotIntensityFit plot_f(res);
    plot_f.save("em_intensity_fit." + setting::figures::format);

    // Residual plot
    plots::PlotIntensityFitResiduals plot_r(res);
    plot_r.save("em_residuals." + setting::figures::format);

    FitReporter::report(res);
}

TEST_CASE("generate_contour", "[em],[files],[slow],[manual]") {
    setting::axes::qmax = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("data/A2M_native/native.pdb");
    hist::ScatteringHistogram hist = protein.get_histogram();

    auto data = image.cutoff_scan_fit({1000, 1, 2}, hist);

    SimpleDataset& scan = data.contour;
    SimpleDataset& fit = data.fit.evaluated_points;
    fit.add_plot_options("markers", {{"color", color::orange}});

    plots::PlotDataset plot(scan);
    plot.plot(fit);
    plot.save("figures/test/em/chi2_landscape.pdf");
}

/**
 * @brief Generate a contour plot of the chi2 landscape for the EM fit. 
 */
TEST_CASE("check_fit", "[em],[files],[manual],[slow]") {
    string mfile = "data/SASDDD3/SASDDD3.dat";
    string mapfile = "data/SASDDD3/emd_0560.map";

    // read custom settings
    std::string path = std::filesystem::path(mfile).parent_path().string();
    if (std::filesystem::exists(path + "/settings.txt")) {
        std::cout << "Using discovered settings file at " << path << "/settings.txt" << std::endl;
        setting::read(path + "/settings.txt");
    }

    setting::em::sample_frequency = 2;
    setting::protein::use_effective_charge = false;
    setting::fit::verbose = true;
    setting::em::hydrate = true;
    em::ImageStack map(mapfile);

    auto data = map.cutoff_scan_fit({1000, 0.025, 0.03}, mfile);
    SimpleDataset& scan = data.contour;
    SimpleDataset& fit = data.fit.evaluated_points;
    fit.add_plot_options("markers", {{"color", color::orange}});

    auto fitted_water_factors = map.get_fitted_water_factors_dataset();
    plots::PlotDataset::quick_plot(fitted_water_factors, "figures/test/em/check_fit_landscape_wf.pdf");

    plots::PlotDataset plot(scan);
    plot.plot(fit);
    plot.save("figures/test/em/check_fit_landscape.pdf");
}

TEST_CASE("check_bound_savings", "[em],[files],[slow]") {
    setting::axes::qmax = 0.4;
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

TEST_CASE("repeat_chi2_contour", "[em],[files],[slow],[manual]") {
    unsigned int repeats = 50;

    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::axes::qmax = 0.4;

    // prepare measured data
    Protein protein("data/A2M_native/native.pdb");
    SimpleDataset data = protein.get_histogram().calc_debye_scattering_intensity();
    data.reduce(setting::fit::N, true);
    data.limit_x(Limit(setting::axes::qmin, setting::axes::qmax));
    data.simulate_errors();

    // prepare fit data
    em::ImageStack image("sim/native_10.ccp4");
    auto hist = protein.get_histogram();

    vector<std::pair<double, double>> optvals;
    Multiset contours;
    Multiset evaluations;

    // comparison function. check if two datasets are exactly equal
    auto compare_contours = [&contours] (const Dataset2D& data) {
        Dataset2D& base = contours[0];
        REQUIRE(base.size() == data.size());
        for (unsigned int i = 0; i < base.size(); i++) {
            if (base.x(i) != data.x(i)) {
                utility::print_warning("Error: x values are not equal.");
                REQUIRE(base.x(i) == data.x(i));
            }

            if (base.y(i) != data.y(i)) {
                utility::print_warning("Error: y values are not equal.");
                REQUIRE(base.y(i) == data.y(i));
            }
        }
    };

    SECTION("check_equality_no_noise") {
        setting::em::sample_frequency = 2;
        setting::em::simulation::noise = false;
        for (unsigned int i = 0; i < repeats; i++) {
            auto contour = image.cutoff_scan({10, 0, 6}, hist);
            contours.push_back(contour.contour);
            compare_contours(contour.contour);
        }
    }

    SECTION("with_noise") {
        setting::em::simulation::noise = true;
        for (unsigned int i = 0; i < repeats; i++) {
            auto data = image.cutoff_scan_fit({100, 1.5, 4.5}, hist);
            SimpleDataset& fit = data.fit.evaluated_points;
            optvals.push_back({fit.x(fit.size()-1), fit.y(fit.size()-1)});

            // chi2 contour plot
            Dataset2D& scan = data.contour;
            fit.add_plot_options("markers", {{"color", kOrange+2}});

            plots::PlotDataset plot_c(scan);
            plot_c.plot(fit);
            plot_c.save("figures/test/em/repeat_chi2_contours/" + std::to_string(i) + ".png");

            scan.add_plot_options("line", {{"color", kBlack}, {"alpha", 1.}});
            contours.push_back(scan);
            evaluations.push_back(fit);
        }

        Dataset2D fit_mins;
        fit_mins.set_plot_options(plots::PlotOptions("markers", {{"color", kOrange+2}, {"ms", 8}, {"s", 0.8}}));
        for (const auto& val : optvals) {
            fit_mins.push_back(val.first, val.second);
            std::cout << "(x, y): " << "(" << val.first << ", " << val.second << ")" << std::endl;
        }

        Dataset2D scan_mins;
        scan_mins.set_plot_options(plots::PlotOptions("markers", {{"color", kBlue+2}, {"ms", 8}, {"s", 0.8}}));
        for (const Dataset2D& contour : contours) {
            scan_mins.push_back(contour.find_minimum());
        }

        plots::PlotDataset plot_c(contours);
        plot_c.plot(fit_mins);
        plot_c.plot(scan_mins);
        plot_c.save("figures/test/em/repeat_chi2_contours.pdf");


        // Difference plot
        REQUIRE(scan_mins.size() == fit_mins.size());

        Dataset2D diff;
        diff.set_plot_options(plots::PlotOptions("markers", {{"color", kOrange+2}, {"ms", 8}, {"s", 0.8}, {"xlabel", "\\Delta cutoff"}, {"ylabel", "\\Delta \\chi^{2}"}}));
        for (unsigned int i = 0; i < scan_mins.size(); i++) {
            double delta_x = fit_mins.x(i) - scan_mins.x(i);
            double delta_y = fit_mins.y(i) - scan_mins.y(i);
            diff.push_back({delta_x, delta_y});
        }
        plots::PlotDataset::quick_plot(diff, "figures/test/em/diff.pdf");        
    }
}

TEST_CASE("plot_images", "[em],[files],[manual],[slow]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::axes::qmax = 0.4;

    setting::plot::image::contour = {-100, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 19, 22, 25};

    string file = "data/A2M_ma/A2M_ma.ccp4";
    em::ImageStack image(file);
    for (unsigned int i = 0; i < image.size(); i++) {
        plots::PlotImage plot(image.image(i));
        // plot.plot_atoms(-1);
        plot.save("figures/test/em/images/" + utility::stem(file) + "/" + std::to_string(i) + ".png");
    }
}

TEST_CASE("get_histogram", "[em],[files],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;

    string file = "data/A2M_ma/A2M_ma.ccp4";
    em::ImageStack image(file);
    auto hist = image.get_histogram(2);
    plots::PlotHistogram::quick_plot(hist, "figures/test/em/histogram.pdf");
}

TEST_CASE("voxelplot", "[em],[files],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::axes::qmax = 0.4;
    em::ImageStack image("data/A2M_ma/A2M_ma.ccp4");

    CHECK(image.image(95).count_voxels(15) == 32);
    CHECK(image.image(100).count_voxels(10) == 208);
    CHECK(image.image(105).count_voxels(10) == 84);

    // The following generates the actual plots where the number of voxels can be manually counted.
    // plots::PlotImage p95(image.image(95));
    // p95.plot_atoms(15);
    // p95.save("figures/test/em/voxels/95.png");
    // std::cout << "95.png: " << image.image(95).count_voxels(15) << std::endl;

    // plots::PlotImage p100(image.image(100));
    // p100.plot_atoms(10);
    // p100.save("figures/test/em/voxels/100.png");
    // plots::PlotImage::quick_plot(image.image(100), "figures/test/em/voxels/100_clean.png");
    // std::cout << "100.png: " << image.image(100).count_voxels(10) << std::endl;

    // plots::PlotImage p105(image.image(105));
    // p105.plot_atoms(10);
    // p105.save("figures/test/em/voxels/105.png");
    // std::cout << "105.png: " << image.image(105).count_voxels(10) << std::endl;
}

TEST_CASE("voxelcount", "[em],[slow],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 1;
    setting::axes::qmax = 0.4;
    em::ImageStack image("data/A2M_ma/A2M_ma.ccp4");

    Dataset2D data;
    Axis range(1000, 0, 15);
    for (const double& val : range.as_vector()) {
        data.push_back({val, double(image.count_voxels(val))});
    }

    data.add_plot_options("markers", {{"xlabel", "cutoff"}, {"ylabel", "number of voxels"}});
    plots::PlotDataset::quick_plot(data, "figures/test/em/voxel_count.pdf"); 
}

TEST_CASE("instability", "[em],[files],[manual]") {
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;
    setting::axes::qmax = 0.4;
    em::ImageStack image("data/A2M_native/emd_12747.ccp4");

    SimpleDataset data;
    Axis range(100, image.level(0.5), image.level(7));
    unsigned int prev = image.count_voxels(range.min);
    for (unsigned int i = 1; i < range.bins; i++) {
        double cutoff = range.min + i * range.step();
        unsigned int count = image.count_voxels(cutoff);
        data.push_back(cutoff, prev - count);
        prev = count;
    }

    data.add_plot_options("markers", {{"xlabel", "cutoff"}, {"ylabel", "change"}});
    plots::PlotDataset::quick_plot(data, "figures/test/em/instability.pdf"); 
}

TEST_CASE("save_as_pdb", "[em],[manual]") {
    em::ImageStack image("data/A2M_ma/A2M_ma.ccp4");
    image.get_protein(image.level(3))->save("figures/test/em/save_as_pdb.pdb");
}

TEST_CASE("plot_pdb_as_points", "[em],[files],[manual]") {
    Protein protein("data/maptest.pdb");

    auto h = protein.get_histogram();
    SimpleDataset data = h.calc_debye_scattering_intensity();
    // data.set_resolution(25); // set the resolution //! ???
    data.reduce(100, true);  // reduce to 100 datapoints
    data.simulate_errors();  // simulate y errors
    // data.scale_errors(1000); // scale all errors so we can actually see them

    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve
    plot.plot_intensity(data);                          // plot simulated data points
    plot.save("figures/test/em/plot_pdb_as_points.pdf");
}

TEST_CASE("check_simulated_errors", "[em],[files],[manual],[broken]") {
    setting::axes::qmax = 0.4;
    setting::protein::use_effective_charge = false;
    setting::em::sample_frequency = 2;

    em::ImageStack image("sim/native_10.ccp4");
    auto hist = image.get_histogram(2);
    auto data = hist.calc_debye_scattering_intensity();
    data.normalize(1.1);
    data.simulate_errors();
    data.save("temp/em/simulated_errors.txt");

    data.add_plot_options("errors", {{"logx", true}, {"logy", true}});
    plots::PlotDataset plot(data);
    plot.save("temp/em/check_errors.pdf");
}

TEST_CASE("staining_and_limits", "[em],[files],[broken]") {
    SECTION("maptest.ccp4") {
        em::ImageStack image("data/maptest.ccp4");
        // CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::axes::qmin, setting::axes::qmax));
    }

    SECTION("A2M_ma.ccp4") {
        em::ImageStack image("data/A2M_ma/A2M_ma.ccp4");
        // CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::axes::qmin, setting::axes::qmax));
    }

    SECTION("native10.ccp4") {
        em::ImageStack image("sim/native_10.ccp4", 10);
        // CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::axes::qmin, 2*M_PI/10));
    }

    SECTION("native23.ccp4") {
        em::ImageStack image("sim/native_23.ccp4", 23);
        // CHECK(image.positively_stained());
        CHECK(image.get_limits() == Limit(setting::axes::qmin, 2*M_PI/23));
    }
}

TEST_CASE("minimum_area", "[em]") {
    SECTION("correct_bounds") {
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        em::Image image(data);

        ObjectBounds2D bounds = image.setup_bounds(1);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 1);
        CHECK(bounds[0].max == 4);
        CHECK(bounds[1].min == 1);
        CHECK(bounds[1].max == 3);
        CHECK(bounds[2].min == 2);
        CHECK(bounds[2].max == 4);
        CHECK(bounds[3].min == 1);
        CHECK(bounds[3].max == 4);
        CHECK(bounds[4].min == 1);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 1);
        CHECK(bounds[5].max == 5);

        bounds = image.setup_bounds(5);
        REQUIRE(bounds.size() == 6);
        CHECK(bounds[0].min == 3);
        CHECK(bounds[0].max == 3);
        CHECK(bounds[1].min == 2);
        CHECK(bounds[1].max == 3);
        CHECK(bounds[2].min == 0);
        CHECK(bounds[2].max == 0);
        CHECK(bounds[3].min == 3);
        CHECK(bounds[3].max == 3);
        CHECK(bounds[4].min == 3);
        CHECK(bounds[4].max == 3);
        CHECK(bounds[5].min == 5);
        CHECK(bounds[5].max == 5);
    }

    SECTION("more bounds") {
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
        bounds = image2.setup_bounds(3);
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
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        em::Image image(data);

        ObjectBounds2D bounds = image.setup_bounds(1);
        CHECK(bounds.total_area() == 6*6);
        CHECK(bounds.bounded_area() == (4 + 3 + 3 + 4 + 3 + 5));

        bounds = image.setup_bounds(2);
        CHECK(bounds.bounded_area() == (2 + 3 + 2 + 3 + 2 + 3));
    }

    SECTION("correct_bounds_imagestack") {
        Matrix data = Matrix<float>{{0, 1, 3, 5, 1, 0}, {0, 3, 5, 5, 0, 0}, {0, 0, 1, 3, 3, 0}, {0, 3, 0, 5, 1, 0}, {0, 1, 3, 5, 0, 0}, {0, 1, 0, 3, 1, 5}};
        Matrix data2 = Matrix<float>{{0, 1, 2, 3, 2, 1}, {0, 3, 2, 1, 3, 0}, {0, 1, 2, 0, 1, 0}, {2, 0, 0, 3, 1, 0}, {0, 1, 2, 1, 1, 0}, {3, 3, 3, 2, 1, 1}};
        em::ImageStack images({data, data2});
        
        ObjectBounds3D bounds = images.minimum_volume(1);
        CHECK(bounds.total_volume() == 2*6*6);
        CHECK(bounds.bounded_volume() == ((4 + 3 + 3 + 4 + 3 + 5) + (5 + 4 + 4 + 5 + 4 + 6)));

        bounds = images.minimum_volume(2);
        CHECK(bounds.bounded_volume() == ((2 + 3 + 2 + 3 + 2 + 3) + (3 + 4 + 1 + 4 + 1 + 4)));
    }
}