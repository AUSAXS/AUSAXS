#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>

#include <em/Image.h>
#include <em/ImageStack.h>
#include <em/ObjectBounds3D.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>
#include <dataset/Multiset.h>
#include <utility/Utility.h>
#include <utility/Console.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Water.h>
#include <data/Atom.h>
#include <em/detail/ExtendedLandscape.h>
#include <io/ExistingFile.h>
#include <settings/All.h>
#include <utility/Constants.h>

using std::vector;

TEST_CASE("extract_image", "[manual]") {
    em::ImageStack image("test/files/A2M_2020_Q4.ccp4"); 

    plots::PlotImage plot(image.image(5));
    // plot.plot_atoms(0.1);
    plot.save("test.pdf");
}

TEST_CASE("test_model", "[slow],[manual]") {
    settings::axes::qmax = 0.4;
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("test/A2M_native/native.pdb");
    auto res = image.fit(protein.get_histogram());

    // set optimal cutoff
    std::cout << "Optimal cutoff is " << res->get_parameter("cutoff").value << std::endl;

    // Fit intensity plot (debug, should be equal to the next one)
    plots::PlotIntensity plot_i(protein.get_histogram(), style::color::black);
    plot_i.plot(res, style::color::blue);
    plot_i.save("em_intensity.pdf");

    // Fit plot
    plots::PlotIntensityFit plot_f(res);
    plot_f.save("em_intensity_fit.pdf");

    // Residual plot
    plots::PlotIntensityFitResiduals plot_r(res);
    plot_r.save("em_residuals.pdf");

    fitter::FitReporter::report(res);
}

TEST_CASE("generate_contour", "[files],[slow],[manual]") {
    settings::axes::qmax = 0.4;
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 2;
    em::ImageStack image("sim/native_10.ccp4");
    Protein protein("data/A2M_native/native.pdb");
    hist::ScatteringHistogram hist = protein.get_histogram();

    auto[r, s] = image.cutoff_scan_fit({1, 2, 1000}, hist);
    auto fit = r.evaluated_points.as_dataset();
    auto scan = s.as_dataset();
    fit.add_plot_options("markers", {{"color", style::color::orange}});

    plots::PlotDataset plot(scan);
    plot.plot(fit);
    plot.save("figures/test/em/chi2_landscape.pdf");
}

/**
 * @brief Generate a contour plot of the chi2 landscape for the EM fit. 
 */
TEST_CASE("check_fit", "[files],[manual],[slow]") {
    std::string mfile = "data/SASDDD3/SASDDD3.dat";
    std::string mapfile = "data/SASDDD3/emd_0560.map";

    // read custom settings
    std::string path = std::filesystem::path(mfile).parent_path().string();
    if (std::filesystem::exists(path + "/settings.txt")) {
        std::cout << "Using discovered settings file at " << path << "/settings.txt" << std::endl;
        settings::read(path + "/settings.txt");
    }

    settings::em::sample_frequency = 2;
    settings::protein::use_effective_charge = false;
    settings::fit::verbose = true;
    settings::em::hydrate = true;
    em::ImageStack map(mapfile);

    auto[r, s] = map.cutoff_scan_fit({0.025, 0.03, 1000}, mfile);
    auto fit = r.evaluated_points.as_dataset();
    auto scan = s.as_dataset();
    fit.add_plot_options(style::draw::points, {{"color", style::color::orange}});

    auto fitted_water_factors = map.get_fitted_water_factors_dataset();
    plots::PlotDataset::quick_plot(fitted_water_factors, "figures/test/em/check_fit_landscape_wf.pdf");

    plots::PlotDataset plot(scan);
    plot.plot(fit);
    plot.save("figures/test/em/check_fit_landscape.pdf");
}

TEST_CASE("check_bound_savings", "[manual],[slow]") {
    settings::axes::qmax = 0.4;
    settings::protein::use_effective_charge = false;
    em::ImageStack image("sim/native_10.ccp4");

    em::ObjectBounds3D bounds = image.minimum_volume(1);
    std::cout << "Cutoff = 1: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(2);
    std::cout << "Cutoff = 2: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(3);
    std::cout << "Cutoff = 3: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;

    bounds = image.minimum_volume(4);
    std::cout << "Cutoff = 4: Using " << bounds.bounded_volume() << " of " << bounds.total_volume() << " voxels." << std::endl;
}

TEST_CASE("repeat_chi2_contour", "[files],[slow],[manual]") {
    unsigned int repeats = 50;

    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;
    settings::axes::qmax = 0.4;

    // prepare measured data
    Protein protein("data/A2M_native/native.pdb");
    SimpleDataset data = protein.get_histogram().calc_debye_scattering_intensity();
    data.reduce(settings::fit::N, true);
    data.limit_x(Limit(settings::axes::qmin, settings::axes::qmax));
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
                console::print_warning("Error: x values are not equal.");
                REQUIRE(base.x(i) == data.x(i));
            }

            if (base.y(i) != data.y(i)) {
                console::print_warning("Error: y values are not equal.");
                REQUIRE(base.y(i) == data.y(i));
            }
        }
    };

    SECTION("check_equality_no_noise") {
        settings::em::sample_frequency = 2;
        settings::em::simulation::noise = false;
        for (unsigned int i = 0; i < repeats; i++) {
            auto landscape = image.cutoff_scan({10, 0, 6}, hist).as_dataset();
            contours.push_back(landscape);
            compare_contours(landscape);
        }
    }

    SECTION("with_noise") {
        settings::em::simulation::noise = true;
        for (unsigned int i = 0; i < repeats; i++) {
            auto[r, s] = image.cutoff_scan_fit({1.5, 4.5, 100}, hist);
            auto fit = r.evaluated_points.as_dataset();
            optvals.push_back({fit.x(fit.size()-1), fit.y(fit.size()-1)});

            // chi2 contour plot
            auto scan = s.as_dataset();
            fit.add_plot_options(style::draw::points, {{"color", style::color::orange}});

            plots::PlotDataset plot_c(scan);
            plot_c.plot(fit);
            plot_c.save("figures/test/em/repeat_chi2_contours/" + std::to_string(i) + ".png");

            scan.add_plot_options(style::draw::line, {{"color", style::color::black}});
            contours.push_back(scan);
            evaluations.push_back(fit);
        }

        Dataset2D fit_mins;
        fit_mins.set_plot_options(plots::PlotOptions(style::draw::points, {{"color", style::color::orange}, {"ms", 8}, {"s", 0.8}}));
        for (const auto& val : optvals) {
            fit_mins.push_back(val.first, val.second);
            std::cout << "(x, y): " << "(" << val.first << ", " << val.second << ")" << std::endl;
        }

        Dataset2D scan_mins;
        scan_mins.set_plot_options(plots::PlotOptions(style::draw::points, {{"color", style::color::blue}, {"ms", 8}, {"s", 0.8}}));
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
        diff.set_plot_options(plots::PlotOptions(style::draw::points, {{"color", style::color::orange}, {"ms", 8}, {"s", 0.8}, {"xlabel", "\\Delta cutoff"}, {"ylabel", "\\Delta \\chi^{2}"}}));
        for (unsigned int i = 0; i < scan_mins.size(); i++) {
            double delta_x = fit_mins.x(i) - scan_mins.x(i);
            double delta_y = fit_mins.y(i) - scan_mins.y(i);
            diff.push_back({delta_x, delta_y});
        }
        plots::PlotDataset::quick_plot(diff, "figures/test/em/diff.pdf");        
    }
}

TEST_CASE("plot_images", "[files],[manual],[slow]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;
    settings::axes::qmax = 0.4;

    io::ExistingFile file = "test/files/A2M_2020_Q4.ccp4";
    settings::plots::contour = {-100, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 19, 22, 25};
    em::ImageStack image(file);
    for (unsigned int i = 0; i < image.size(); i++) {
        plots::PlotImage plot(image.image(i));
        // plot.plot_atoms(-1);
        plot.save("figures/test/em/images/" + file.stem() + "/" + std::to_string(i) + ".png");
    }
}

TEST_CASE("get_histogram", "[manual]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;

    std::string file = "test/files/A2M_2020_Q4.ccp4";
    std::cout << "DEBUG" << std::endl;
    em::ImageStack image(file);
    std::cout << "DEBUG" << std::endl;
    auto hist = image.get_histogram(2);
    std::cout << "DEBUG" << std::endl;
    plots::PlotHistogram::quick_plot(hist, "figures/test/em/histogram.pdf");
}

TEST_CASE("voxelplot", "[manual]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;
    settings::axes::qmax = 0.4;
    em::ImageStack image("test/files/A2M_2020_Q4.ccp4");

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

TEST_CASE("voxelcount", "[manual]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;
    settings::axes::qmax = 0.4;
    em::ImageStack image("data/emd_24889/emd_24889.map");

    Dataset2D data;
    Axis range(image.from_level(1), image.from_level(20), 1000);
    for (const double& val : range.as_vector()) {
        data.push_back({val, double(image.count_voxels(val))});
    }

    data.add_plot_options("markers", {{"xlabel", "cutoff"}, {"ylabel", "number of voxels"}, {"logy", true}});
    plots::PlotDataset::quick_plot(data, "temp/test/em/voxel_count.png"); 
}

TEST_CASE("mass_cutoff_plot", "[manual]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 1;
    settings::axes::qmax = 0.4;
    em::ImageStack image("data/emd_24889/emd_24889.map");

    Dataset2D data;
    Axis range(image.from_level(0.5), image.from_level(7), 100);
    for (const double& val : range.as_vector()) {
        auto protein = image.get_protein(val);
        protein->generate_new_hydration();
        data.push_back({val, protein->get_volume_grid()*constants::mass::density::protein});
        // data.push_back({val, image.get_protein(val)->get_volume_grid()*constants::mass::density::protein});
    }

    data.add_plot_options("markers", {{"xlabel", "cutoff"}, {"ylabel", "mass"}});
    plots::PlotDataset::quick_plot(data, "temp/em/mass_cutoff.png"); 
}

TEST_CASE("instability", "[files],[manual]") {
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 2;
    settings::axes::qmax = 0.4;
    em::ImageStack image("data/emd_12747/emd_12747.map");

    SimpleDataset data;
    Axis range(image.from_level(0.5), image.from_level(7), 100);
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

TEST_CASE("save_as_pdb", "[manual]") {
    em::ImageStack image("data/A2M_2020_Q4/A2M_2020_Q4.ccp4");
    image.get_protein(image.from_level(3))->save("figures/test/em/save_as_pdb.pdb");
}

TEST_CASE("plot_pdb_as_points", "[files],[manual]") {
    Protein protein("data/maptest.pdb");

    auto h = protein.get_histogram();
    SimpleDataset data = h.calc_debye_scattering_intensity();
    // data.set_resolution(25); // set the resolution //! ???
    data.reduce(100, true);  // reduce to 100 datapoints
    data.simulate_errors();  // simulate y errors
    // data.scale_errors(1000); // scale all errors so we can actually see them

    plots::PlotIntensity plot(protein.get_histogram()); // plot actual curve
    plot.plot(data);                          // plot simulated data points
    plot.save("figures/test/em/plot_pdb_as_points.pdf");
}

TEST_CASE("check_simulated_errors", "[files],[manual],[broken]") {
    settings::axes::qmax = 0.4;
    settings::protein::use_effective_charge = false;
    settings::em::sample_frequency = 2;

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