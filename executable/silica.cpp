
#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <data/Protein.h>
#include <hist/ScatteringHistogram.h>
#include <mini/detail/Parameter.h>
#include <plots/all.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <fitter/FitReporter.h>
#include <crystal/Fval.h>
#include <mini/dlibMinimizer.h>
#include <settings/All.h>

int main(int argc, char const *argv[]) {
    //###############################################//
    //###                 setup                   ###//
    //###############################################//
    io::File crystal;
    CLI::App app{"Silica project"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/silica/");
    CLI11_PARSE(app, argc, argv);
    settings::general::output += crystal.stem() + "/";
    settings::axes::qmin = 1e-2;
    settings::axes::qmax = 0.5;
    settings::axes::bins = 100;
    settings::crystal::grid_expansion = 1;
    settings::crystal::h = 10000; settings::crystal::k = 10000; settings::crystal::l = 0;
    settings::crystal::max_q = 0.5;
    settings::protein::use_effective_charge = false;
    settings::axes::max_distance = 14000;
    // settings::crystal::miller_generation_strategy = settings::crystal::MillerGenerationChoice::Reduced;
    // settings::crystal::reduced::basis_q = 5;

    //###############################################//
    //###               fourier plot              ###//
    //###############################################//
    // xy slice
    crystal::CrystalScattering cs(crystal);
    auto fourierxy = cs.calculate();
    fourierxy.limit_y(1e-4, std::numeric_limits<double>::max());
    fourierxy.add_plot_options({{plots::option::legend, "xy"}});
    plots::PlotIntensity::quick_plot(fourierxy, settings::general::output + "fourierxy.png");

    // complete box
    // settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 10;
    // settings::crystal::max_q = 100;
    // settings::axes::bins = 1000;
    // crystal::CrystalScattering cs2(crystal);
    // auto fourier = cs2.calculate();
    // fourier.limit_y(1e-4, std::numeric_limits<double>::max());
    // fourier.add_plot_options({{plots::option::legend, "xyz"}});

    // plots::PlotIntensity plot(fourierxy, style::color::orange);
    // plot.plot(fourier, style::color::blue);
    // plot.save(settings::general::output + "xy_xyz_comparison.png");
    // plots::PlotIntensity::quick_plot(fourier, settings::general::output + "fourier.png");

    //###############################################//
    //###               histogram                 ###//
    //###############################################//
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    Protein protein = crystal::Fval::as_protein();
    auto hist = protein.get_histogram();
    auto debye = hist.calc_debye_scattering_intensity();
    auto hist_dataset = hist.as_dataset();
    hist_dataset.add_plot_options({{"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
    plots::PlotDataset::quick_plot(hist_dataset, settings::general::output + "distances.png");
    plots::PlotIntensity::quick_plot(debye, settings::general::output + "debye.png");

    // comparison plot
    debye.add_plot_options({{"color", style::color::red}});
    plots::PlotIntensity(fourierxy)
        .plot(debye, style::color::red)
        .save(settings::general::output + "debye_comparison.png");

    // single-particle analysis
    if (std::filesystem::exists(settings::general::output + "hmin.txt")) {
        std::cout << "hmin.txt found, performing single-particle analysis" << std::endl;

        std::ifstream input(settings::general::output + "hmin.txt");
        if (!input.is_open()) {throw except::io_error("Could not open file " + settings::general::output + "hmin.txt");}
        std::string line;
        std::getline(input, line);
        unsigned int x0 = std::stoi(line);
        std::getline(input, line);
        unsigned int xmax = std::stoi(line);
        hist_dataset.limit_x(0, xmax);

        double y0 = hist_dataset.interpolate_y(x0);

        // model is y = a*x²
        // this is constrained at the location 'x0' to the value 'y0', such that y0 = a*x0² => a = y0/x0²
        double a = y0 / (x0 * x0);
        auto f = [a] (double x) {return a*x*x;};

        // plot model on top of histogram
        auto model_dataset = hist_dataset;
        for (unsigned int i = 0; i < model_dataset.size(); i++) {
            model_dataset.y(i) = f(model_dataset.x(i));
        }

        model_dataset.add_plot_options({{"color", style::color::red}});
        plots::PlotDataset(hist_dataset)
            .plot(model_dataset)
            .save(settings::general::output + "subtracted_curve.png");

        // subtract model from histogram
        auto subtracted_hist = hist_dataset;
        for (unsigned int i = 0; i < subtracted_hist.size(); i++) {
            subtracted_hist.y(i) = std::max(subtracted_hist.y(i) - f(subtracted_hist.x(i)), 0.);
        }
        unsigned int imax = 0;
        while (subtracted_hist.x(imax) < x0) {imax++;}
        subtracted_hist.limit_x({0, subtracted_hist.x(imax)});

        Axis new_axis(0, subtracted_hist.x(imax), subtracted_hist.size());
        std::vector<double> empty(subtracted_hist.size(), 0);
        hist::ScatteringHistogram new_hist(subtracted_hist.y(), empty, empty, subtracted_hist.y(), new_axis);
        auto new_hist_dataset = new_hist.as_dataset();
        new_hist_dataset.add_plot_options({{"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
        plots::PlotDataset::quick_plot(new_hist_dataset, settings::general::output + "distances_single.png");

        auto new_intensity = new_hist.calc_debye_scattering_intensity();
        plots::PlotIntensity::quick_plot(new_intensity, settings::general::output + "intensity_single.png");
    } else {
        std::cout << "No hmin.txt file found, skipping single-particle analysis." << std::endl;
        hist_dataset.limit_x({0, 300});
        plots::PlotDataset::quick_plot(hist_dataset, settings::general::output + "distances_limited.png");
    }

    return 0;
}