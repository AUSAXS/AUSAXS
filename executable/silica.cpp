
#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
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
    settings::crystal::h = 200; settings::crystal::k = 200; settings::crystal::l = 0;
    settings::crystal::max_q = 0.5;
    settings::protein::use_effective_charge = false;
    settings::axes::max_distance = 14000;
    settings::crystal::miller_generation_strategy = settings::crystal::MillerGenerationChoice::Reduced;
    settings::crystal::reduced::basis_q = 5;

    //###############################################//
    //###               fourier plot              ###//
    //###############################################//
    // xy slice
    crystal::CrystalScattering cs(crystal);
    auto fourierxy = cs.calculate();
    fourierxy.limit_y(1e-4, std::numeric_limits<double>::max());
    // fourierxy.limit_x(1e-2, 1);
    fourierxy.add_plot_options({{plots::option::legend, "xy"}});

    // complete box
    // settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 10;
    // settings::crystal::max_q = 100;
    // settings::axes::bins = 1000;
    // crystal::CrystalScattering cs2(crystal);
    // auto fourier = cs2.calculate();
    // fourier.limit_y(1e-4, 1e10);
    // fourier.limit_x(1e-2, 1);
    // fourier.add_plot_options({{plots::option::legend, "xyz"}});

    // plots::PlotIntensity plot(fourierxy, style::color::orange);
    // plot.plot(fourier, style::color::blue);
    // plot.save(settings::general::output + "xy_xyz_comparison.png");
    plots::PlotIntensity::quick_plot(fourierxy, settings::general::output + "fourierxy.png");
    // plots::PlotIntensity::quick_plot(fourier, settings::general::output + "fourier.png");

    //###############################################//
    //###               histogram                 ###//
    //###############################################//
    // Protein protein = crystal::Fval::as_protein();
    // protein.set_histogram_manager<hist::HistogramManagerMT>();
    // auto hist = protein.get_histogram();
    // auto debye = hist.calc_debye_scattering_intensity();
    // auto hist_dataset = hist.as_dataset();
    // hist_dataset.add_plot_options({{"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
    // plots::PlotDataset::quick_plot(hist_dataset, settings::general::output + "distances.png");
    // hist_dataset.limit_x({0, 300});
    // plots::PlotDataset::quick_plot(hist_dataset, settings::general::output + "distances_limited.png");
    // plots::PlotIntensity::quick_plot(debye, settings::general::output + "debye.png");

    // // comparison plot
    // debye.add_plot_options({{"color", style::color::red}});
    // plots::PlotIntensity intensity_compare(fourier);
    // intensity_compare.plot(debye, style::color::red);
    // intensity_compare.save(settings::general::output + "debye_comparison.png");

    // // cut histogram if hmin exists
    // if (std::filesystem::exists(settings::general::output + "hmin.txt")) {
    //     std::ifstream input(settings::general::output + "hmin.txt");
    //     if (!input.is_open()) {throw except::io_error("Could not open file " + settings::general::output + "hmin.txt");}
    //     std::string line;
    //     std::getline(input, line);
    //     unsigned int min = std::stoi(line);
    //     std::getline(input, line);
    //     unsigned int max = std::stoi(line);

    //     std::vector<double> xc, yc;
    //     unsigned int i = 0;
    //     while (hist.d[i] < min) {i++;}
    //     while (hist.d[i] < max) {
    //         xc.push_back(hist.d[i]);
    //         yc.push_back(hist.p[i++]);
    //     }

    //     auto model = [] (double x, double x0, double a) {
    //         return a*std::pow(x-x0, 2);
    //     };

    //     auto chi2 = [&xc, &yc, &model] (std::vector<double> p) {
    //         double chi2 = 0;
    //         for (unsigned int i = 0; i < xc.size(); i++) {
    //             chi2 += std::pow(yc[i] - model(xc[i], p[0], p[1]), 2)/std::max(yc[i], 1.0);
    //         }
    //         return chi2;
    //     };

    //     auto mini = mini::dlibMinimizer<mini::type::DLIB_GLOBAL>(chi2, {
    //         mini::Parameter("x0", {-30, xc.front()+30}), 
    //         mini::Parameter("b", {1, 20000})
    //     });
    //     mini.set_max_evals(1000);
    //     auto res = mini.minimize();

    //     SimpleDataset fit;
    //     {
    //         std::vector<double> y(xc.size());
    //         for (unsigned int i = 0; i < xc.size(); i++) {
    //             y[i] = model(xc[i], res[0], res[1]);
    //         }
    //         fit = SimpleDataset(xc, y);
    //         fit.add_plot_options({{"color", style::color::red}});
    //     }
    //     SimpleDataset fit_extrapolated;
    //     {
    //         std::vector<double> xn, yn;
    //         unsigned int i = 0;
    //         while (hist.d[i] <= min) {
    //             xn.push_back(hist.d[i]);
    //             double temp = hist.d[i] < res[0] ? 0 : model(hist.d[i], res[0], res[1]);
    //             yn.push_back(temp);
    //             i++;
    //         }
    //         fit_extrapolated = SimpleDataset(xn, yn);
    //         fit_extrapolated.add_plot_options({{"color", style::color::red}, {"linestyle", style::line::dashed}});
    //     }
    //     auto hist_dataset = hist.as_dataset();
    //     hist_dataset.limit_x(0, 2*max);
    //     hist_dataset.add_plot_options({{"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
    //     plots::PlotDataset plot(hist_dataset);
    //     plot.plot(fit);
    //     plot.plot(fit_extrapolated);
    //     plot.save(settings::general::output + "fit.png");
    //     std::cout << "Fit vals: " << res[0] << ", " << res[1] << std::endl;

    //     i = 0;
    //     std::vector<double> yn;
    //     while (hist.d[i] < min) {
    //         double temp = hist.d[i] < res[0] ? 0 : model(hist.d[i], res[0], res[1]);
    //         yn.push_back(hist.p[i] - temp);
    //         i++;
    //     }
    
    //     Axis new_axis(yn.size(), 0, hist.d[i]);
    //     std::vector<double> empty(yn.size(), 0);
    //     hist::ScatteringHistogram new_hist(yn, empty, empty, yn, new_axis);
    //     auto new_hist_dataset = new_hist.as_dataset();
    //     new_hist_dataset.add_plot_options({{"xlabel", "Distance [$\\AA$]"}, {"ylabel", "Count"}});
    //     plots::PlotDataset::quick_plot(new_hist_dataset, settings::general::output + "distances_single.png");

    //     auto new_intensity = new_hist.calc_debye_scattering_intensity();
    //     plots::PlotIntensity::quick_plot(new_intensity, settings::general::output + "intensity_single.png");
    // }

    return 0;
}