#include <CLI/CLI.hpp>

#include <plots/all.h>
#include <math/SimpleLeastSquares.h>

#include <string>

int main(int argc, char** argv) {
    CLI::App app{"Plotting program. Plots the specified dataset and all .fit files in its directory."};

    std::string path, out, settings;
    app.add_option("mfile", path, "Path to the measured SAXS curve.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", out, "Path to save the generated figures at.");
    app.add_option("--qlow", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qhigh", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    CLI11_PARSE(app, argc, argv);

    if (out.empty()) {out = std::filesystem::path(path).parent_path().string() + "/";}

    // if a settings file was provided
    if (p_settings->count() != 0) {
        setting::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        setting::discover(out);
    }

    // search for .fit files in the same directory as the measurement file
    std::vector<std::string> fitfiles;
    for (auto& p: std::filesystem::directory_iterator(std::filesystem::path(path).parent_path())) {
        if (p.path().extension() == ".fit") {
            fitfiles.push_back(p.path().string());
        }
    }

    // make the plot
    SimpleDataset data(path);
    // data.normalize(1);
    data.add_plot_options("error", {{"color", kBlack}, {"xlabel", "q"}, {"ylabel", "I(q)"}, {"s", 0.3}});

    plots::PlotIntensity plot(data);
    std::vector<ELineStyle> styles = {kSolid};//, kDashed, kDotted, kDashDotted};
    std::vector<int> colors = {kOrange+2, kAzure+1, kGreen+1, kYellow+1, kViolet+1, kRed-4, kMagenta-7, kSpring+10};
    std::vector<std::string> names = {"orange", "blue", "green", "yellow", "purple", "red", "magenta", "light green"};
    for (unsigned int i = 0; i < fitfiles.size(); i++) {
        Dataset2D fit(fitfiles[i]);
        if ((fitfiles[i].find("foxs") != std::string::npos) || (fitfiles[i].find("crysol") != std::string::npos)) {fit.y() = fit.xerr();}
        if (fitfiles[i].find("waxsis") != std::string::npos) {
            std::cout << "\tWAXSiS fits are not handled correctly!" << std::endl;
            fit.normalize(data.y(0));
        }
        // fit.normalize(1);
        fit.add_plot_options("lines", {{"color", colors[i % colors.size()]}, {"linestyle", styles[i % styles.size()]}});
        plot.plot_intensity(fit);

        // calculate chi2
        double chi2 = 0;
        for (unsigned int j = 0; j < data.x().size(); j++) {
            chi2 += std::pow((data.y()[j] - fit.y()[j]) / data.yerr()[j], 2);
        }
        std::cout << "\tChi2 for file is " << chi2 << std::endl;
        std::cout << "\tPlotted file in " + names[i % names.size()] + "." << std::endl;
    }
    plot.save(out + "plot." + setting::figures::format);

    return 0;
}