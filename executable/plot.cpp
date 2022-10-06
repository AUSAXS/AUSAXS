#include <CLI/CLI.hpp>

#include <plots/all.h>

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
    data.add_plot_options("error", {{"color", kBlack}, {"xlabel", "q"}, {"ylabel", "I(q)"}});

    plots::PlotIntensity plot(data);
    std::vector<ELineStyle> styles = {kSolid, kDashed, kDotted, kDashDotted};
    std::vector<int> colors = {kOrange+2, kAzure+1, kGreen+1, kYellow+1, kViolet+1, kRed-4, kMagenta-7, kSpring+10};
    for (unsigned int i = 0; i < fitfiles.size(); i++) {
        SimpleDataset fit(fitfiles[i]);
        fit.add_plot_options("lines", {{"color", colors[i % colors.size()]}, {"linestyle", styles[i % styles.size()]}});
        plot.plot_intensity(fit);
    }
    plot.save(out + "plot." + setting::figures::format);

    return 0;
}