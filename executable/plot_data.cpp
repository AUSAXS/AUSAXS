#include <CLI/CLI.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <utility/SimpleDataset.h>
#include <utility/Settings.h>
#include <plots/all.h>

int main(int argc, char const *argv[]) {
    CLI::App app{"Plot a SAXS dataset."};

    std::string input_measurement, output = "figures/plot_data/", placement_strategy, settings;
    app.add_option("input_m", input_measurement, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--qlow", setting::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qhigh", setting::axes::qmax, "Upper limit on used q values from measurement file.");
    CLI11_PARSE(app, argc, argv);

    if (p_settings->count() != 0) {
        setting::read(settings);
        CLI11_PARSE(app, argc, argv);
    } else {
        std::string path = std::filesystem::path(input_measurement).parent_path().string();
        if (std::filesystem::exists(path + "/settings.txt")) {
            std::cout << "Using discovered settings file at " << path << "/settings.txt" << std::endl;
            setting::read(path + "/settings.txt");
        }
    }

    SimpleDataset data(input_measurement);
    data.limit_x(setting::axes::qmin, setting::axes::qmax);
    auto span = data.span_x();
    data.add_plot_options("errors", {{"color", kOrange+2}, {"xlimits", span}});

    plots::PlotIntensity plot(data);
    plot.save(output + "plot.pdf");
}