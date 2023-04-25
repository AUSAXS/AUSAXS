#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <plots/all.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>

int main(int argc, char const *argv[]) {
    std::string crystal;
    CLI::App app{"Crystal Scattering"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/crystal/");
    CLI11_PARSE(app, argc, argv);
    settings::general::output += utility::stem(crystal) + "/";

    settings::axes::qmin = 1e-4;
    settings::axes::bins = 1000;
    settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 0;
    settings::crystal::mgc = settings::crystal::MillerGenerationChoice::All;
    settings::crystal::max_q = 100;
    // crystal::CrystalScattering cs(crystal);
    // auto fourier = cs.calculate();
    // fourier.limit_y(1e-4, 1e10);
    // fourier.limit_x(1e-2, 1);
    // fourier.add_plot_options({{plots::option::color, style::color::black}, {plots::option::legend, "fourier"}});
    // plots::PlotIntensity plot(fourier);

    // check if the input is a structure file (e.g. pdb). if so, we also calculate the scattering with the Debye method.
    // if (constants::filetypes::structure.validate(crystal)) {
    //     settings::protein::use_effective_charge = false;
    //     Protein protein(crystal);
    //     auto debye = protein.get_histogram().calc_debye_scattering_intensity();
    //     debye.scale_y(fourier.y(0)/debye.y(0));
    //     debye.add_plot_options(plots::option::draw_markers, {{plots::option::color, style::color::red}, {plots::option::legend, "debye"}});
    //     plot.plot(debye, style::color::red);
    // }
    // plot.save(settings::general::output + "comparison.png");

    settings::axes::qmin = 1e-4;
    settings::axes::bins = 1000;
    settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 0;
    settings::crystal::mgc = settings::crystal::MillerGenerationChoice::Reduced;
    settings::crystal::max_q = 100;
    crystal::CrystalScattering cs(crystal);
    auto fourier = cs.calculate();
    fourier.limit_y(1e-4, 1e10);
    fourier.limit_x(1e-2, 1);
    fourier.add_plot_options({{plots::option::color, style::color::orange}, {plots::option::legend, "xy"}});
    plots::PlotIntensity plot(fourier);

    settings::crystal::l = 100;
    settings::crystal::mgc = settings::crystal::MillerGenerationChoice::All;
    auto fourier2 = cs.calculate();
    fourier2.limit_y(1e-4, 1e10);
    fourier2.limit_x(1e-2, 1);
    fourier2.add_plot_options({{plots::option::color, style::color::blue}, {plots::option::legend, "xyz"}});
    plot.plot(fourier2);
    plot.save(settings::general::output + "comparison.png");

    return 0;
}