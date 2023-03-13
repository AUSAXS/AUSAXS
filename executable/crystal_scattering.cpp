#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <plots/all.h>

int main(int argc, char const *argv[]) {
    std::string crystal;
    CLI::App app{"Crystal Scattering"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", setting::general::output, "Path to save the generated figures at.")->default_val("output/crystal/");
    CLI11_PARSE(app, argc, argv);
    setting::general::output += utility::stem(crystal) + "/";

    setting::protein::use_effective_charge = false;
    Protein protein(crystal);
    auto debye = protein.get_histogram().calc_debye_scattering_intensity();

    setting::axes::qmin = 1e-4;
    setting::axes::bins = 1000;
    setting::crystal::h = 100; setting::crystal::k = 100; setting::crystal::l = 100;
    setting::crystal::mgc = setting::crystal::MillerGenerationChoice::All;
    setting::crystal::reduced::hkl_limit = 20;
    crystal::CrystalScattering cs(crystal);
    auto fourier = cs.calculate();
    fourier.limit_y(1e-4, 1e10);
    fourier.scale_y(debye.y(0)/fourier.y(0));

    fourier.add_plot_options({{plots::option::color, style::color::black}, {plots::option::legend, "fourier"}});
    debye.add_plot_options(plots::option::draw_markers, {{plots::option::color, style::color::red}, {plots::option::legend, "debye"}});
    plots::PlotIntensity plot(fourier);
    plot.plot(debye, style::color::red);
    plot.save(setting::general::output + "comparison.png");

    return 0;
}