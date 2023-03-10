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

    setting::axes::qmin = 0.02;
    setting::axes::bins = 100;
    setting::crystal::h = 100; setting::crystal::k = 100; setting::crystal::l = 10;
    setting::crystal::mgc = setting::crystal::MillerGenerationChoice::Reduced;
    setting::crystal::reduced::hkl_limit = 20;
    crystal::CrystalScattering cs(crystal);
    auto fourier = cs.calculate();

    plots::PlotIntensity::quick_plot(fourier, setting::general::output + "fourier.png");

    return 0;
}