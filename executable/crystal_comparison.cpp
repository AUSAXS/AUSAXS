#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <plots/all.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <fitter/FitReporter.h>
#include <crystal/Settings.h>

int main(int argc, char const *argv[]) {
    std::string crystal;
    CLI::App app{"Crystal Scattering"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", setting::general::output, "Path to save the generated figures at.")->default_val("output/crystal_compare/");
    CLI11_PARSE(app, argc, argv);
    setting::general::output += utility::stem(crystal) + "/";

    setting::protein::use_effective_charge = false;
    Protein protein(crystal);
    auto debye = protein.get_histogram().calc_debye_scattering_intensity();

    setting::axes::qmin = 1e-4;
    setting::axes::bins = 1000;
    setting::crystal::h = 200; setting::crystal::k = 200; setting::crystal::l = 200;

    // compare with fibonnaci
    // setting::crystal::mgc = setting::crystal::MillerGenerationChoice::Reduced;
    // setting::crystal::max_q = 100;
    // for (unsigned int i = 3; i < 11; i++) {
    //     setting::crystal::reduced::basis_q = i;
    //     crystal::CrystalScattering cs(crystal);
    //     auto fourier = cs.calculate();
    //     fourier.limit_y(1e-4, 1e10);
    //     fourier.scale_y(debye.y(0)/fourier.y(0));
    //     fourier.save(setting::general::output + "reduced_" + std::to_string(i) + ".fit");

    //     fourier.add_plot_options({{plots::option::color, style::color::black}, {plots::option::legend, "max_basis_hkl = " + std::to_string(i)}});
    //     debye.add_plot_options(plots::option::draw_markers, {{plots::option::color, style::color::red}});
    //     plots::PlotIntensity plot(fourier);
    //     plot.plot(debye, style::color::red);
    //     plot.save(setting::general::output + "reduced_" + std::to_string(i) + ".png");
    // }
    // debye.save(setting::general::output + "debye.dat");

    // compare with fibonnaci
    setting::crystal::mgc = setting::crystal::MillerGenerationChoice::Reduced;
    setting::crystal::max_q = 100;
    std::vector<std::vector<double>> lsqs(11, std::vector<double>(20, 0));
    for (unsigned int i = 3; i < 11; i++) {
        for (unsigned int j = 1; j < 20; j++) {
            setting::crystal::grid_expansion = 1 + 0.2*j;

            setting::crystal::reduced::basis_q = i;
            crystal::CrystalScattering cs(crystal);
            auto fourier = cs.calculate();
            fourier.limit_y(1e-4, 1e10);
            fourier.scale_y(debye.y(0)/fourier.y(0));
            fourier.save(setting::general::output + std::to_string(setting::crystal::grid_expansion) + "/reduced_" + std::to_string(i) + ".fit");

            // calculate square difference
            double lsq = 0;
            for (unsigned int k = 0; k < fourier.size(); k++) {
                lsq += std::pow(fourier.y(k) - debye.y(k), 2);
            }
            lsqs[i][j] = lsq;
        }
    }
    debye.save(setting::general::output + "debye.dat");

    return 0;
}