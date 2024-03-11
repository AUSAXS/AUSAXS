#include <CLI/CLI.hpp>

#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <plots/All.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <fitter/FitReporter.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>

int main(int argc, char const *argv[]) {
    io::ExistingFile crystal;
    CLI::App app{"Crystal Scattering"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/crystal_compare/");
    CLI11_PARSE(app, argc, argv);
    settings::general::output += crystal.stem() + "/";
    settings::crystal::detail::use_checkpointing = false;
    settings::general::detail::job_size = 10000;

    settings::molecule::use_effective_charge = false;
    data::Molecule protein(crystal);
    auto debye = protein.get_histogram()->debye_transform().as_dataset();

    settings::axes::qmin = 1e-4;
    settings::crystal::h = 200; settings::crystal::k = 200; settings::crystal::l = 200;

    // compare with fibonnaci
    // settings::crystal::mgc = settings::crystal::MillerGenerationChoice::Reduced;
    // settings::crystal::max_q = 100;
    // for (unsigned int i = 3; i < 11; i++) {
    //     settings::crystal::reduced::basis_q = i;
    //     crystal::CrystalScattering cs(crystal);
    //     auto fourier = cs.calculate();
    //     fourier.limit_y(1e-4, 1e10);
    //     fourier.scale_y(debye.y(0)/fourier.y(0));
    //     fourier.save(settings::general::output + "reduced_" + std::to_string(i) + ".fit");

    //     fourier.add_plot_options({{plots::option::color, style::color::black}, {plots::option::legend, "max_basis_hkl = " + std::to_string(i)}});
    //     debye.add_plot_options(plots::option::draw_markers, {{plots::option::color, style::color::red}});
    //     plots::PlotIntensity plot(fourier);
    //     plot.plot(debye, style::color::red);
    //     plot.save(settings::general::output + "reduced_" + std::to_string(i) + ".png");
    // }
    // debye.save(settings::general::output + "debye.dat");

    // compare with fibonnaci
    settings::crystal::miller_generation_strategy = settings::crystal::MillerGenerationChoice::Reduced;
    settings::crystal::max_q = 100;
    std::vector<std::vector<double>> lsqs(11, std::vector<double>(20, 0));
    for (unsigned int i = 3; i < 11; i++) {
        for (unsigned int j = 1; j < 20; j++) {
            settings::crystal::grid_expansion = 1 + 0.2*j;

            settings::crystal::reduced::basis_q = i;
            crystal::CrystalScattering cs(crystal);
            auto fourier = cs.calculate();
            fourier.limit_y(1e-4, 1e10);
            fourier.scale_y(debye.y(0)/fourier.y(0));
            fourier.save(settings::general::output + std::to_string(settings::crystal::grid_expansion) + "/reduced_" + std::to_string(i) + ".fit");

            plots::PlotIntensity plot(debye, plots::PlotOptions(plots::option::draw_markers, {{plots::option::color, style::color::red}, {plots::option::legend, "debye"}}));
            plot.plot(fourier, plots::PlotOptions({{plots::option::color, style::color::black}, {plots::option::legend, "max_basis_hkl = " + std::to_string(i)}}));
            plot.save(settings::general::output + std::to_string(settings::crystal::grid_expansion) + "/reduced_" + std::to_string(i) + ".png");

            // calculate square difference
            double lsq = 0;
            for (unsigned int k = 0; k < fourier.size(); k++) {
                lsq += std::pow(fourier.y(k) - debye.y(k), 2);
            }
            lsqs[i][j] = lsq;
        }
    }
    debye.save(settings::general::output + "debye.dat");

    return 0;
}