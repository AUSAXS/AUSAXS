#include <CLI/CLI.hpp>

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <crystal/CrystalScattering.h>
#include <hydrate/GridReader.h>
#include <math/SimpleLeastSquares.h>
#include <math/CubicSpline.h>
#include <fitter/FitReporter.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <settings/All.h>
#include <plots/All.h>

int main(int argc, char const *argv[]) {
    io::File crystal;
    CLI::App app{"Crystal Scattering"};
    app.add_option("input", crystal, "File containing the crystal data.")->required();
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/crystal/");
    CLI11_PARSE(app, argc, argv);
    settings::general::output += crystal.stem() + "/";

    settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 100;
    settings::crystal::miller_generation_strategy = settings::crystal::MillerGenerationChoice::Reduced;
    settings::crystal::grid_expansion = 3;
    settings::crystal::max_q = 0.5;
    crystal::CrystalScattering cs(crystal);
    auto fourier = cs.rotational_average(100);
    plots::PlotIntensity plot(fourier, plots::PlotOptions({{plots::option::color, style::color::black}, {plots::option::legend, "fourier"}}));

    // check if the input is a structure file (e.g. pdb). if so, we also calculate the scattering with the Debye method.
    if (constants::filetypes::structure.validate(crystal)) {
        settings::molecule::use_effective_charge = false;
        settings::molecule::implicit_hydrogens = false;
        data::Molecule protein(crystal);
        auto debye = protein.get_histogram()->debye_transform().as_dataset();

        // fit the debye function to the fourier function
        fourier = fourier.interpolate(debye.x());
        fitter::SimpleLeastSquares fitter(debye.y(), fourier.y());
        auto res = fitter.fit();
        std::transform(debye.y().begin(), debye.y().end(), debye.y().begin(), [a = res->get_parameter("a"), b = res->get_parameter("b")](double y) {return y*a + b;});
        fitter::FitReporter::report(res.get());
        debye.scale_y(fourier.y(0)/debye.y(0));
        plot.plot(debye, plots::PlotOptions(plots::option::draw_markers, {{plots::option::color, style::color::red}, {plots::option::legend, "debye"}}));
    }
    plot.save(settings::general::output + "comparison.png");

    //#########################################//
    //### CUSTOM SECTION FOR SILICA PROJECT ###//
    //#########################################//
    if (false) {
        settings::crystal::grid_expansion = 1;
        settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 0;
        settings::crystal::max_q = 0.2;
        crystal::CrystalScattering cs1(crystal);
        auto fourier = cs1.calculate();
        fourier.limit_y(1e-4, 1e10);
        fourier.limit_x(1e-2, 1);
        plots::PlotIntensity plot(fourier, plots::PlotOptions({{plots::option::color, style::color::orange}, {plots::option::legend, "xy"}}));
        plot.save(settings::general::output + "fourier.png");

        // settings::axes::qmin = 1e-4;
        // settings::axes::bins = 100;
        // settings::crystal::h = 100; settings::crystal::k = 0; settings::crystal::l = 10;
        // crystal::CrystalScattering cs2(crystal);
        // auto fourier2 = cs2.calculate();
        // fourier2.limit_y(1e-4, 1e10);
        // fourier2.limit_x(1e-2, 1);
        // fourier2.add_plot_options({{plots::option::legend, "xz"}});
        // plot.plot(fourier2, style::color::purple);

        // settings::axes::qmin = 1e-4;
        // settings::axes::bins = 100;
        // settings::crystal::h = 0; settings::crystal::k = 100; settings::crystal::l = 10;
        // crystal::CrystalScattering cs3(crystal);
        // auto fourier3 = cs3.calculate();
        // fourier3.limit_y(1e-4, 1e10);
        // fourier3.limit_x(1e-2, 1);
        // fourier3.add_plot_options({{plots::option::legend, "yz"}});
        // plot.plot(fourier3, style::color::green);

        settings::crystal::h = 100; settings::crystal::k = 100; settings::crystal::l = 10;
        settings::crystal::max_q = 100;
        crystal::CrystalScattering cs4(crystal);
        auto fourier4 = cs4.calculate();
        fourier4.limit_y(1e-4, 1e10);
        fourier4.limit_x(1e-2, 1);
        plot.plot(fourier4, plots::PlotOptions({{plots::option::color, style::color::blue}, {plots::option::legend, "xyz"}}));
        plot.save(settings::general::output + "comparison.png");
    }

    return 0;
}