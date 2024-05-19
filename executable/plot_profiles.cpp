#include "fitter/HydrationFitter.h"
#include "fitter/FitReporter.h"
#include <CLI/CLI.hpp>

#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <grid/Grid.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <dataset/SimpleDataset.h>
#include <settings/All.h>
#include <plots/All.h>

#include <sstream>

int main(int argc, char const *argv[]) {
    io::ExistingFile pdb;
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/saxs_fitter/")->group("General options");
    CLI11_PARSE(app, argc, argv);

    //### GENERATE INTERNAL PLOT ###//
    settings::general::output += pdb.stem() + "/";
    settings::axes::qmin = 1e-4;
    settings::axes::qmax = 1;
    settings::molecule::use_effective_charge = false;
    settings::general::threads = 1;

    SimpleDataset ausaxs_crysol_aa, ausaxs_crysol_xx, ausaxs_foxs_aa, ausaxs_foxs_xx, ausaxs_aa, ausaxs_xx;
    {   // crysol mimic
        settings::molecule::implicit_hydrogens = true;
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::CrysolManager;
        data::Molecule molecule(pdb);

        auto hist = molecule.get_histogram();
        ausaxs_crysol_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_crysol_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }
    {   // foxs mimic
        settings::molecule::implicit_hydrogens = true;
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::FoXSManager;
        data::Molecule molecule(pdb);

        auto hist = molecule.get_histogram();
        ausaxs_foxs_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_foxs_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }
    {   // ausaxs
        settings::molecule::implicit_hydrogens = true;
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit;
        data::Molecule molecule(pdb);

        auto hist = molecule.get_histogram();
        ausaxs_aa = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_aa().as_dataset();
        ausaxs_xx = static_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())->get_profile_xx().as_dataset();
    }

    //### CHECK FOR PRESENCE OF EXTERNAL DATA IN OUTPUT FOLDER ###//
    io::File crysol(settings::general::output + "crysol.dat");
    io::File foxs(settings::general::output + "foxs_xx.dat");

    if (!crysol.exists() || !foxs.exists()) {
        std::cout << "External data not found in output folder. Skipping external comparison." << std::endl;
        return 0;
    }

    // generate plot
    SimpleDataset crysol_data;
    {
        Dataset tmp(crysol);
        crysol_data = SimpleDataset(tmp.x(), tmp.col(3));
    }
    auto foxs_data = SimpleDataset(foxs);

    plots::PlotIntensity()
        .plot(crysol_data, plots::PlotOptions({{"legend", "CRYSOL"}, {"xlabel", "q (Å⁻¹)"}, {"ylabel", "I(q)"}, {"color", style::color::cyan}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"xrange", Limit(1e-2, 1)}}))
        .plot(foxs_data, plots::PlotOptions({{"legend", "FoXS"}, {"color", style::color::orange}}))
        .plot(ausaxs_crysol_xx, plots::PlotOptions({{"legend", "AUSAXS$_{crysol}$"}, {"color", style::color::blue}, {"linestyle", style::line::dashed}}))
        .plot(ausaxs_foxs_xx, plots::PlotOptions({{"legend", "AUSAXS$_{foxs}$"}, {"color", style::color::red}, {"linestyle", style::line::dashed}}))
        .plot(ausaxs_xx, plots::PlotOptions({{"legend", "AUSAXS"}, {"color", style::color::green}, {"linestyle", style::line::dashed}}))
    .save(settings::general::output + "profiles_xx.png");

    SimpleDataset crysol_diff = crysol_data;
    for (size_t i = 0; i < crysol_diff.size(); ++i) {
        crysol_diff.y(i) /= ausaxs_crysol_xx.interpolate_x(crysol_diff.x(i), 1);
    }

    SimpleDataset foxs_diff = foxs_data;
    for (size_t i = 0; i < foxs_diff.size(); ++i) {
        foxs_diff.y(i) /= ausaxs_foxs_xx.interpolate_x(foxs_diff.x(i), 1);
    }

    plots::PlotDataset()
        .plot(crysol_diff, plots::PlotOptions({{"xlabel", "q (Å⁻¹)"}, {"ylabel", "CRYSOL / AUSAXS"}, {"title", pdb.stem() + " $I_{xx}$ profiles"}, {"yrange", Limit(0.5, 1.5)}, {"xrange", Limit(1e-2, 1)}, {"color", style::color::cyan}}))
        .plot(foxs_diff, plots::PlotOptions({{"ylabel", "FoXS / AUSAXS"}, {"color", style::color::orange}}))
    .save(settings::general::output + "profiles_xx_diff.png");
}