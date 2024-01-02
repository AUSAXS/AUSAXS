#include <CLI/CLI.hpp>

#include <iostream>

#include <plots/All.h>
#include <em/ImageStack.h>
#include <data/Molecule.h>
#include <utility/Utility.h>
#include <fitter/FitReporter.h>
#include <fitter/HydrationFitter.h>
#include <settings/All.h>
#include <constants/Constants.h>
#include <em/manager/ProteinManager.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

int main(int argc, char const *argv[]) {
    settings::hist::weighted_bins = false;
    settings::molecule::use_effective_charge = false;
    settings::em::mass_axis = true;
    settings::em::hydrate = true;
    settings::fit::verbose = true;
    settings::em::alpha_levels = {1, 10};

    io::ExistingFile mfile, mapfile, settings;
    CLI::App app{"Fit an EM map to a SAXS measurement."};
    app.add_option("input-map", mapfile, "Path to the EM map.")->required()->check(CLI::ExistingFile);
    app.add_option("input-exp", mfile, "Path to the SAXS measurement.")->required()->check(CLI::ExistingFile);
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from measurement file.");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from measurement file.");
    app.add_option("--levelmin", settings::em::alpha_levels.min, "Lower limit on the alpha levels to use for the EM map. Note that lowering this limit severely impacts the performance.");
    app.add_option("--levelmax", settings::em::alpha_levels.max, "Upper limit on the alpha levels to use for the EM map. Increasing this limit improves the performance.");
    app.add_option("--frequency", settings::em::sample_frequency, "Sampling frequency of the EM map.");
    app.add_option("--max-iterations", settings::fit::max_iterations, "Maximum number of iterations to perform. This is only approximate.");
    app.add_flag("--mass-axis,!--no-mass-axis", settings::em::mass_axis, "Whether to use a mass axis in place of the threshold axis.");
    app.add_flag("--hydrate,!--no-hydrate", settings::em::hydrate, "Whether to hydrate the protein before fitting.");
    app.add_flag("--fixed-weight,!--no-fixed-weight", settings::em::fixed_weights, "Whether to use a fixed weight for the fit.");
    CLI11_PARSE(app, argc, argv);

    //###################//
    //### PARSE INPUT ###//
    //###################//
    std::cout << "Performing EM fit with map " << mapfile << " and measurement " << mfile << std::endl;
    em::ImageStack map(mapfile); 

    fitter::HydrationFitter fitter(mfile);
    Axis axis(settings::em::alpha_levels, settings::fit::max_iterations);

    std::vector<double> charge_levels;
    double sigma = map.rms();
    Limit cutoff_lims = {map.from_level(axis.min), map.from_level(axis.max)};
    for (double level = sigma; level < cutoff_lims.min; level += sigma) {charge_levels.push_back(level);}
    for (double level = cutoff_lims.min; level < cutoff_lims.max; level += sigma/10) {charge_levels.push_back(level);}
    for (double level = cutoff_lims.max; level < 20*sigma; level += sigma) {charge_levels.push_back(level);}
    map.get_protein_manager()->set_charge_levels(charge_levels);

    unsigned int counter = 0;
    std::vector<unsigned int> atoms;
    for (auto& level : axis.as_vector()) {
        auto protein = map.get_protein(level);
        protein->generate_new_hydration();
        fitter.set_scattering_hist(protein->get_histogram());
        auto res = fitter.fit();
        if (settings::fit::verbose) {
            std::cout << "Step " << utility::print_element(counter++, 4) << ": Evaluated cutoff value " << utility::print_element(level, 8) << " with chi2 " << utility::print_element(res->fval, 8) << std::flush << "\r";
        }
        atoms.push_back(protein->atom_size());
    }
    std::cout << std::endl;
    std::cout << "Range of atoms: [" << atoms.front() << ", " << atoms.back() << "]" << std::endl;
    return 0;
}