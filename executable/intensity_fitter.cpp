#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/Water.h>
#include <data/Protein.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/FitReporter.h>
#include <fitter/Fit.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <utility/Constants.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <hist/CompositeDistanceHistogram.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::string s_pdb, s_mfile, s_settings, placement_strategy = "Radial";
    bool use_existing_hydration = false, fit_excluded_volume = false;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT; // smaller overhead

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/intensity_fitter/");
    app.add_option("--reduce,-r", settings::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")->default_val(settings::grid::percent_water);
    app.add_option("--grid_width,--gw", settings::grid::width, "The distance between each grid point in Ångström. Lower widths increase the precision.")->default_val(settings::grid::width);
    app.add_option("--bin_width,--bw", settings::axes::distance_bin_width, "Bin width for the distance histograms.")->default_val(settings::axes::distance_bin_width);
    app.add_option("--placement_strategy,--ps", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.")->default_val(placement_strategy);
    app.add_option("--radius_a,--ra", settings::grid::ra, "Radius of the protein atoms.")->default_val(settings::grid::ra);
    app.add_option("--radius_h,--rh", settings::grid::rh, "Radius of the hydration atoms.")->default_val(settings::grid::rh);
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from the measurement file.")->default_val(settings::axes::qmin);
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from the measurement file.")->default_val(settings::axes::qmax);
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads);
    auto p_settings = app.add_option("-s,--settings", s_settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_flag("--center,!--no-center", settings::protein::center, "Decides whether the protein will be centered.")->default_val(settings::protein::center);
    app.add_flag("--effective-charge,!--no-effective-charge", settings::protein::use_effective_charge, "Decides whether the effective atomic charge will be used.")->default_val(settings::protein::use_effective_charge);
    app.add_flag("--use-existing-hydration,!--no-use-existing-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used.")->default_val(use_existing_hydration);
    app.add_flag("--fit-excluded-volume,!--no-fit-excluded-volume", fit_excluded_volume, "Decides whether the excluded volume will be fitted.")->default_val(fit_excluded_volume);
    CLI11_PARSE(app, argc, argv);

    //###################//
    //### PARSE INPUT ###//
    //###################//
    io::ExistingFile pdb(s_pdb), mfile(s_mfile), settings(s_settings);
    settings::general::output += mfile.stem() + "/";

    // if a settings file was provided
    if (p_settings->count() != 0) {
        settings::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (settings::discover(mfile.directory())) {
            CLI11_PARSE(app, argc, argv);
        }
    }

    if (settings::general::threads == 1) {settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManager;}
    else {settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;}

    // validate input
    if (!constants::filetypes::structure.validate(pdb)) {
        // check if the two inputs are switched
        if (constants::filetypes::structure.validate(mfile)) {
            // if so, silently swap them and proceed
            std::swap(pdb, mfile);
        } else {
            throw except::invalid_argument("Unknown PDB extensions: " + pdb + " and " + mfile);
        }
    }
    if (!constants::filetypes::saxs_data.validate(mfile)) {
        throw except::invalid_argument("Unknown SAXS data extension: " + mfile);
    }

    // parse strategy
    if (placement_strategy == "Radial") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {settings::grid::placement_strategy = settings::grid::PlacementStrategy::JanStrategy;}

    //######################//
    //### ACTUAL PROGRAM ###//
    //######################//
    Protein protein(pdb);
    if (!use_existing_hydration || protein.water_size() == 0) {
        protein.generate_new_hydration();
    }
    // hist::ScatteringHistogram h = protein.get_histogram();
    // plots::PlotDistance::quick_plot(h, output + "p(r)." + settings::plot::format);

    std::shared_ptr<fitter::HydrationFitter> fitter;
    if (fit_excluded_volume) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(mfile, protein);}
    else {fitter = std::make_shared<fitter::HydrationFitter>(mfile, protein.get_histogram());}
    std::shared_ptr<fitter::Fit> result = fitter->fit();
    fitter::FitReporter::report(result);
    fitter::FitReporter::save(result, settings::general::output + "report.txt");

    // save fit
    auto fit = fitter->get_model_dataset();
    auto data = fitter->get_dataset();
    fit.save(settings::general::output + "fit.fit");
    data.save(settings::general::output + mfile.stem() + ".dat");

    // Fit plot
    plots::PlotIntensityFit::quick_plot(result, settings::general::output + "fit." + settings::plots::format);
    plots::PlotIntensityFitResiduals::quick_plot(result, settings::general::output + "residuals." + settings::plots::format);

    // calculate rhoM
    double rhoM = protein.absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
    std::cout << "RhoM is " << rhoM << " g/cm³" << std::endl;

    // std::vector<double> q;
    // for (double qv = 0; qv < 1; qv+=0.01) {q.push_back(qv);}

    // SimpleDataset data = fitter.get_model_dataset(q);
    // plots::PlotDataset plot_d(data);
    // plot_d.save(output + "fitted_model." + settings::figures::format);

    // // double intercept = result->params["b"];
    // double I0 = fitter.get_intercept();
    // double DrhoV2 = std::pow(protein.get_relative_charge(), 2);
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2);
    // double m = protein.get_absolute_mass()*constants::unit::mg;

    // cout << "concentration is: " << I0*m/(DrhoV2*re2) << endl;

    return 0;
}