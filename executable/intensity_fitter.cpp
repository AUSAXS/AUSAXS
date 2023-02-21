#include <CLI/CLI.hpp>

#include <vector>
#include <string>
#include <iostream>

#include <data/Body.h>
#include <data/Protein.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/Fit.h>
#include <plots/all.h>
#include <fitter/FitReporter.h>

using std::cout, std::endl;

int main(int argc, char const *argv[]) {
    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};

    std::string pdb, mfile, output, settings, placement_strategy = "Radial";
    bool use_existing_hydration = false, fit_excluded_volume = false;
    app.add_option("input_s", pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", output, "Path to save the generated figures at.");
    app.add_option("--reduce,-r", setting::grid::percent_water, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.");
    app.add_option("--grid_width,--gw", setting::grid::width, "The distance between each grid point in Ångström (default: 1). Lower widths increase the precision.");
    app.add_option("--bin_width,--bw", setting::axes::distance_bin_width, "Bin width for the distance histograms. Default: 1.");
    app.add_option("--placement_strategy,--ps", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.");
    app.add_option("--radius_a,--ra", setting::grid::ra, "Radius of the protein atoms.");
    app.add_option("--radius_h,--rh", setting::grid::rh, "Radius of the hydration atoms.");
    app.add_option("--qmin", setting::axes::qmin, "Lower limit on used q values from the measurement file.");
    app.add_option("--qmax", setting::axes::qmax, "Upper limit on used q values from the measurement file.");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile);
    app.add_flag("--center,!--no-center", setting::protein::center, "Decides whether the protein will be centered. Default: true.");
    app.add_flag("--effective-charge,!--no-effective-charge", setting::protein::use_effective_charge, "Decides whether the effective atomic charge will be used. Default: true.");
    app.add_flag("--use-existing-hydration,!--no-use-existing-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used. Default: false.");
    app.add_flag("--fit-excluded-volume,!--no-fit-excluded-volume", fit_excluded_volume, "Decides whether the excluded volume will be fitted. Default: true.");
    CLI11_PARSE(app, argc, argv);

    //####################//
    //### PARSE INPUT ###//
    //####################//
    // if a settings file was provided
    if (p_settings->count() != 0) {
        setting::read(settings);        // read it
        CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
    } else {                            // otherwise check if there is a settings file in the same directory
        if (setting::discover(std::filesystem::path(mfile).parent_path().string())) {
            CLI11_PARSE(app, argc, argv);
        }
    }

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
    if (placement_strategy == "Radial") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::RadialStrategy;}
    else if (placement_strategy == "Axes") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::AxesStrategy;}
    else if (placement_strategy == "Jan") {setting::grid::placement_strategy = setting::grid::PlacementStrategy::JanStrategy;}

    if (output.empty()) {
        output = "figures/intensity_fitter/" + utility::stem(mfile) + "/";
    }

    //######################//
    //### ACTUAL PROGRAM ###//
    //######################//
    Protein protein(pdb);
    if (!use_existing_hydration || protein.hydration_atoms.empty()) {
        protein.generate_new_hydration();
    }
    // hist::ScatteringHistogram h = protein.get_histogram();
    // plots::PlotDistance::quick_plot(h, output + "p(r)." + setting::plot::format);

    std::shared_ptr<fitter::HydrationFitter> fitter;
    if (fit_excluded_volume) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(mfile, protein);}
    else {fitter = std::make_shared<fitter::HydrationFitter>(mfile, protein.get_histogram());}
    std::shared_ptr<Fit> result = fitter->fit();
    FitReporter::report(result);
    FitReporter::save(result, output + "report.txt");

    // save fit
    auto fit = fitter->get_model_dataset();
    auto data = fitter->get_dataset();
    fit.save(output + "fit.fit");
    data.save(output + utility::stem(mfile) + ".dat");

    // Fit plot
    plots::PlotIntensityFit::quick_plot(result, output + "fit." + setting::plot::format);
    plots::PlotIntensityFitResiduals::quick_plot(result, output + "residuals." + setting::plot::format);

    // calculate rhoM
    double rhoM = protein.absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
    std::cout << "RhoM is " << rhoM << " g/cm³" << std::endl;

    // std::vector<double> q;
    // for (double qv = 0; qv < 1; qv+=0.01) {q.push_back(qv);}

    // SimpleDataset data = fitter.get_model_dataset(q);
    // plots::PlotDataset plot_d(data);
    // plot_d.save(output + "fitted_model." + setting::figures::format);

    // // double intercept = result->params["b"];
    // double I0 = fitter.get_intercept();
    // double DrhoV2 = std::pow(protein.get_relative_charge(), 2);
    // double re2 = pow(constants::radius::electron*constants::unit::cm, 2);
    // double m = protein.get_absolute_mass()*constants::unit::mg;

    // cout << "concentration is: " << I0*m/(DrhoV2*re2) << endl;

    return 0;
}