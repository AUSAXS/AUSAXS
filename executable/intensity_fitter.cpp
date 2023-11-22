#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/FitReporter.h>
#include <fitter/Fit.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string s_pdb, s_mfile, s_settings, placement_strategy = "radial", histogram_manager = "hmmt"; // not using partial histograms has a slightly smaller overhead
    bool use_existing_hydration = false, fit_excluded_volume = false;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/intensity_fitter/")->group("General options");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads)->group("General options");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from the measurement file.")->default_val(settings::axes::qmax)->group("General options");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from the measurement file.")->default_val(settings::axes::qmin)->group("General options");
    auto p_settings = app.add_option("-s,--settings", s_settings, "Path to the settings file.")->check(CLI::ExistingFile)->group("General options");

    app.add_flag("--center,!--no-center", settings::molecule::center, "Decides whether the protein will be centered.")->default_val(settings::molecule::center)->group("Protein options");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::molecule::use_effective_charge, "Decides whether the effective atomic charge will be used.")->default_val(settings::molecule::use_effective_charge)->group("Protein options");
    app.add_flag("--use-existing-hydration,!--no-use-existing-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used.")->default_val(use_existing_hydration)->group("Protein options");
    app.add_flag("--fit-excluded-volume,!--no-fit-excluded-volume", fit_excluded_volume, "Decides whether the excluded volume will be fitted.")->default_val(fit_excluded_volume)->group("Protein options");

    app.add_option("--reduce,-r", settings::grid::water_scaling, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")->default_val(settings::grid::water_scaling)->group("Advanced options");
    app.add_option("--grid_width,--gw", settings::grid::width, "The distance between each grid point in Ångström. Lower widths increase the precision.")->default_val(settings::grid::width)->group("Advanced options");
    app.add_option("--placement_strategy,--ps", placement_strategy, "The placement strategy to use. Options: Radial, Axes, Jan.")->default_val(placement_strategy)->group("Advanced options");
    app.add_option("--exv_radius,--er", settings::grid::exv_radius, "The radius of the excluded volume sphere used for the grid-based excluded volume calculations in Ångström.")->default_val(settings::grid::exv_radius)->group("Advanced options");
    auto p_hm = app.add_option("--histogram-manager,--hm", histogram_manager, "The histogram manager to use. Options: HM, HMMT, HMMTFF, PHM, PHMMT, PHMMTFF.")->group("Advanced options");

    app.add_flag("--foxs", settings::hist::use_foxs_method, "Decides whether the FOXS method will be used.")->default_val(settings::hist::use_foxs_method)->group("Hidden");
    app.add_flag("--weighted_bins", settings::hist::weighted_bins, "Decides whether the weighted bins will be used.")->default_val(settings::hist::weighted_bins)->group("Hidden");
    app.add_option("--rvol", settings::grid::rvol, "The radius of the excluded volume sphere around each atom.")->default_val(settings::grid::rvol)->group("Hidden");
    app.add_flag("--save_exv", settings::grid::save_exv, "Decides whether the excluded volume will be saved.")->default_val(settings::grid::save_exv)->group("Hidden");
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

    if (p_hm->count() == 0) {
        if (settings::general::threads == 1) {settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManager;}
        else {settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;}
    } else {
        settings::detail::parse_option("histogram_manager", {histogram_manager});
        if (fit_excluded_volume) {
            switch(settings::hist::histogram_manager) {
                case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
                case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
                case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: break;
                default: throw except::invalid_argument("The histogram manager must be either HMMTFF, HMMTFFExplicit or HMMTFFGrid when fitting the excluded volume.");
            }
        }
    }
    switch (settings::hist::histogram_manager) {
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit: settings::molecule::use_effective_charge = false;
        default: break;
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

    //######################//
    //### ACTUAL PROGRAM ###//
    //######################//
    data::Molecule protein(pdb);
    if (!use_existing_hydration || protein.water_size() == 0) {
        protein.generate_new_hydration();
    }

    std::shared_ptr<fitter::HydrationFitter> fitter;
    if (fit_excluded_volume) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(mfile, protein.get_histogram());}
    else {fitter = std::make_shared<fitter::HydrationFitter>(mfile, protein.get_histogram());}
    std::shared_ptr<fitter::Fit> result = fitter->fit();
    fitter::FitReporter::report(result);
    fitter::FitReporter::save(result, settings::general::output + "report.txt");

    plots::PlotDistance::quick_plot(fitter->get_scattering_hist(), settings::general::output + "p(r)." + settings::plots::format);
    plots::PlotProfiles::quick_plot(fitter->get_scattering_hist(), settings::general::output + "profiles." + settings::plots::format);

    // save fit
    fitter->get_model_dataset().save(settings::general::output + "fit.fit");
    fitter->get_dataset().save(settings::general::output + mfile.stem() + ".scat");

    // calculate rhoM
    double rhoM = protein.absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
    std::cout << "RhoM is " << rhoM << " g/cm³" << std::endl;

    protein.save(settings::general::output + "model.pdb");

    // auto h = fitter->get_scattering_hist();
    // h->get_profile_aa().as_dataset().save(settings::general::output + "ausaxs_aa.dat");
    // h->get_profile_aw().as_dataset().save(settings::general::output + "ausaxs_aw.dat");
    // h->get_profile_ww().as_dataset().save(settings::general::output + "ausaxs_ww.dat");
    // if (auto cast = dynamic_cast<hist::ICompositeDistanceHistogramExv*>(h.get())) {
    //     if (fit_excluded_volume) {cast->apply_excluded_volume_scaling_factor(result->get_parameter("d"));}
    //     cast->get_profile_ax().as_dataset().save(settings::general::output + "ausaxs_ax.dat");
    //     cast->get_profile_wx().as_dataset().save(settings::general::output + "ausaxs_wx.dat");
    //     cast->get_profile_xx().as_dataset().save(settings::general::output + "ausaxs_xx.dat");
    // }

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