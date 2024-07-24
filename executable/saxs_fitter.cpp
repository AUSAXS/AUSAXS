#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <fitter/HydrationFitter.h>
#include <fitter/ExcludedVolumeFitter.h>
#include <fitter/FitReporter.h>
#include <plots/All.h>
#include <settings/All.h>
#include <io/ExistingFile.h>
#include <constants/Constants.h>
#include <mini/detail/Evaluation.h>
#include <mini/detail/FittedParameter.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Console.h>

#include <vector>
#include <string>
#include <iostream>

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::string s_pdb, s_mfile, s_settings, histogram_manager = "hmmt"; // not using partial histograms has a slightly smaller overhead
    bool use_existing_hydration = false;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.add_option("input_s", s_pdb, "Path to the structure file.")->required()->check(CLI::ExistingFile);
    app.add_option("input_m", s_mfile, "Path to the measured data.")->required()->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Path to save the generated figures at.")->default_val("output/saxs_fitter/")->group("General options");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads)->group("General options");
    app.add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, "The unit of the q values in the measurement file. Options: A, nm.")->group("General options");
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from the measurement file.")->default_val(settings::axes::qmax)->group("General options");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from the measurement file.")->default_val(settings::axes::qmin)->group("General options");
    auto p_settings = app.add_option("-s,--settings", s_settings, "Path to the settings file.")->check(CLI::ExistingFile)->group("General options");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");

    // protein options group
    app.add_flag("--center,!--no-center", settings::molecule::center, "Decides whether the protein will be centered.")->default_val(settings::molecule::center)->group("Protein options");
    app.add_flag("--effective-charge,!--no-effective-charge", settings::molecule::use_effective_charge, "Decides whether the effective atomic charge will be used.")->default_val(settings::molecule::use_effective_charge)->group("Protein options");
    app.add_flag("--keep-hydration,!--discard-hydration", use_existing_hydration, "Decides whether the hydration layer will be generated from scratch or if the existing one will be used.")->default_val(use_existing_hydration)->group("Protein options");
    app.add_flag("--fit-excluded-volume,!--no-fit-excluded-volume", settings::hist::fit_excluded_volume, "Decides whether the excluded volume will be fitted.")->default_val(settings::hist::fit_excluded_volume)->group("Protein options");
    app.add_flag("--use-occupancy,!--ignore-occupancy", settings::molecule::use_occupancy, "Decides whether the atomic occupancies from the file will be used.")->default_val(settings::molecule::use_occupancy)->group("Protein options");

    // advanced options group
    app.add_option("--reduce,-r", settings::grid::water_scaling, "The desired number of water molecules as a percentage of the number of atoms. Use 0 for no reduction.")->default_val(settings::grid::water_scaling)->group("Advanced options");
    app.add_option("--grid_width,--gw", settings::grid::width, "The distance between each grid point in Ångström. Lower widths increase the precision.")->default_val(settings::grid::width)->group("Advanced options");
    app.add_option_function<std::string>("--hydration-strategy,--hs", [] (const std::string& s) {settings::detail::parse_option("hydration_strategy", {s});}, "The hydration model to use. Options: Radial, Axes, Jan.")->group("Advanced options");
    app.add_option("--exv_radius,--er", settings::grid::exv_radius, "The radius of the excluded volume sphere used for the grid-based excluded volume calculations in Ångström.")->default_val(settings::grid::exv_radius)->group("Advanced options");
    app.add_option_function<std::string>("--histogram-manager,--hm", [] (const std::string& s) {settings::detail::parse_option("histogram_manager", {s});}, "The histogram manager to use. Options: HM, HMMT, HMMTFF, PHM, PHMMT, PHMMTFF.")->group("Advanced options");
    app.add_flag("--exit-on-unknown-atom,!--no-exit-on-unknown-atom", settings::molecule::throw_on_unknown_atom, "Decides whether the program will exit if an unknown atom is encountered.")->default_val(settings::molecule::throw_on_unknown_atom)->group("Advanced options");
    app.add_flag("--implicit-hydrogens,!--no-implicit-hydrogens", settings::molecule::implicit_hydrogens, "Decides whether implicit hydrogens will be added to the structure.")->default_val(settings::molecule::implicit_hydrogens)->group("Advanced options");
    app.add_flag("--keep-hydrogens,!--discard-hydrogens", settings::general::keep_hydrogens, "Decides whether hydrogens will be kept in the structure.")->default_val(settings::general::keep_hydrogens)->group("Advanced options");

    // hidden options group
    app.add_option("--rvol", settings::grid::rvol, "The radius of the excluded volume sphere around each atom.")->default_val(settings::grid::rvol)->group("");
    app.add_option("--surface-thickness", settings::grid::surface_thickness, "The thickness of the surface layer in Ångström.")->default_val(settings::grid::surface_thickness)->group("");
    app.add_flag("--weighted-bins", settings::hist::weighted_bins, "Decides whether the weighted bins will be used.")->default_val(settings::hist::weighted_bins)->group("");
    app.add_flag("--save-exv", settings::grid::save_exv, "Decides whether the excluded volume will be saved.")->default_val(settings::grid::save_exv)->group("");
    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
        io::ExistingFile pdb(s_pdb), mfile(s_mfile), settings(s_settings);
        settings::general::output += mfile.stem() + "/";

        // if a settings file was provided
        if (p_settings->count() != 0) {
            settings::read(settings);       // read it
            CLI11_PARSE(app, argc, argv);   // re-parse the command line arguments so they take priority
        } else {                            // otherwise check if there is a settings file in the same directory
            if (settings::discover(mfile.directory())) {
                CLI11_PARSE(app, argc, argv);
            }
        }
        settings::validate_settings();

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
        if (settings::molecule::implicit_hydrogens) {protein.add_implicit_hydrogens();}
        if (!use_existing_hydration || protein.size_water() == 0) {
            if (protein.size_water() != 0) {console::print_text("\tDiscarding existing hydration atoms.");}
            protein.generate_new_hydration();
        }

        std::shared_ptr<fitter::HydrationFitter> fitter;
        if (settings::hist::fit_excluded_volume) {fitter = std::make_shared<fitter::ExcludedVolumeFitter>(mfile, protein.get_histogram());}
        else {fitter = std::make_shared<fitter::HydrationFitter>(mfile, protein.get_histogram());}
        auto result = fitter->fit();
        fitter::FitReporter::report(result.get());
        fitter::FitReporter::save(result.get(), settings::general::output + "report.txt", argc, argv);

        plots::PlotDistance::quick_plot(fitter->get_scattering_hist(), settings::general::output + "p(r)." + settings::plots::format);
        plots::PlotProfiles::quick_plot(fitter->get_scattering_hist(), settings::general::output + "profiles." + settings::plots::format);

        // save fit
        fitter->get_model_dataset().save(settings::general::output + "ausaxs.fit");
        fitter->get_dataset().save(settings::general::output + mfile.stem() + ".scat");

        // calculate rhoM
        double rhoM = protein.get_absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
        std::cout << "RhoM is " << rhoM << " g/cm³" << std::endl;

        protein.save(settings::general::output + "model.pdb");
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}