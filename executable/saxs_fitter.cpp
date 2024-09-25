#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/record/Water.h>
#include <data/Molecule.h>
#include <fitter/SmartFitter.h>
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

using namespace ausaxs;

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    io::ExistingFile pdb, mfile, settings;
    settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
    bool use_existing_hydration = false, save_settings = false;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.fallthrough();
    auto input_s = app.add_option("input_structure", pdb, "Path to the structure file.")->check(CLI::ExistingFile);
    auto input_m = app.add_option("input_measurement", mfile, "Path to the measured SAXS data.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Output folder to write the results to.")->default_val("output/saxs_fitter/")->group("General options");
    auto p_settings = app.add_option("-s,--settings", settings, "Path to the settings file.")->check(CLI::ExistingFile)->group("General options");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    app.add_flag_callback("-v,--version", [] () {std::cout << constants::version << std::endl; exit(0);}, "Print the AUSAXS version.");

    // dataset options group
    app.add_option("--qmax", settings::axes::qmax, "Upper limit on used q values from the measurement file.")->default_val(settings::axes::qmax)->group("Dataset options");
    app.add_option("--qmin", settings::axes::qmin, "Lower limit on used q values from the measurement file.")->default_val(settings::axes::qmin)->group("Dataset options");
    app.add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, 
        "The unit of the q values in the measurement file. Options: A, nm.")->group("Dataset options");
    app.add_option("--skip", settings::axes::skip, "Number of points to skip in the measurement file.")->default_val(settings::axes::skip)->group("Dataset options");

    // advanced options
    app.add_flag("!--ignore-unknown-atom", settings::molecule::throw_on_unknown_atom, 
        "Do not exit upon encountering an unknown atom. This is not enabled by default to ensure you are aware of the issue."
    )->default_val(settings::molecule::throw_on_unknown_atom)->group("Advanced options");
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads)->group("Advanced options");
    app.add_flag("--save-settings", save_settings, "Save the settings to a file.")->default_val(save_settings)->group("Advanced options");

    // molecule subcommands
    auto sub_mol = app.add_subcommand("molecule", "See and set additional options for the molecular structure file.");
    sub_mol->add_flag("--center,!--no-center", settings::molecule::center, 
        "Decides whether the protein will be centered.")->default_val(settings::molecule::center);
    sub_mol->add_flag("--effective-charge,!--no-effective-charge", settings::molecule::use_effective_charge, 
        "Decides whether the effective atomic charge will be used.")->default_val(settings::molecule::use_effective_charge);
    sub_mol->add_flag("--use-occupancy,!--ignore-occupancy", settings::molecule::use_occupancy, 
        "Decides whether the atomic occupancies from the file will be used.")->default_val(settings::molecule::use_occupancy);

    // exv subcommands
    auto sub_exv = app.add_subcommand("exv", "See and set additional options for the excluded volume calculations.");
    sub_exv->add_option_function<std::string>("--model,-m", [] (const std::string& s) 
        {
            settings::detail::parse_option("histogram_manager", {s});
            if (s == "fraser" || s == "grid") {settings::molecule::use_effective_charge = false;} // ensure correct setup for standard args
        }, 
        "The excluded volume model to use. Options: Simple, Fraser, Grid.");
    sub_exv->add_flag("--fit", settings::fit::fit_excluded_volume, 
        "Fit the excluded volume.")->default_val(settings::fit::fit_excluded_volume);
    sub_exv->add_flag("--fit-density", settings::fit::fit_solvent_density, 
        "Fit the solvent density for the excluded volume.")->default_val(settings::fit::fit_solvent_density);
    sub_exv->add_option_function<std::string>("--table", [] (const std::string& s) {settings::detail::parse_option("exv_volume", {s});}, 
        "The set of displaced volume tables to use. Options: Traube, vdw, Voronoi, Minimum fluctuation (mf)."
    );
    sub_exv->add_option("--surface-thickness", settings::grid::exv::surface_thickness, 
        "The thickness of the surface layer in Ångström."
    )->default_val(settings::grid::exv::surface_thickness)->group("");

    auto sub_exv_w = sub_exv->add_option("--width,-w", settings::grid::exv::width, 
        "The width of the excluded volume dummy atoms used for the grid-based excluded volume calculations in Ångström."
    )->default_val(settings::grid::exv::width);
    sub_exv->add_flag("--save", settings::grid::exv::save, 
        "Write a PDB representation of the excluded volume to disk."
    )->default_val(settings::grid::exv::save);

    // solvation subcommands
    auto sub_water = app.add_subcommand("solv", "See and set additional options for the solvation calculations.");
    sub_water->add_option_function<std::string>("--model,-m", [] (const std::string& s) {settings::detail::parse_option("hydration_strategy", {s});}, 
        "The hydration model to use. Options: Radial, Axes, None.");
    sub_water->add_flag("--keep,!--discard", use_existing_hydration, 
        "Keep or discard water molecules from the structure file. "
        "If they are discarded, a new solvation shell is generated."
    )->default_val(use_existing_hydration);
    sub_water->add_flag("--fit", settings::fit::fit_hydration, 
        "Fit the hydration shell.")->default_val(settings::fit::fit_hydration);
    sub_water->add_option("--reduce,-r", settings::grid::water_scaling, 
        "Reduce the number of generated water molecules to a percentage of the number of atoms. "
        "Use 0 for no reduction."
    )->default_val(settings::grid::water_scaling);

    // hydrogen subcommands
    auto sub_hydrogen = app.add_subcommand("hydrogens", "See and set additional options for the handling of hydration atoms.");
    sub_hydrogen->add_flag("--keep,!--discard", settings::general::keep_hydrogens, 
        "Keep or discard hydrogens from the structure file.")->default_val(settings::general::keep_hydrogens);
    sub_hydrogen->add_flag("--implicit,!--explicit", settings::molecule::implicit_hydrogens, 
        "Add implicit hydrogens to the structure. "
        "This should only be disabled if they are explicitly provided in the structure file. "
        "This option also switches between the implicit and explicit hydrogen variants of the excluded volume tables."
    )->default_val(settings::molecule::implicit_hydrogens)->group("Model options");

    // grid subcommands
    auto sub_grid = app.add_subcommand("grid", "See and set additional options for the grid calculations.");
    sub_grid->add_option("--rvol", settings::grid::min_exv_radius, 
        "The radius of the excluded volume sphere around each atom."
    )->default_val(settings::grid::min_exv_radius)->group("");
    auto sub_grid_w = sub_grid->add_option("--width,-w", settings::grid::cell_width, 
        "The distance between each grid point in Ångström. Lower widths increase the precision."
    )->default_val(settings::grid::cell_width);

    // hidden options group
    app.add_flag("--weighted-bins", settings::hist::weighted_bins, 
        "Decides whether the weighted bins will be used."
    )->default_val(settings::hist::weighted_bins)->group("");

    app.final_callback([&] () {
        // required args (not marked ->required() since that interferes with the help flag for subcommands)
        if (!input_s->count() || !input_m->count()) {
            std::cout << "Error: Both input_structure and input_measurement are required." << std::endl;
            exit(1);
        }

        // adjust grid width to support user-specified excluded volume width
        if (sub_exv_w->count() && !sub_grid_w->count()) {
            settings::grid::cell_width = settings::grid::exv::width;
        }

        // save settings if requested
        if (save_settings) {
            settings::write("settings.txt");
            console::print_info("Settings saved to settings.txt in current directory.");
        }
    });

    CLI11_PARSE(app, argc, argv);

    console::print_info("Running AUSAXS " + std::string(constants::version));

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
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
        if (!constants::filetypes::structure.check(pdb)) {
            // check if the two inputs are switched
            if (constants::filetypes::structure.check(mfile)) {
                // if so, silently swap them and proceed
                std::swap(pdb, mfile);
            } else {
                throw except::invalid_argument("Unknown PDB extensions: " + pdb + " and " + mfile);
            }
        }
        if (!constants::filetypes::saxs_data.check(mfile)) {
            throw except::invalid_argument("Unknown SAXS data extension: " + mfile);
        }

        settings::general::output += mfile.stem() + "/";

        //######################//
        //### ACTUAL PROGRAM ###//
        //######################//
        data::Molecule protein(pdb);
        if (settings::molecule::implicit_hydrogens) {protein.add_implicit_hydrogens();}
        if (!use_existing_hydration || protein.size_water() == 0) {
            if (protein.size_water() != 0) {console::print_text("\tDiscarding existing hydration atoms.");}
            protein.generate_new_hydration();
        }

        fitter::SmartFitter fitter(mfile, protein.get_histogram());
        auto result = fitter.fit();
        fitter::FitReporter::report(result.get());
        fitter::FitReporter::save(result.get(), settings::general::output + "report.txt", argc, argv);

        plots::PlotDistance::quick_plot(fitter.get_model(), settings::general::output + "p(r)." + settings::plots::format);
        plots::PlotProfiles::quick_plot(fitter.get_model(), settings::general::output + "profiles." + settings::plots::format);
        result->curves.save(
            settings::general::output + "ausaxs.fit", 
            "chi2=" + std::to_string(result->fval/result->dof) + " dof=" + std::to_string(result->dof)
        );

        // calculate extra stuff
        console::print_info("\nExtra informaton");
        console::indent();
        double rhoM = protein.get_absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
        double d = settings::fit::fit_excluded_volume ? result->get_parameter("d") : 1;
        console::print_text("Volume (vdW):  " + std::to_string((int) std::round(protein.get_volume_vdw()))  + " Å^3");
        console::print_text("Volume (grid): " + std::to_string((int) std::round(protein.get_volume_grid())) + " Å^3");
        console::print_text("Volume (exv):  " + std::to_string((int) std::round(protein.get_volume_exv(d)))  + " Å^3");
        console::print_text("RhoM:          " + std::to_string(rhoM) + " g/cm^3");
        console::unindent();

        protein.save(settings::general::output + "model.pdb");
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}