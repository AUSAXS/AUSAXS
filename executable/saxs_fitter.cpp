#include <CLI/CLI.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
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
#include <utility/Logging.h>

#include <vector>
#include <string>
#include <iostream>

using namespace ausaxs;

int main(int argc, char const *argv[]) {
    std::ios_base::sync_with_stdio(false);
    io::ExistingFile pdb, mfile, exv_ref_file, settings;
    bool use_existing_hydration = false, save_settings = false, save_grid = false, save_exv = false;

    CLI::App app{"Generate a new hydration layer and fit the resulting scattering intensity histogram for a given input data file."};
    app.fallthrough();
    auto input_s = app.add_option("input_structure", pdb, "Path to the structure file.")->check(CLI::ExistingFile);
    auto input_m = app.add_option("input_measurement", mfile, "Path to the measured SAXS data.")->check(CLI::ExistingFile);
    app.add_option("--output,-o", settings::general::output, "Output folder to write the results to.")->default_val("output/saxs_fitter/")->group("General options");
    app.add_flag_callback("--licence", [] () {std::cout << constants::licence << std::endl; exit(0);}, "Print the licence.");
    app.add_flag_callback("-v,--version", [] () {std::cout << constants::version << std::endl; exit(0);}, "Print the AUSAXS version.");
    app.add_flag("!--ignore-unknown-atom", settings::molecule::throw_on_unknown_atom, 
        "Do not exit upon encountering an unknown atom. This is not enabled by default to ensure awareness of potential issues.")
        ->default_val(settings::molecule::throw_on_unknown_atom);
    app.add_flag("--offline", settings::general::offline, "Run the program in offline mode. This will prevent any network requests.")
        ->default_val(settings::general::offline);
    app.add_option("--threads,-t", settings::general::threads, "Number of threads to use.")->default_val(settings::general::threads);

    // config subcommands
    auto sub_config = app.add_subcommand("config", "See and set additional options for the configuration.");
    auto p_settings = sub_config->add_option("--file,-f", settings, "The configuration file to use.")->check(CLI::ExistingFile);
    sub_config->add_flag("--save", save_settings, "Save the settings to a file.");
    sub_config->add_flag_callback("--log", [] () {logging::start("saxs_fitter");}, "Enable logging to a file.");

    // data subcommands
    auto sub_data = app.add_subcommand("data", "See and set additional options for the SAXS data.");
    sub_data->add_option(
        "--qmax", 
        settings::axes::qmax, 
        "Upper limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmax)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    sub_data->add_option(
        "--qmin", 
        settings::axes::qmin, 
        "Lower limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmin)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    sub_data->add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});}, 
        "The unit of the q values in the measurement file. Options: A, nm.");
    sub_data->add_option("--skip", settings::axes::skip, "Number of points to skip in the measurement file.")->default_val(settings::axes::skip);
    sub_data->add_flag("--rebin", settings::flags::data_rebin, "Rebin the data to increase the information content of each data point.")->default_val(settings::flags::data_rebin);

    // molecule subcommands
    auto sub_mol = app.add_subcommand("molecule", "See and set additional options for the molecular structure file.");
    sub_mol->add_flag("--center,!--no-center", settings::molecule::center, 
        "Decides whether the protein will be centered.")->default_val(settings::molecule::center);
    sub_mol->add_flag("--use-occupancy,!--ignore-occupancy", settings::molecule::use_occupancy, 
        "Decides whether the atomic occupancies from the file will be used.")->default_val(settings::molecule::use_occupancy);

    // exv subcommands
    auto sub_exv = app.add_subcommand("exv", "See and set additional options for the excluded volume calculations.");
    sub_exv->add_option_function<std::string>("--model,-m", [] (const std::string& s) 
        {settings::detail::parse_option("exv_model", {s});}, 
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
    sub_exv->add_flag("--save", save_exv, 
        "Write a PDB representation of the excluded volume to disk."
    )->default_val(save_exv);
    sub_exv->add_option("--ref,--reference", exv_ref_file, 
        "Path to the excluded volume reference file."
    )->check(CLI::ExistingFile);

    // solvation subcommands
    auto sub_water = app.add_subcommand("solv", "See and set additional options for the solvation calculations.");
    sub_water->add_option_function<std::string>("--model,-m", [] (const std::string& s) {settings::detail::parse_option("hydration_strategy", {s});}, 
        "The hydration model to use. Options: Radial, Axes, None.");
    sub_water->add_flag("--keep,!--discard", use_existing_hydration, 
        "Keep or discard water molecules from the structure file. "
        "If they are discarded, a new solvation shell is generated."
    )->default_val(use_existing_hydration);
    sub_water->add_flag("--fit,!--no-fit", settings::fit::fit_hydration, 
        "Fit the hydration shell.")->default_val(settings::fit::fit_hydration);

    // hydrogen subcommands
    auto sub_hydrogen = app.add_subcommand("hydrogens", "See and set additional options for the handling of hydration atoms.");
    sub_hydrogen->add_flag("--keep,!--discard", settings::general::keep_hydrogens, "Keep or discard hydrogens from the structure file.")->default_val(settings::general::keep_hydrogens);

    // grid subcommands
    auto sub_grid = app.add_subcommand("grid", "See and set additional options for the grid calculations.");
    sub_grid->add_option("--rvol", settings::grid::min_exv_radius, 
        "The radius of the excluded volume sphere around each atom."
    )->default_val(settings::grid::min_exv_radius);
    auto sub_grid_w = sub_grid->add_option("--width,-w", settings::grid::cell_width, 
        "The distance between each grid point in Ångström. Lower widths increase the precision."
    )->default_val(settings::grid::cell_width);
    sub_grid->add_flag("--save", save_grid, 
        "Write a PDB representation of the grid to disk."
    )->default_val(save_grid);

    // fit subcommands
    // auto sub_fit = app.add_subcommand("fit", "See and set additional options for the fitting process.");
    // sub_fit->add_flag("--atomic-debye-waller", settings::fit::fit_atomic_debye_waller, 
    //     "Fit the atomic form factor debye-waller factor."
    // )->default_val(settings::fit::fit_atomic_debye_waller);
    // sub_fit->add_flag("--exv-debye-waller", settings::fit::fit_exv_debye_waller, 
    //     "Fit the excluded volume form factor debye-waller factor."
    // )->default_val(settings::fit::fit_exv_debye_waller);

    // hidden options group
    app.add_flag("--weighted-bins", settings::hist::weighted_bins, 
        "Decides whether the weighted bins will be used."
    )->default_val(settings::hist::weighted_bins)->group("");

    app.final_callback([&] () {
        // save settings if requested
        if (save_settings) {
            settings::write("settings.txt");
            console::print_info("Settings saved to settings.txt in current directory.");
            if (!input_s->count() || !input_m->count()) { // gracefully exit if no input files are provided
                exit(0);
            }
        }

        // required args (not marked ->required() since that interferes with the help flag for subcommands)
        if (!input_s->count()) {
            std::cout << "Error: input_structure is required." << std::endl;
            exit(1);
        }

        // adjust grid width to support user-specified excluded volume width
        if (sub_exv_w->count() && !sub_grid_w->count()) {
            settings::grid::cell_width = settings::grid::exv::width;
        }

        // adjust excluded volume width to be at least as large as the grid width
        if (sub_grid_w->count() && !sub_exv_w->count() && settings::grid::exv::width < settings::grid::cell_width) {
            settings::grid::exv::width = settings::grid::cell_width;
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
                throw except::invalid_argument("Unknown PDB extensions: " + pdb.str() + " and " + mfile.str());
            }
        }

        if (!mfile.empty() && !constants::filetypes::saxs_data.check(mfile)) {
            throw except::invalid_argument("Unknown SAXS data extension: " + mfile.str());
        }

        //######################//
        //### ACTUAL PROGRAM ###//
        //######################//

        data::Molecule protein(pdb);
        if (!exv_ref_file.empty()) {protein.set_grid(grid::Grid::create_from_reference(exv_ref_file, protein));}
        if (!use_existing_hydration || protein.size_water() == 0) {
            if (protein.size_water() != 0) {console::print_text("\tDiscarding existing hydration atoms.");}
            protein.generate_new_hydration();
        }
        std::string msg_exv_vol, msg_solv_dens;

        // simulation mode
        if (mfile.empty()) {
            console::print_info("\nSimulation mode enabled.");
            console::print_text("Please note that the evaluated hydration shell contribution will be quite poor for most molecules in this mode. For more information, refer to the documentation.");
            settings::general::output += "simulated/" + pdb.stem() + "/";
            auto hist = protein.get_histogram();

            plots::PlotDistance::quick_plot(hist.get(), settings::general::output + "p(r)." + settings::plots::format);
            plots::PlotProfiles::quick_plot(hist.get(), settings::general::output + "profiles." + settings::plots::format);
        }

        // fitting mode
        else {
            settings::general::output += mfile.stem() + "/";
            SimpleDataset saxs_data(mfile);
            if (settings::flags::data_rebin) {console::indent(); saxs_data.rebin(); console::unindent();}
    
            fitter::SmartFitter fitter(std::move(saxs_data), protein.get_histogram());
            auto result = fitter.fit();
            fitter::FitReporter::report(result.get());
            fitter::FitReporter::save(result.get(), settings::general::output + "report.txt", argc, argv);
                
            plots::PlotDistance::quick_plot(fitter.get_model(), settings::general::output + "p(r)." + settings::plots::format);
            plots::PlotProfiles::quick_plot(fitter.get_model(), settings::general::output + "profiles." + settings::plots::format);
            result->curves.select_columns({0, 1, 2, 3}).save(
                settings::general::output + "ausaxs.fit", 
                "chi2=" + std::to_string(result->fval/result->dof) + " dof=" + std::to_string(result->dof)
            );
            if (settings::fit::fit_excluded_volume) {
                msg_exv_vol = 
                    "\tExcluded:        " 
                    + std::to_string((int) std::round(protein.get_volume_exv(settings::fit::fit_excluded_volume ? result->get_parameter(constants::fit::Parameters::SCALING_EXV) : 1))) 
                    + " A^3"
                ;
            }
            if (settings::fit::fit_solvent_density) {
                msg_solv_dens = 
                    "\tSolvent density: " 
                    + std::to_string(constants::charge::density::water*result->get_parameter(constants::fit::Parameters::SCALING_RHO)) 
                    + " e/A^3"
                ;
            }
        }

        // calculate extra stuff
        console::print_info("\nExtra informaton");
        console::print_text("Volume:");
        double rhoM = protein.get_absolute_mass()/protein.get_volume_grid()*constants::unit::gm/(std::pow(constants::unit::cm, 3));
        double exv_vol = protein.get_volume_grid();
        double mol_charge = protein.get_total_atomic_charge();
        console::print_text("\tvan der Waals:   " + std::to_string((int) std::round(protein.get_volume_vdw()))  + " A^3");
        console::print_text("\tGrid:            " + std::to_string((int) std::round(protein.get_volume_grid())) + " A^3");
        if (settings::fit::fit_excluded_volume) {console::print_text(msg_exv_vol);}
        console::print_text("\nCharge:");
        console::print_text("\tMolecular:       " + utility::round_double(mol_charge, 1) + " e");
        console::print_text("\tExcluded volume: " + utility::round_double(exv_vol*constants::charge::density::water, 1) + " e");
        console::print_text("\tExcess density:  " + utility::round_double((mol_charge - exv_vol*constants::charge::density::water)/exv_vol, 3) + " e/A^3");
        if (settings::fit::fit_solvent_density) {console::print_text(msg_solv_dens);}
        console::print_text("\nOther properties:");
        console::print_text("\tRhoM:            " + utility::round_double(rhoM, 3) + " g/cm^3");

        protein.save(settings::general::output + "model.pdb");
        if (save_grid) {protein.get_grid()->save(settings::general::output + "grid.pdb");}
        if (save_exv) {protein.get_grid()->generate_excluded_volume().save(settings::general::output + "exv.pdb");}
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        throw e;
    }
    return 0;
}