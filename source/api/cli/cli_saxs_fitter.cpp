// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_saxs_fitter.h>

#include <api/cli/cli_sans_fitter.h>
#include <api/cli/detail/CLICommon.h>
#include <constants/ConstantsFitParameters.h>
#include <constants/ConstantsUnits.h>
#include <constants/ValidFileExtensions.h>
#include <data/Molecule.h>
#include <fitter/FitReporter.h>
#include <fitter/SmartFitter.h>
#include <grid/Grid.h>
#include <hist/scattering_model/IScatteringModel.h>  // IWYU pragma: keep
#include <io/ExistingFile.h>
#include <plots/PlotDistance.h>
#include <plots/PlotProfiles.h>
#include <settings/All.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/StringUtils.h>

#include <CLI/CLI.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace ausaxs;

namespace {
    using Radiation = settings::scattering::Radiation;

    /**
     * @brief Fit a structure to a measurement, or simulate its profiles if no measurement is given.
     *        The probe is decided by the tool, after any settings file has been read.
     */
    int run_fitter(int argc, char const* const* argv, Radiation radiation) {
        const bool neutron = radiation == Radiation::Neutron;
        const std::string probe_name = neutron ? "SANS" : "SAXS";
        const std::string tool_name = neutron ? "sans_fitter" : "saxs_fitter";

        io::ExistingFile pdb, mfile, exv_ref_file;
        cli::ConfigOptions config;
        bool keep_hydration = false, save_grid = false, save_exv = false;

        CLI::App app{"Perform " + probe_name + " fitting for a given structure and measurement."};
        app.fallthrough();
        auto* input_s = app.add_option("input_structure", pdb, "Path to the structure file.")->check(CLI::ExistingFile);
        auto* input_m = app.add_option("input_measurement", mfile, "Path to the measured " + probe_name + " data.")->check(CLI::ExistingFile);
        cli::add_general_options(app, "output/" + tool_name + "/");
        cli::add_structure_options(app);
        cli::add_config_subcommand(app, config, tool_name);
        cli::add_data_subcommand(app, probe_name);
        cli::add_molecule_subcommand(app);

        // exv subcommands
        auto* sub_exv = app.add_subcommand("exv", "See and set additional options for the excluded volume calculations.");
        sub_exv->add_option_function<std::string>("--model,-m", [] (const std::string& s)
            {settings::detail::parse_option("exv_model", {s});},
            "The excluded volume model to use. Options: Simple, Fraser, Grid.");
        sub_exv->add_flag("--fit", settings::fit::fit_excluded_volume.value,
            "Fit the excluded volume.")->default_val(settings::fit::fit_excluded_volume.value);
        sub_exv->add_flag("--fit-density", settings::fit::fit_solvent_density.value,
            "Fit the solvent density for the excluded volume.")->default_val(settings::fit::fit_solvent_density.value);
        sub_exv->add_option_function<std::string>("--table", [] (const std::string& s) {settings::detail::parse_option("exv_volume", {s});},
            "The set of displaced volume tables to use. Options: Traube, vdw, Voronoi, Minimum fluctuation (mf)."
        );
        sub_exv->add_option("--surface-thickness", settings::grid::exv::surface_thickness.value,
            "The thickness of the surface layer in Ångström."
        )->default_val(settings::grid::exv::surface_thickness.value)->group("");
        auto* sub_exv_w = sub_exv->add_option("--width,-w", settings::grid::exv::width.value,
            "The width of the excluded volume dummy atoms used for the grid-based excluded volume calculations in Ångström."
        )->default_val(settings::grid::exv::width.value);
        sub_exv->add_flag("--save", save_exv,
            "Write a PDB representation of the excluded volume to disk."
        )->default_val(save_exv);
        sub_exv->add_option("--ref,--reference", exv_ref_file,
            "Path to the excluded volume reference file."
        )->check(CLI::ExistingFile);

        // solvation subcommands
        auto* sub_water = app.add_subcommand("solv", "See and set additional options for the solvation calculations.");
        sub_water->add_option_function<std::string>("--model,-m", [] (const std::string& s) {settings::detail::parse_option("hydration_strategy", {s});},
            "The hydration model to use. Options: Radial, Axes, None.");
        sub_water->add_flag("--keep,!--discard", keep_hydration,
            "Keep or discard water molecules from the structure file. "
            "If they are discarded, a new solvation shell is generated."
        )->default_val(keep_hydration);
        sub_water->add_flag("--fit,!--no-fit", settings::fit::fit_hydration.value,
            "Fit the hydration shell.")->default_val(settings::fit::fit_hydration.value);
        if (neutron) {cli::add_neutron_solvent_options(sub_water);}

        cli::add_hydrogens_subcommand(app);

        // grid subcommands
        auto* sub_grid = app.add_subcommand("grid", "See and set additional options for the grid calculations.");
        sub_grid->add_option("--rvol", settings::grid::min_exv_radius.value,
            "The radius of the excluded volume sphere around each atom."
        )->default_val(settings::grid::min_exv_radius.value);
        auto* sub_grid_w = sub_grid->add_option("--width,-w", settings::grid::cell_width.value,
            "The distance between each grid point in Ångström. Lower widths increase the precision."
        )->default_val(settings::grid::cell_width.value);
        sub_grid->add_flag("--save", save_grid,
            "Write a PDB representation of the grid to disk."
        )->default_val(save_grid);

        app.final_callback([&] () {
            // save settings if requested
            if (config.save) {
                settings::write("settings.txt");
                console::print_info("Settings saved to settings.txt in current directory.");
                if (!input_s->count() || !input_m->count()) { // gracefully exit if no input files are provided
                    throw CLI::Success();
                }
            }

            // required args (not marked ->required() since that interferes with the help flag for subcommands)
            if (!input_s->count()) {
                console::print_warning("Error: input_structure is required.");
                throw CLI::RuntimeError(1);
            }

            // adjust grid width to support user-specified excluded volume width
            if (sub_exv_w->count() && !sub_grid_w->count()) {
                settings::grid::cell_width = static_cast<double>(settings::grid::exv::width);
            }

            // adjust excluded volume width to be at least as large as the grid width
            if (sub_grid_w->count() && !sub_exv_w->count() && settings::grid::exv::width < settings::grid::cell_width) {
                settings::grid::exv::width = static_cast<double>(settings::grid::cell_width);
            }
        });

        //###################//
        //### PARSE INPUT ###//
        //###################//
        try {
            if (auto code = cli::parse(app, argc, argv, config, [&] () {return mfile.directory();})) {return *code;}
            settings::scattering::radiation = radiation; // the tool decides the probe, also over a settings file
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
                throw except::invalid_argument("Unknown " + probe_name + " data extension: " + mfile.str());
            }

            //######################//
            //### ACTUAL PROGRAM ###//
            //######################//

            if (neutron) {cli::print_neutron_solvent();}
            data::Molecule molecule(pdb);
            if (!exv_ref_file.empty()) {molecule.set_grid(grid::Grid::create_from_reference(exv_ref_file, molecule));}
            if (!keep_hydration) {
                if (molecule.size_water() != 0) {console::print_text("\tDiscarding existing hydration atoms.");}
                molecule.generate_new_hydration();
            } else if (molecule.size_water() == 0) {
                console::print_warning("No hydration atoms were found in the structure file, but --keep was specified. No hydration atoms will be used for this calculation.");
            }
            std::string msg_exv_vol, msg_solv_dens;
            const std::string charge_unit = neutron ? "fm" : "e";

            // simulation mode
            if (mfile.empty()) {
                console::print_info("\nSimulation mode enabled.");
                console::print_text("Please note that the evaluated hydration shell contribution will be quite poor for most molecules in this mode. For more information, refer to the documentation.");
                settings::general::output += "simulated/" + pdb.stem() + "/";
                auto hist = molecule.get_histogram();

                plots::PlotDistance::quick_plot(hist.get(), settings::general::output + "p(r)." + settings::plots::format);
                plots::PlotProfiles::quick_plot(hist.get(), settings::general::output + "profiles." + settings::plots::format);
            }

            // fitting mode
            else {
                settings::general::output += mfile.stem() + "/";
                SimpleDataset data(mfile);
                if (settings::axes::rebin) {
                    console::indent();
                    data.rebin();
                    console::unindent();
                    data.save(settings::general::output + "rebinned.dat");
                }

                fitter::SmartFitter fitter(data, molecule.get_histogram());
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
                        + std::to_string((int) std::round(molecule.get_volume_exv(result->get_parameter(constants::fit::Parameters::SCALING_EXV))))
                        + " A^3"
                    ;
                }
                if (settings::fit::fit_solvent_density) {
                    msg_solv_dens =
                        "\tSolvent density: "
                        + std::to_string(settings::internal_state::solvent_density()*result->get_parameter(constants::fit::Parameters::SCALING_RHO))
                        + " " + charge_unit + "/A^3"
                    ;
                }
            }

            // calculate extra stuff
            console::print_info("\nExtra informaton");
            console::print_text("Volume:");
            // the grid volume is dry, so the mass must be dry too
            double rhoM = molecule.get_absolute_mass(false)/molecule.get_volume_grid()*constants::unit::external::gm/constants::unit::external::mL;
            double exv_vol = molecule.get_volume_grid();
            double mol_charge = molecule.get_total_atomic_charge();
            console::print_text("\tvan der Waals:   " + std::to_string((int) std::round(molecule.get_volume_vdw()))  + " A^3");
            console::print_text("\tGrid:            " + std::to_string((int) std::round(molecule.get_volume_grid())) + " A^3");
            if (settings::fit::fit_excluded_volume) {console::print_text(msg_exv_vol);}
            console::print_text("\nCharge:");
            console::print_text("\tMolecular:       " + utility::round_double(mol_charge, 1) + " " + charge_unit);
            console::print_text("\tExcluded volume: " + utility::round_double(exv_vol*settings::internal_state::solvent_density(), 1) + " " + charge_unit);
            console::print_text("\tExcess density:  " + utility::round_double((mol_charge - exv_vol*settings::internal_state::solvent_density())/exv_vol, 3) + " " + charge_unit + "/A^3");
            if (settings::fit::fit_solvent_density) {console::print_text(msg_solv_dens);}
            console::print_text("\nOther properties:");
            console::print_text("\tRhoM:            " + utility::round_double(rhoM, 3) + " g/cm^3");

            molecule.save(settings::general::output + "model.pdb");
            if (save_grid) {molecule.get_grid()->save(settings::general::output + "grid.pdb");}
            if (save_exv) {molecule.get_grid()->generate_excluded_volume().save(settings::general::output + "exv.pdb");}
        } catch (const except::base&) {
            return 1; // our own exceptions already report themselves when thrown
        } catch (const std::exception& e) {
            console::print_warning(e.what());
            return 1;
        }
        return 0;
    }
}

int cli_saxs_fitter(int argc, char const *argv[]) {
    return run_fitter(argc, argv, Radiation::XRay);
}

int cli_sans_fitter(int argc, char const *argv[]) {
    return run_fitter(argc, argv, Radiation::Neutron);
}
