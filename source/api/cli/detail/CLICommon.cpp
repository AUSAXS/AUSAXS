// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/detail/CLICommon.h>

#include <constants/ConstantsAxes.h>
#include <constants/Version.h>
#include <form_factor/neutron/NeutronSolvent.h>
#include <gpu/GPUInstaller.h>
#include <settings/All.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <utility/StringUtils.h>

using namespace ausaxs;

void cli::add_general_options(CLI::App& app, const std::string& default_output) {
    app.add_option("--output,-o", settings::general::output.value, "Output folder to write the results to.")->default_val(default_output)->group("General options");
    app.add_flag_callback("--licence",    [] () {console::print_text(constants::licence); throw CLI::Success();}, "Print the licence.");
    app.add_flag_callback("-v,--version", [] () {console::print_text(constants::version); throw CLI::Success();}, "Print the AUSAXS version.");
    app.add_option("--threads,-t", settings::general::threads.value, "Number of threads to use.")->default_val(settings::general::threads.value);
    settings::general::gpu = gpu::installer::enabled(); // choosing a GPU with "ausaxs gpu install" is the opt-in
    app.add_flag("--gpu,!--no-gpu", settings::general::gpu.value,
        R"(Use GPU acceleration if available. On by default once a GPU is chosen with "ausaxs gpu install", unless turned off with "ausaxs gpu disable".)")
        ->default_val(settings::general::gpu.value);
}

void cli::add_structure_options(CLI::App& app) {
    app.add_flag("--allow-unknown-atoms", settings::molecule::allow_unknown_atoms.value,
        "Allow processing files with unknown atoms. Use only if you understand the implications.")
        ->default_val(settings::molecule::allow_unknown_atoms.value);
    app.add_flag("--allow-unknown-residues", settings::molecule::allow_unknown_residues.value,
        "Allow processing files with unknown residues. Use only if you understand the implications.")
        ->default_val(settings::molecule::allow_unknown_residues.value);
    app.add_flag("--offline", settings::general::offline.value, "Run the program in offline mode. This will prevent any network requests.")
        ->default_val(settings::general::offline.value);
}

void cli::add_config_subcommand(CLI::App& app, ConfigOptions& config, const std::string& log_name) {
    auto* sub_config = app.add_subcommand("config", "See and set additional options for the configuration.");
    config.file_option = sub_config->add_option("--file,-f", config.file, "The configuration file to use.")->check(CLI::ExistingFile);
    sub_config->add_flag("--save", config.save, "Save the settings to a file.");
    sub_config->add_flag_callback("--log", [log_name] () {logging::start(log_name);}, "Enable logging to a file.");
}

CLI::App* cli::add_data_subcommand(CLI::App& app, const std::string& data_name) {
    auto* sub_data = app.add_subcommand("data", "See and set additional options for the " + data_name + " data.");
    sub_data->add_option(
        "--qmax",
        settings::axes::qmax.value,
        "Upper limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmax.value)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    sub_data->add_option(
        "--qmin",
        settings::axes::qmin.value,
        "Lower limit on used q values from the measurement file.")
        ->default_val(settings::axes::qmin.value)
        ->check(CLI::Range(constants::axes::q_axis.min, constants::axes::q_axis.max))
    ;
    add_measurement_options(sub_data);
    sub_data->add_flag("--rebin", settings::axes::rebin.value, "Rebin the data to increase the information content of each data point.")->default_val(settings::axes::rebin.value);
    sub_data->add_flag("--weighted-bins", settings::hist::weighted_bins.value, "Use weighted bins (true/false). By default, this is automatically decided for each molecule.")->default_val(settings::hist::weighted_bins.value);
    return sub_data;
}

void cli::add_measurement_options(CLI::App* sub) {
    sub->add_option_function<std::string>("--unit,-u", [] (const std::string& s) {settings::detail::parse_option("unit", {s});},
        "The unit of the q values in the measurement file. Options: A, nm.");
    sub->add_option("--skip", settings::axes::skip.value, "Number of points to skip in the measurement file.")->default_val(settings::axes::skip.value);
}

CLI::App* cli::add_hydrogens_subcommand(CLI::App& app) {
    auto* sub_hydrogen = app.add_subcommand("hydrogens", "See and set additional options for the handling of hydrogen atoms.");
    sub_hydrogen->add_flag("--keep,!--discard", settings::general::keep_hydrogens.value, "Keep or discard hydrogens from the structure file.")->default_val(settings::general::keep_hydrogens.value);
    sub_hydrogen->add_flag("--force-field-names", settings::molecule::force_field_residue_names.value,
        "Parse forcefield residue names (protonation states HID, HSE, CYX, ...; water SOL) by their forcefield meaning instead of their usual CCD definitions.")
        ->default_val(settings::molecule::force_field_residue_names.value);
    return sub_hydrogen;
}

void cli::add_neutron_solvent_options(CLI::App* sub) {
    // assigned through the settings, since their change callbacks rebuild the form factor tables
    sub->add_option_function<double>("--d2o", [] (double x) {settings::scattering::neutron::set_d2o_fraction(x);},
        "The D2O fraction of the H2O/D2O solvent, from 0 to 1. The default is 1, pure D2O.")
        ->check(CLI::Range(0., 1.));
    sub->add_option_function<double>("--exchange", [] (double x) {settings::scattering::neutron::exchange_efficiency = x;},
        "The fraction of the exchangeable implicit hydrogens (bound to N, O or S) that takes the isotope composition of the solvent, from 0 to 1. The default is 1, full exchange.")
        ->check(CLI::Range(0., 1.));
}

std::optional<int> cli::parse(
    CLI::App& app, int argc, char const* const* argv, const ConfigOptions& config, const std::function<io::Folder()>& discover_dir
) {
    auto parse_once = [&] () -> std::optional<int> {
        try {
            app.parse(argc, argv);
        } catch (const CLI::ParseError& e) {
            return app.exit(e);
        }
        return std::nullopt;
    };

    if (auto code = parse_once()) {return code;}
    console::print_info("Running AUSAXS " + std::string(constants::version));

    if (config.file_option->count() != 0) { // a settings file was given
        settings::read(config.file);
        return parse_once();                // the command line takes priority over the file
    }
    if (discover_dir && settings::discover(discover_dir())) {
        return parse_once();
    }
    return std::nullopt;
}

void cli::add_grid_width_option(CLI::App* sub) {
    sub->add_option_function<double>("--width,-w", [] (double w) {
            settings::grid::cell_width = w;
            settings::grid::exv::width = w;
        },
        "The distance between each grid point in Ångström, used for both collision detection and grid-based excluded volume calculations."
        "Lower widths increase the precision."
    )->default_str(std::to_string(settings::grid::cell_width.value))->check(CLI::PositiveNumber);
}

void cli::print_neutron_solvent() {
    const double density = settings::scattering::neutron::solvent_density;
    console::print_text(
        "Neutron solvent: " + utility::round_double(100*form_factor::neutron::d2o_fraction(density), 1) + "% D2O "
        "(" + utility::round_double(density, 4) + " fm/A^3), exchange efficiency "
        + utility::round_double(100*settings::scattering::neutron::exchange_efficiency, 1) + "%."
    );
}
