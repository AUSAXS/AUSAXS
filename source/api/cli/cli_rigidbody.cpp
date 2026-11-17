// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_rigidbody.h>
#include <api/cli/detail/CLICommon.h>

#include <CLI/CLI.hpp>

#include <constants/ValidFileExtensions.h>
#include <fitter/FitReporter.h>
#include <io/File.h>
#include <rigidbody/sequencer/detail/SequenceParser.h>
#include <settings/All.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>
#include <utility/Logging.h>

#include <string>

using namespace ausaxs;

int cli_rigidbody(int argc, char const *argv[]) {
    settings::grid::scaling = 2;
    settings::grid::cubic = true;
    settings::general::verbose = true;
    io::File script;
    cli::ConfigOptions config;
    CLI::App app{"Perform rigid-body optimization."};
    app.fallthrough();
    auto* p_script = app.add_option("input_script", script, "Path to the rigid-body configuration script.")->check(CLI::ExistingFile);
    cli::add_general_options(app, "output/rigidbody/");
    cli::add_structure_options(app);
    cli::add_config_subcommand(app, config, "rigidbody");
    cli::add_data_subcommand(app, "SAXS");
    cli::add_molecule_subcommand(app);
    cli::add_hydrogens_subcommand(app);

    // grid subcommands
    auto* sub_grid = app.add_subcommand("grid", "See and set additional options for the grid calculations.");
    sub_grid->add_option("--width,-w", settings::grid::cell_width.value, 
        "The distance between each grid point in Ångström. Lower widths increase the precision."
    )->default_val(settings::grid::cell_width.value);

    app.final_callback([&] () {if (config.save) {
        settings::write("settings.txt");
        console::print_info("Settings saved to settings.txt in current directory.");
    }});

    //###################//
    //### PARSE INPUT ###//
    //###################//
    try {
        if (auto code = cli::parse(app, argc, argv, config)) {return *code;} // no settings file is looked for next to the script
        logging::start("rigidbody");

        // rigid-body optimization is only available through the sequencer configuration format
        if (p_script->count() == 0) {
            throw except::missing_option("rigidbody: A configuration script must be supplied.");
        }
        if (!constants::filetypes::config.check(script)) {
            throw except::invalid_argument(
                "rigidbody: The input file \"" + script.str() + "\" is not a configuration script. "
                "Rigid-body optimization must be described in a configuration script."
            );
        }

        auto res = rigidbody::sequencer::SequenceParser().parse_file(script)->execute();
        fitter::FitReporter::report(res.get()); // what lands on disk is up to the script's own 'save' elements
    } catch (const except::base&) {
        return 1; // our own exceptions already report themselves when thrown
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        return 1;
    }
    return 0;
}
