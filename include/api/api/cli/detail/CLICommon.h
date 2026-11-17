// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/ExistingFile.h>
#include <io/Folder.h>

#include <CLI/CLI.hpp>

#include <functional>
#include <optional>
#include <string>

/**
 * @brief The command-line options and setup shared by the structure-based tools (fit, sans, debye, ...).
 *        Each function registers one group of options on a tool's CLI::App; the tool decides which groups it offers.
 */
namespace ausaxs::cli {
    /**
     * @brief The state of the config subcommand. It must outlive the CLI::App it is registered on.
     */
    struct ConfigOptions {
        io::ExistingFile file;
        bool save = false;
        CLI::Option* file_option = nullptr;
    };

    /**
     * @brief Register the top-level options: --output, --licence, --version, --allow-unknown-atoms, --allow-unknown-residues,
     *        --offline, --threads and --gpu.
     */
    void add_general_options(CLI::App& app, const std::string& default_output);

    /**
     * @brief Register the config subcommand: --file, --save and --log.
     */
    void add_config_subcommand(CLI::App& app, ConfigOptions& config, const std::string& log_name);

    /**
     * @brief Register the data subcommand for a measurement file: --qmin, --qmax, --unit, --skip, --rebin and --weighted-bins.
     *
     * @param data_name The kind of data the tool fits, e.g. "SAXS", used in the help text.
     */
    CLI::App* add_data_subcommand(CLI::App& app, const std::string& data_name);

    /**
     * @brief Register the molecule subcommand: --use-occupancy.
     */
    CLI::App* add_molecule_subcommand(CLI::App& app);

    /**
     * @brief Register the hydrogens subcommand: --keep and --force-field-names.
     */
    CLI::App* add_hydrogens_subcommand(CLI::App& app);

    /**
     * @brief Register the neutron solvent options --d2o and --exchange on a subcommand.
     */
    void add_neutron_solvent_options(CLI::App* sub);

    /**
     * @brief Parse the command line, read a settings file if one was given with config --file or is discovered, and then parse
     *        the command line again so that it takes priority over the file.
     *
     * @param discover_dir The folder to look for a settings file in when none was given. Called after the first parse.
     * @return The exit code if the program should stop (e.g. after --help or a parse error), otherwise nothing.
     */
    std::optional<int> parse(
        CLI::App& app, int argc, char const* argv[], const ConfigOptions& config, const std::function<io::Folder()>& discover_dir
    );

    /**
     * @brief Print the neutron solvent and the exchange efficiency the calculation assumes.
     */
    void print_neutron_solvent();
}
