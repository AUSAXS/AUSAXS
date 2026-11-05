// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_gpu.h>

#include <CLI/CLI.hpp>

#include <constants/Version.h>
#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#ifdef _WIN32
    #include <io.h>
#else
    #include <unistd.h>
#endif

using namespace ausaxs;

namespace {
    // whether stdin is an interactive terminal, i.e. whether there is anyone to answer a prompt
    bool interactive() {
        #ifdef _WIN32
            return _isatty(_fileno(stdin)) != 0;
        #else
            return isatty(fileno(stdin)) != 0;
        #endif
    }

    /**
     * @brief Let the user pick one of the GPUs the backend found, or none of them. Without a terminal to ask on, no choice is made.
     */
    std::optional<std::size_t> choose_device(const std::vector<std::string>& devices) {
        console::print_text("");
        console::print_text("");
        if (devices.size() == 1) {console::print_text("Found 1 GPU:");}
        else {console::print_text("Found " + std::to_string(devices.size()) + " GPUs, the one most likely to be dedicated first:");}
        for (int i = 0; i < static_cast<int>(devices.size()); ++i) {
            console::print_text("  [" + std::to_string(i + 1) + "] " + devices[i]);
        }
        console::print_text("  [0] None, use the CPU");
        console::print_text("");
        if (!interactive()) {
            console::print_text("Not running in a terminal, so none is chosen.");
            return std::nullopt;
        }

        std::cout << "Please pick a dedicated GPU from the list [1]: " << std::flush;
        std::string answer;
        if (!std::getline(std::cin, answer)) {return std::nullopt;}
        if (answer.empty()) {return 0;}
        if (answer == "0") {return gpu::GPULoader::no_device;} // checked by hand, since atoi also gives 0 for anything that is not a number
        int choice = std::atoi(answer.c_str());
        if (choice < 1 || static_cast<int>(devices.size()) < choice) {
            throw except::invalid_argument(
                "\"" + answer + "\" is not a number between 0 and " + std::to_string(devices.size()) + ". Run \"ausaxs gpu install\" again to choose one."
            );
        }
        return choice - 1;
    }

    // the exit status: 0 if a backend is loaded and has a usable device, 1 otherwise
    int status() {
        auto folder = gpu::installer::folder();
        console::print_text("Install folder:  " + folder.path());
        if (auto source = gpu::installer::installed_source(); !source.empty()) {
            console::print_text("Installed from:  " + source);
        } else if (!folder.empty() && io::File(folder.path() + "/" + std::string(gpu::GPULoader::library_name())).exists()) {
            console::print_text("Installed from:  incomplete installation; run \"ausaxs gpu install\" again");
        } else {
            console::print_text("Installed from:  not installed");
        }
        if (auto device = gpu::installer::installed_device(); !device.empty()) {
            console::print_text("Chosen GPU:      " + device);
        }
        if (!gpu::installer::installed_library().empty()) {
            console::print_text(std::string("Used by default: ") + (gpu::installer::enabled() ? "yes; turn off with \"ausaxs gpu disable\"" : "no; turn on with \"ausaxs gpu enable\""));
        }
        if (!gpu::installer::supported()) {
            console::print_text("Downloader:      not available in this build");
        }

        // loading the backend reports the device it found, or why there is none
        bool available = gpu::GPULoader::available();
        if (auto path = gpu::GPULoader::path(); !path.empty()) {console::print_text("Loaded from:     " + path);}
        return available ? 0 : 1;
    }
}

int cli_gpu(int argc, char const *argv[]) {
    settings::general::verbose = true;

    std::string from;
    CLI::App app{"Install, check, or remove the GPU backend."};
    app.require_subcommand(1);
    app.add_flag_callback("-v,--version", [] () {console::print_text(constants::version); throw CLI::Success();}, "Print the AUSAXS version.");

    auto* sub_install = app.add_subcommand("install", "Install the GPU backend for this platform, unless a compatible one is already installed. If there are several GPUs, asks which one to use, so run it again to choose another.");
    sub_install->add_option("--from", from, "Install from this URL or already downloaded zip file instead of the latest release, replacing any existing installation.");
    sub_install->add_flag("--offline", settings::general::offline.value, "Prevent any network requests. Only --from with a local file works.");

    auto* sub_status = app.add_subcommand("status", "Show the installed GPU backend and the device it finds.");
    auto* sub_remove = app.add_subcommand("remove", "Remove the installed GPU backend.");
    auto* sub_enable = app.add_subcommand("enable", "Use the installed GPU backend by default. Choosing a usable GPU with \"ausaxs gpu install\" already does this.");
    auto* sub_disable = app.add_subcommand("disable", "Stop using the installed GPU backend by default, without removing it. It can still be used with --gpu.");
    CLI11_PARSE(app, argc, argv);

    try {
        if (*sub_install) {
            if (from.empty()) {
                if (gpu::installer::select_device(choose_device)) {return 0;}
                if (gpu::installer::asset_name().empty()) {
                    console::print_warning("No prebuilt GPU backend exists for this platform.");
                    return 1;
                }
                from = gpu::installer::asset_url();
            }
            gpu::installer::install(from, choose_device);
            return 0;
        }
        if (*sub_status) {
            return status();
        }
        if (*sub_enable || *sub_disable) {
            bool enable = sub_enable->parsed();
            gpu::installer::set_enabled(enable);
            console::print_success(std::string("The GPU is now ") + (enable ? "used by default." : "only used when asked for with --gpu."));
            return 0;
        }
        if (*sub_remove) {
            if (gpu::installer::remove()) {console::print_success("Removed the GPU backend from " + gpu::installer::folder().path());}
            else {console::print_text("No GPU backend is installed in " + gpu::installer::folder().path());}
            return 0;
        }
    } catch (const except::base&) {
        return 1; // already printed when it was thrown
    } catch (const std::exception& e) {
        console::print_warning(e.what());
        return 1;
    }
    return 0;
}
