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
#include <iostream>
#include <string>
#include <vector>

#ifdef _WIN32
    #include <io.h>
    #define AUSAXS_ISATTY(fd) _isatty(fd)
    #define AUSAXS_FILENO(f) _fileno(f)
#else
    #include <unistd.h>
    #define AUSAXS_ISATTY(fd) isatty(fd)
    #define AUSAXS_FILENO(f) fileno(f)
#endif

using namespace ausaxs;

namespace {
    /**
     * @brief Let the user pick one of the GPUs the backend found. Without a terminal to ask on, the first one is picked.
     */
    std::size_t choose_device(const std::vector<std::string>& devices) {
        console::print_text("Found " + std::to_string(devices.size()) + " GPUs, listed with the one most likely to be a dedicated GPU first:");
        for (std::size_t i = 0; i < devices.size(); ++i) {
            console::print_text("  [" + std::to_string(i + 1) + "] " + devices[i]);
        }
        console::print_text("Integrated GPUs are often slower than the CPU, so choose a dedicated one if there is one.");
        if (!AUSAXS_ISATTY(AUSAXS_FILENO(stdin))) {
            console::print_text("Not running in a terminal, so using [1]. Run \"ausaxs gpu install\" in a terminal to choose another.");
            return 0;
        }

        while (true) {
            std::cout << "Which GPU should AUSAXS use? [1]: " << std::flush;
            std::string answer;
            if (!std::getline(std::cin, answer) || answer.empty()) {return 0;}
            try {
                std::size_t used = 0;
                auto choice = std::stoul(answer, &used);
                if (used == answer.size() && 1 <= choice && choice <= devices.size()) {return choice - 1;}
            } catch (const std::exception&) {}
            console::print_text("Please enter a number between 1 and " + std::to_string(devices.size()) + ".");
        }
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

    auto* sub_install = app.add_subcommand("install", "Download and install the GPU backend for this platform, replacing any existing one. If there are several GPUs, asks which one to use.");
    sub_install->add_option("--from", from, "Install from this URL or already downloaded zip file instead of the latest release.");
    sub_install->add_flag("--offline", settings::general::offline.value, "Prevent any network requests. Only --from with a local file works.");

    auto* sub_status = app.add_subcommand("status", "Show the installed GPU backend and the device it finds.");
    auto* sub_remove = app.add_subcommand("remove", "Remove the installed GPU backend.");
    CLI11_PARSE(app, argc, argv);

    try {
        if (*sub_install) {
            if (from.empty()) {
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
