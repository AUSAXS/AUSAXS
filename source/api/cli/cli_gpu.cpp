// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/cli/cli_gpu.h>

#include <CLI/CLI.hpp>

#include <constants/Version.h>
#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <string>

using namespace ausaxs;

namespace {
    // the exit status: 0 if a backend is loaded and has a usable device, 1 otherwise
    int status() {
        auto folder = gpu::installer::folder();
        console::print_text("Install folder:  " + folder.string());
        if (auto source = gpu::installer::installed_source(); !source.empty()) {
            console::print_text("Installed from:  " + source);
        } else if (std::filesystem::exists(folder/gpu::GPULoader::library_name())) {
            console::print_text("Installed from:  incomplete installation; run \"ausaxs gpu install\" again");
        } else {
            console::print_text("Installed from:  not installed");
        }
        if (!gpu::installer::supported()) {
            console::print_text("Downloader:      not available in this build (GPU_DOWNLOADER=OFF)");
        }

        // loading the backend reports the device it found, or why there is none
        bool available = gpu::GPULoader::available();
        if (auto path = gpu::GPULoader::path(); !path.empty()) {console::print_text("Loaded from:     " + path);}
        return available ? 0 : 1;
    }
}

int cli_gpu(int argc, char const *argv[]) {
    settings::general::verbose = true;

    std::string tag = gpu::installer::default_tag();
    std::string from;
    CLI::App app{"Install, check, or remove the GPU backend."};
    app.require_subcommand(1);
    app.add_flag_callback("-v,--version", [] () {console::print_text(constants::version); exit(0);}, "Print the AUSAXS version.");

    auto* sub_install = app.add_subcommand("install", "Download and install the GPU backend for this platform, replacing any existing one.");
    auto* p_tag = sub_install->add_option("--tag", tag, "The release to download the backend from.")->default_val(tag);
    sub_install->add_option("--from", from, "Install from this URL or already downloaded zip file instead.")->excludes(p_tag);
    sub_install->add_flag("--offline", settings::general::offline, "Prevent any network requests. Only --from with a local file works.");

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
                from = gpu::installer::asset_url(tag);
            }
            gpu::installer::install(from);
            return 0;
        }
        if (*sub_status) {
            return status();
        }
        if (*sub_remove) {
            if (gpu::installer::remove()) {console::print_success("Removed the GPU backend from " + gpu::installer::folder().string());}
            else {console::print_text("No GPU backend is installed in " + gpu::installer::folder().string());}
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
