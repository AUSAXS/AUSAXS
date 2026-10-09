// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <filesystem>
#include <string>

/**
 * @brief Downloading and installing a prebuilt GPU backend.
 *
 * Backends are published as one zip per platform on the AUSAXS-GPU releases page, each holding the backend and the AdaptiveCpp runtime it
 * needs. An installation is unpacked into settings::general::gpu_folder, and marked complete only once it has been loaded successfully,
 * so the loader never picks up a half-written or broken one.
 *
 * Unpacking needs miniz. Builds without it (GPU_DOWNLOADER=OFF) can still report on and remove an installation, but not make one.
 */
namespace ausaxs::gpu::installer {
    /**
     * @brief Whether this build can download and install a backend.
     */
    bool supported();

    /**
     * @brief Name of the release asset for this platform, or empty if no prebuilt backend exists for it.
     */
    std::string asset_name();

    /**
     * @brief The release a backend is downloaded from by default: the one matching this version of AUSAXS.
     */
    std::string default_tag();

    /**
     * @brief URL of the release asset for this platform from release @a tag.
     */
    std::string asset_url(const std::string& tag);

    /**
     * @brief The folder a backend is installed to.
     */
    std::filesystem::path folder();

    /**
     * @brief Path of the installed backend library, or empty if there is no complete installation.
     */
    std::filesystem::path installed_library();

    /**
     * @brief Where the installed backend came from, or empty if there is no complete installation.
     */
    std::string installed_source();

    /**
     * @brief Download and install a backend, replacing any existing installation.
     *
     * The backend is test-loaded before the installation is marked complete. A backend that loads but finds no device is still installed,
     * since it may be installed on one machine and used on another, e.g. a login and a compute node.
     *
     * @param source A URL, or the path of an already downloaded zip.
     * @throws except::runtime_error if the backend could not be installed. Nothing is left behind in that case.
     */
    void install(const std::string& source);

    /**
     * @brief Remove an installed backend.
     *
     * @return False if there was nothing to remove.
     * @throws except::runtime_error if the folder does not look like an installation, or could not be removed.
     */
    bool remove();
}
