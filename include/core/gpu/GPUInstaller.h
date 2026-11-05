// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <gpu/GPULoader.h>
#include <io/File.h>

#include <string>

/**
 * @brief Downloading and installing a prebuilt GPU backend.
 *
 * Backends are published as one zip per platform on the AUSAXS-GPU releases page, each holding the backend and the AdaptiveCpp runtime it
 * needs. An installation is unpacked into the gpu/ subfolder of the cache (settings::general::helper::gpu_folder), and marked complete only
 * once it has been loaded successfully, so the loader never picks up a half-written or broken one. Every path it unpacks is recorded beside
 * that folder, and only recorded paths are ever replaced or deleted, so anything else in the folder is left alone.
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
     * @brief URL of the release asset for this platform from the latest release.
     */
    std::string asset_url();

    /**
     * @brief The folder a backend is installed to.
     */
    io::Folder folder();

    /**
     * @brief The installed backend library, or empty if there is no complete installation.
     */
    io::File installed_library();

    /**
     * @brief Where the installed backend came from, or empty if there is no complete installation.
     */
    std::string installed_source();

    /**
     * @brief Identifier of the GPU chosen during the installation, or empty if there is no complete installation, or it found no GPU to choose.
     */
    std::string installed_device();

    /**
     * @brief Download and install a backend, replacing any existing installation.
     *
     * The backend is test-loaded before the installation is marked complete. A backend that loads but finds no device is still installed,
     * since it may be installed on one machine and used on another, e.g. a login and a compute node. If it finds a GPU, the one tested is
     * recorded as the one to use from then on.
     *
     * @param source A URL, or the path of an already downloaded zip.
     * @param choose Asked to pick a GPU whenever the backend finds one. Only a GPU it picks is recorded, and only if that GPU is usable is it
     *               used by default, as by set_enabled(true); the first may well be integrated graphics, which is slower than the CPU. Without
     *               it, or if it makes no choice, the backend is still installed, just not used by default. The same happens if it picks
     *               GPULoader::no_device, in which case nothing is tested, if it throws, in which case the exception is rethrown, or if the
     *               process is interrupted while it is asked.
     * @throws except::runtime_error if the backend could not be installed, or the folder holds something that is not from an earlier
     *         installation. Nothing is left behind in that case.
     */
    void install(const std::string& source, const GPULoader::DeviceChooser& choose = {});

    /**
     * @brief Choose the GPU the installed backend uses, testing it like install() does.
     *
     * @param choose As for install(). If it makes no choice or throws, the installation is left as it was. If it picks
     *               GPULoader::no_device, the GPU chosen earlier is forgotten and the GPU is no longer used by default, like remove(), except
     *               that the backend is kept, so a GPU can be chosen again without downloading it.
     * @return False, with nothing changed, if there is no installation or it cannot be loaded, e.g. because it is from a release that
     *         implements a different ABI version. A fresh installation is needed then.
     */
    bool select_device(const GPULoader::DeviceChooser& choose = {});

    /**
     * @brief Whether the GPU is used by default, i.e. a backend is installed and the GPU has not been turned off with set_enabled(false).
     */
    bool enabled();

    /**
     * @brief Turn the default use of the GPU on or off. Only a file in the cache is changed, never the installation.
     *
     * @throws except::runtime_error if turning it on while no backend is installed, or no GPU has been chosen for it.
     */
    void set_enabled(bool enable);

    /**
     * @brief Remove an installed backend.
     *
     * Only the recorded paths are deleted. The folder itself is deleted only if nothing else is left in it.
     *
     * @return False if there was nothing to remove.
     * @throws except::runtime_error if the folder exists but there is no record of an installation, or a recorded path could not be deleted.
     */
    bool remove();
}
