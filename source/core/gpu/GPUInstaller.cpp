// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>

#include <constants/Version.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Curl.h>
#include <utility/Exceptions.h>

#ifdef AUSAXS_GPU_DOWNLOADER
    #include <miniz.h>
#endif

#include <fstream>
#include <memory>

using namespace ausaxs;
using namespace ausaxs::gpu;
namespace fs = std::filesystem;

namespace {
    constexpr const char* release_url = "https://github.com/AUSAXS/AUSAXS-GPU/releases/download/";

    // Written last, once the backend has loaded. Its absence marks an installation that did not finish, which the loader then ignores.
    constexpr const char* marker_name = "installed.txt";

    /**
     * @brief Whether @a folder holds an installation, finished or not.
     *
     * The folder is the user's to configure, so nothing is ever replaced or deleted unless it passes this check.
     */
    bool looks_installed(const fs::path& folder) {
        std::error_code ec;
        return fs::exists(folder/marker_name, ec) || fs::exists(folder/GPULoader::library_name(), ec);
    }

    fs::path folder_or_throw() {
        fs::path folder = installer::folder();
        if (folder.empty()) {throw except::runtime_error("gpu::installer: no install folder is configured (settings: gpu_folder)");}
        return folder;
    }

    // Beside the install folder rather than inside it, so it can be moved into place in one rename.
    fs::path staging_folder(const fs::path& folder) {
        return folder.parent_path()/(folder.filename().string() + ".staging");
    }

    #ifdef AUSAXS_GPU_DOWNLOADER
        std::string zip_error(mz_zip_archive& archive) {
            return mz_zip_get_error_string(mz_zip_get_last_error(&archive));
        }

        /**
         * @brief Unpack @a zip into @a destination. Each entry's CRC is checked as it is unpacked, which catches a corrupt or truncated download.
         */
        void extract(const fs::path& zip, const fs::path& destination) {
            mz_zip_archive archive{};
            if (!mz_zip_reader_init_file(&archive, zip.string().c_str(), 0)) {
                throw except::runtime_error("gpu::installer: \"" + zip.string() + "\" is not a readable zip archive: " + zip_error(archive));
            }
            std::unique_ptr<mz_zip_archive, mz_bool(*)(mz_zip_archive*)> guard(&archive, &mz_zip_reader_end);

            mz_uint count = mz_zip_reader_get_num_files(&archive);
            for (mz_uint i = 0; i < count; ++i) {
                mz_zip_archive_file_stat stat;
                if (!mz_zip_reader_file_stat(&archive, i, &stat)) {
                    throw except::runtime_error("gpu::installer: could not read entry " + std::to_string(i) + " of the archive: " + zip_error(archive));
                }

                // an entry must not reach outside the destination, whatever the archive says
                fs::path relative = fs::path(stat.m_filename).lexically_normal();
                if (relative.empty() || relative.has_root_path() || *relative.begin() == "..") {
                    throw except::runtime_error("gpu::installer: the archive contains an unsafe path: \"" + std::string(stat.m_filename) + "\"");
                }
                fs::path target = destination/relative;
                if (mz_zip_reader_is_file_a_directory(&archive, i)) {
                    fs::create_directories(target);
                    continue;
                }
                fs::create_directories(target.parent_path());

                std::ofstream out(target, std::ios::binary);
                if (!out) {throw except::runtime_error("gpu::installer: could not create \"" + target.string() + "\"");}
                auto write = [] (void* opaque, mz_uint64, const void* buffer, size_t n) -> size_t {
                    auto& out = *static_cast<std::ofstream*>(opaque);
                    out.write(static_cast<const char*>(buffer), static_cast<std::streamsize>(n));
                    return out ? n : 0;
                };
                if (!mz_zip_reader_extract_to_callback(&archive, i, write, &out, 0)) {
                    throw except::runtime_error("gpu::installer: could not unpack \"" + std::string(stat.m_filename) + "\": " + zip_error(archive));
                }
                out.close();
                if (!out) {throw except::runtime_error("gpu::installer: could not write \"" + target.string() + "\"");}

                #ifndef _WIN32
                    // the runtime runs the bundled LLVM tools to compile kernels, so their executable bits must survive
                    constexpr int unix_host = 3;
                    if ((stat.m_version_made_by >> 8) == unix_host && ((stat.m_external_attr >> 16) & 0111)) {
                        fs::permissions(target, fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec, fs::perm_options::add);
                    }
                #endif
            }
        }
    #endif
}

bool installer::supported() {
    #ifdef AUSAXS_GPU_DOWNLOADER
        return true;
    #else
        return false;
    #endif
}

std::string installer::asset_name() {
    #if defined(_WIN32) && (defined(_M_X64) || defined(__x86_64__))
        return "ausaxs-gpu-windows-x64.zip";
    #elif defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
        return "ausaxs-gpu-macos-arm64.zip";
    #elif defined(__linux__) && defined(__x86_64__)
        return "ausaxs-gpu-linux-x64.zip";
    #else
        return "";
    #endif
}

std::string installer::default_tag() {
    return std::string(constants::version);
}

std::string installer::asset_url(const std::string& tag) {
    return release_url + tag + "/" + asset_name();
}

fs::path installer::folder() {
    if (settings::general::gpu_folder.empty()) {return {};}
    fs::path folder = fs::path(settings::general::gpu_folder).lexically_normal();
    if (!folder.has_filename()) {folder = folder.parent_path();} // trailing separator
    return folder;
}

fs::path installer::installed_library() {
    fs::path folder = installer::folder();
    if (folder.empty()) {return {};}
    std::error_code ec;
    fs::path library = folder/GPULoader::library_name();
    if (!fs::exists(folder/marker_name, ec) || !fs::exists(library, ec)) {return {};}
    return library;
}

std::string installer::installed_source() {
    if (installed_library().empty()) {return {};}
    std::ifstream in(folder()/marker_name);
    std::string source;
    std::getline(in, source);
    return source;
}

void installer::install([[maybe_unused]] const std::string& source) {
    #ifndef AUSAXS_GPU_DOWNLOADER
        throw except::runtime_error("gpu::installer: this build of AUSAXS cannot install GPU backends, since it was built without the downloader (GPU_DOWNLOADER=OFF)");
    #else
        fs::path folder = folder_or_throw();
        std::error_code ec;
        if (fs::exists(folder) && !fs::is_empty(folder) && !looks_installed(folder)) {
            throw except::runtime_error(
                "gpu::installer: refusing to install into \"" + folder.string() + "\": it is not empty, and does not hold a GPU backend"
            );
        }

        fs::path zip = source;
        bool downloaded = false;
        if (!fs::is_regular_file(zip, ec)) {
            if (settings::general::offline) {throw except::runtime_error("gpu::installer: cannot download \"" + source + "\" in offline mode");}
            fs::create_directories(folder.parent_path());
            zip = folder.parent_path()/(folder.filename().string() + ".download.zip");
            console::print_text("Downloading " + source);
            if (!curl::download(source, zip.string(), true)) {
                throw except::runtime_error("gpu::installer: could not download the GPU backend from \"" + source + "\"");
            }
            downloaded = true;
        }

        fs::path staging = staging_folder(folder);
        fs::remove_all(staging, ec);
        try {
            console::print_text("Unpacking into " + folder.string());
            extract(zip, staging);
            if (!fs::exists(staging/GPULoader::library_name())) {
                throw except::runtime_error(
                    "gpu::installer: the archive does not contain " + std::string(GPULoader::library_name()) + "; it is not a GPU backend for this platform"
                );
            }
        } catch (...) {
            fs::remove_all(staging, ec);
            if (downloaded) {fs::remove(zip, ec);}
            throw;
        }
        if (downloaded) {fs::remove(zip, ec);}

        // The new backend is tested where it will live, since a loaded library cannot be moved on Windows. That means the old installation
        // goes first, and a failed test leaves none at all.
        if (fs::exists(folder)) {
            fs::remove_all(folder, ec);
            if (ec) {
                fs::remove_all(staging, ec);
                throw except::runtime_error(
                    "gpu::installer: could not remove the existing installation at \"" + folder.string() + "\": " + ec.message() +
                    ". Is another AUSAXS process using it?"
                );
            }
        }
        fs::rename(staging, folder, ec);
        if (ec) {
            fs::remove_all(staging, ec);
            throw except::runtime_error("gpu::installer: could not move the backend into \"" + folder.string() + "\": " + ec.message());
        }

        console::print_text("Testing the backend. This may take a while the first time, as the kernels are compiled.");
        auto probe = GPULoader::probe((folder/GPULoader::library_name()).string());
        if (!probe.loaded) {
            fs::remove_all(folder, ec);
            throw except::runtime_error("gpu::installer: the GPU backend could not be loaded, and has been removed again: " + probe.error);
        }
        std::ofstream marker(folder/marker_name);
        marker << source << "\n";
        if (!marker) {throw except::runtime_error("gpu::installer: could not write \"" + (folder/marker_name).string() + "\"");}

        if (probe.available) {console::print_success("Installed the GPU backend. Device: " + probe.device);}
        else {
            console::print_warning(
                "Installed the GPU backend, but it found no usable device on this machine (" + probe.error + "). "
                "It will be used on machines that have one."
            );
        }
    #endif
}

bool installer::remove() {
    fs::path folder = folder_or_throw();
    std::error_code ec;
    fs::remove_all(staging_folder(folder), ec);
    if (!fs::exists(folder)) {return false;}
    if (!looks_installed(folder)) {
        throw except::runtime_error("gpu::installer: refusing to remove \"" + folder.string() + "\": it does not hold a GPU backend");
    }
    fs::remove_all(folder, ec);
    if (ec) {
        throw except::runtime_error(
            "gpu::installer: could not remove \"" + folder.string() + "\": " + ec.message() + ". Is another AUSAXS process using it?"
        );
    }
    return true;
}
