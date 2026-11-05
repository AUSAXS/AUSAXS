// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <gpu/GPUInstaller.h>

#include <gpu/GPULoader.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Curl.h>
#include <utility/Exceptions.h>

#ifdef AUSAXS_GPU_DOWNLOADER
    #include <miniz.h> // NOLINT
#endif

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::gpu;
namespace fs = std::filesystem;

namespace {
    constexpr const char* release_url = "https://github.com/AUSAXS/AUSAXS-GPU/releases/";

    // Written last, once the backend has loaded. Its absence marks an installation that did not finish, which the loader then ignores.
    // Holds the source on the first line, and the identifier of the chosen GPU on the second, which is empty if none was found.
    constexpr const char* marker_name = "installed.txt";

    io::File marker_file(const io::Folder& folder) {
        return {folder.path() + "/" + marker_name};
    }

    io::File library_file(const io::Folder& folder) {
        return {folder.path() + "/" + std::string(GPULoader::library_name())};
    }

    // Beside the install folder rather than inside it, so it is never part of what it unpacks.
    io::Folder staging_folder(const io::Folder& folder) {
        return {folder.path() + ".staging"};
    }

    // Its presence turns the GPU on by default. Created once a usable GPU has been chosen, and kept in the cache rather than the install folder
    // so turning the GPU off and on never touches the installation itself.
    io::File toggle_file(const io::Folder& folder) {
        return {(fs::path(folder.path()).parent_path()/".autotoggle_gpu").string()};
    }

    // Lists every path an installation put into the install folder, relative to it, one per line. Only what is listed is ever deleted, so
    // anything else in the folder is left alone. It is kept beside the folder rather than inside it, so it is never part of what it lists.
    io::File record_file(const io::Folder& folder) {
        return {folder.path() + ".record"};
    }

    /**
     * @brief The paths listed in @a record, or none if there is no record.
     */
    std::vector<fs::path> read_record(const io::File& record) {
        std::vector<fs::path> entries;
        std::ifstream in(record.path());
        for (std::string line; std::getline(in, line);) {
            if (line.empty()) {continue;}
            entries.push_back(fs::path(line).lexically_normal());
        }
        return entries;
    }

    void write_record(const io::File& record, const std::vector<fs::path>& entries) {
        std::string contents;
        for (const auto& entry : entries) {contents += entry.generic_string() + "\n";}
        record.create(contents);
    }

    /**
     * @brief Every path inside @a folder, relative to it. Symbolic links are listed, but not followed.
     */
    std::vector<fs::path> list_contents(const io::Folder& folder) {
        std::vector<fs::path> entries;
        fs::path base = folder.path();
        for (const auto& entry : fs::recursive_directory_iterator(base)) {
            entries.push_back(entry.path().lexically_relative(base));
        }
        return entries;
    }

    /**
     * @brief Delete the paths @a entries lists from @a folder, and nothing else.
     *
     * A listed directory is deleted only once it is empty, so one the user has since added files to is kept, along with those files.
     * The folder itself is always kept.
     *
     * @throws except::runtime_error if a listed path could not be deleted.
     */
    void uninstall(const io::Folder& folder, std::vector<fs::path> entries) {
        if (!folder.exists()) {return;} // nothing to delete
        fs::path base = folder.path();

        // a directory sorts before everything inside it, so in reverse order its contents are deleted first
        std::ranges::sort(entries, std::greater<>());
        std::error_code ec;
        for (const auto& entry : entries) {
            fs::path target = base/entry;
            auto status = fs::symlink_status(target, ec);
            if (ec || !fs::exists(status)) {continue;} // already gone
            if (fs::is_directory(status) && !fs::is_empty(target, ec)) {continue;}
            fs::remove(target, ec);
            if (ec) {
                throw except::runtime_error(
                    "gpu::installer: could not delete \"" + target.string() + "\": " + ec.message() + "."
                );
            }
        }
    }

    /**
     * @brief Delete @a folder if it is empty.
     *
     * @return Whether it is gone.
     */
    bool remove_if_empty(const io::Folder& folder) {
        if (!folder.exists()) {return true;}
        std::error_code ec;
        if (!fs::is_empty(folder.path(), ec) || ec) {return false;}
        return fs::remove(folder.path(), ec);
    }

    // Best effort, for cleaning up after a failure that is already being reported.
    template<typename T>
    void discard(const T& path) noexcept {
        try {path.remove();} catch (...) {}
    }

    // Best effort: delete an installation that failed, along with its record.
    void discard_installation(const io::Folder& folder, const std::vector<fs::path>& entries, const io::File& record) noexcept {
        try {
            uninstall(folder, entries);
            remove_if_empty(folder);
            record.remove();
        } catch (...) {}
    }

    #ifdef AUSAXS_GPU_DOWNLOADER
        std::string zip_error(mz_zip_archive& archive) {
            return mz_zip_get_error_string(mz_zip_get_last_error(&archive));
        }

        /**
         * @brief Unpack @a zip into @a destination. Each entry's CRC is checked as it is unpacked, which catches a corrupt or truncated download.
         */
        void extract(const io::File& zip, const io::Folder& destination) {
            mz_zip_archive archive{};
            if (!static_cast<bool>(mz_zip_reader_init_file(&archive, zip.path().c_str(), 0))) {
                throw except::runtime_error("gpu::installer: \"" + zip.path() + "\" is not a readable zip archive: " + zip_error(archive));
            }
            std::unique_ptr<mz_zip_archive, mz_bool(*)(mz_zip_archive*)> guard(&archive, &mz_zip_reader_end);

            mz_uint count = mz_zip_reader_get_num_files(&archive);
            for (mz_uint i = 0; i < count; ++i) {
                mz_zip_archive_file_stat stat;
                if (!static_cast<bool>(mz_zip_reader_file_stat(&archive, i, &stat))) {
                    throw except::runtime_error("gpu::installer: could not read entry " + std::to_string(i) + " of the archive: " + zip_error(archive));
                }

                // an entry must not reach outside the destination, whatever the archive says
                fs::path relative = fs::path(stat.m_filename).lexically_normal();
                if (relative.empty() || relative.has_root_path() || *relative.begin() == "..") {
                    throw except::runtime_error("gpu::installer: the archive contains an unsafe path: \"" + std::string(stat.m_filename) + "\"");
                }
                fs::path target = fs::path(destination.path())/relative;
                if (static_cast<bool>(mz_zip_reader_is_file_a_directory(&archive, i))) {
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
                if (!static_cast<bool>(mz_zip_reader_extract_to_callback(&archive, i, write, &out, 0))) {
                    throw except::runtime_error("gpu::installer: could not unpack \"" + std::string(stat.m_filename) + "\": " + zip_error(archive));
                }
                out.close();
                if (!out) {throw except::runtime_error("gpu::installer: could not write \"" + target.string() + "\"");}

                #ifndef _WIN32
                    // the runtime runs the bundled LLVM tools to compile kernels, so their executable bits must survive
                    constexpr int unix_host = 3;
                    if ((stat.m_version_made_by >> 8) == unix_host && static_cast<bool>(((stat.m_external_attr >> 16) & 0111))) {
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

std::string installer::asset_url() {
    // GitHub redirects latest/download/ to the newest release, which saves a query of its API
    return release_url + std::string("latest/download/") + asset_name();
}

io::Folder installer::folder() {
    return {settings::general::helper::gpu_folder()};
}

io::File installer::installed_library() {
    io::Folder folder = installer::folder();
    io::File library = library_file(folder);
    if (!marker_file(folder).exists() || !library.exists()) {return {};}
    return library;
}

namespace {
    // Marks the installation in @a folder complete, which must only be done once its backend has loaded.
    void write_marker(const io::Folder& folder, const std::string& source, const std::string& device) {
        marker_file(folder).create(source + "\n" + device + "\n");
    }

    // Report the outcome of testing a backend that loaded. @a what starts the message, e.g. "Installed the GPU backend".
    void report(const std::string& what, const GPULoader::Probe& probe) {
        if (probe.chosen && probe.selected.empty()) {
            console::print_text("No GPU will be used. The backend is kept, so run \"ausaxs gpu install\" again if you change your mind.");
            return;
        }
        if (probe.available) {console::print_success(what + ". Device: " + probe.device);}
        else {
            console::print_warning(
                what + ", but the backend found no usable device on this machine" +
                (probe.error.empty() ? std::string() : " (" + probe.error + ")") + ".\n"
                "If this machine has a GPU, it is usually missing a driver. See " + std::string(GPULoader::troubleshooting_url) + " for common fixes."
            );
        }
        if (!installer::enabled()) {
            console::print_text("The GPU is not used by default. Choose a usable one with \"ausaxs gpu install\" in a terminal, or pass --gpu.");
        }
    }

    // Line @a index of the marker of a complete installation, or empty if there is none.
    std::string marker_line(int index) {
        if (installer::installed_library().empty()) {return {};}
        std::ifstream in(marker_file(installer::folder()).path());
        std::string line;
        for (int i = 0; i <= index; ++i) {
            if (!std::getline(in, line)) {return {};}
        }
        return line;
    }
}

std::string installer::installed_source() {
    return marker_line(0);
}

std::string installer::installed_device() {
    return marker_line(1);
}

void installer::install([[maybe_unused]] const std::string& source, [[maybe_unused]] const GPULoader::DeviceChooser& choose) {
    #ifndef AUSAXS_GPU_DOWNLOADER
        throw except::runtime_error("gpu::installer: this build of AUSAXS cannot install GPU backends, since it was built without the downloader (GPU_DOWNLOADER=OFF)");
    #else
        io::Folder folder = installer::folder();
        io::File record = record_file(folder);
        std::vector<fs::path> previous = read_record(record);

        io::File zip(source);
        bool downloaded = false;
        if (!zip.exists()) {
            if (settings::general::offline) {throw except::runtime_error("gpu::installer: cannot download \"" + source + "\" in offline mode");}
            zip = io::File(folder.path() + ".download.zip");
            zip.directory().create();
            console::print_text("Downloading " + source);
            if (!curl::download(source, zip, true)) {
                throw except::runtime_error("gpu::installer: could not download the GPU backend from \"" + source + "\"");
            }
            downloaded = true;
        }

        io::Folder staging = staging_folder(folder);
        std::vector<fs::path> entries;
        try {
            staging.remove(); // left behind by an interrupted install
            console::print_text("Unpacking into " + folder.path());
            extract(zip, staging);
            if (!library_file(staging).exists()) {
                throw except::runtime_error(
                    "gpu::installer: the archive does not contain " + std::string(GPULoader::library_name()) + "; it is not a GPU backend for this platform"
                );
            }
            entries = list_contents(staging);
            entries.emplace_back(marker_name);

            // Only what an earlier installation put there may be replaced. Checked before anything is changed, so a refusal changes nothing.
            for (const auto& entry : entries) {
                if (entry.has_parent_path() || std::ranges::find(previous, entry) != previous.end()) {continue;}
                std::error_code ec;
                if (fs::exists(fs::symlink_status(fs::path(folder.path())/entry, ec))) {
                    throw except::runtime_error(
                        "gpu::installer: refusing to install into \"" + folder.path() + "\": it already holds \"" + entry.string() + "\", "
                        "which is not from an earlier installation"
                    );
                }
            }
        } catch (...) {
            discard(staging);
            if (downloaded) {discard(zip);}
            throw;
        }
        if (downloaded) {discard(zip);}

        // The new backend is tested where it will live, since a loaded library cannot be moved on Windows. That means the old installation
        // goes first, and a failed test leaves none at all.
        try {
            uninstall(folder, previous);
        } catch (...) {
            discard(staging);
            throw;
        }
        discard(toggle_file(folder)); // the GPU chosen for the old installation no longer applies

        // Recorded before anything is moved in, so an interrupted install can still be removed.
        write_record(record, entries);
        try {
            folder.create();
            for (const auto& entry : entries) {
                if (entry.has_parent_path() || entry == marker_name) {continue;}
                fs::rename(fs::path(staging.path())/entry, fs::path(folder.path())/entry);
            }
        } catch (const fs::filesystem_error& e) {
            discard(staging);
            throw except::runtime_error(
                "gpu::installer: could not move the backend into \"" + folder.path() + "\": " + e.code().message() +
                ". If an earlier installation left files of yours behind, move them elsewhere and try again."
            );
        }
        discard(staging);

        console::print_text("Testing the backend. This may take a few seconds while the kernel is compiled.");
        // The chooser is only asked once the backend has loaded and listed its GPUs. From then on the installation is kept without a chosen
        // GPU, even if no choice is made because the chooser throws or the process is interrupted at a prompt, so choosing one later does
        // not need it downloaded again.
        GPULoader::DeviceChooser keep_then_choose;
        if (choose) {
            keep_then_choose = [&] (const std::vector<std::string>& devices) {
                write_marker(folder, source, "");
                return choose(devices);
            };
        }
        auto probe = GPULoader::probe(library_file(folder).path(), keep_then_choose);
        if (!probe.loaded) {
            discard_installation(folder, entries, record);
            throw except::runtime_error(
                "gpu::installer: the GPU backend could not be loaded, and has been removed again: " + probe.error +
                "\nSee " + std::string(GPULoader::troubleshooting_url) + " for common fixes."
            );
        }
        write_marker(folder, source, probe.chosen ? probe.selected : "");
        set_enabled(probe.chosen && probe.available); // only a GPU someone chose, since the first may be integrated graphics
        report("Installed the GPU backend", probe);
    #endif
}

bool installer::select_device(const GPULoader::DeviceChooser& choose) {
    io::File library = installed_library();
    if (library.empty()) {return false;}
    std::string source = installed_source();

    console::print_text("Testing the installed backend. This may take a few seconds while the kernel is compiled.");
    auto probe = GPULoader::probe(library.path(), choose); // a throwing chooser leaves the installation as it was
    if (!probe.loaded) {
        console::print_text("The installed GPU backend cannot be used: " + probe.error);
        return false;
    }
    if (probe.chosen) { // otherwise the earlier choice stands
        write_marker(folder(), source, probe.selected); // empty if none was chosen
        set_enabled(probe.available);
    }
    report(probe.chosen ? "Chose the GPU" : "Tested the installed GPU backend", probe);
    return true;
}

bool installer::enabled() {
    return toggle_file(folder()).exists() && !installed_library().empty();
}

void installer::set_enabled(bool enable) {
    io::File toggle = toggle_file(folder());
    if (!enable) {
        if (toggle.exists()) {toggle.remove();}
        return;
    }
    if (installed_library().empty()) {throw except::runtime_error("gpu::installer: cannot turn the GPU on, since no GPU backend is installed");}
    if (installed_device().empty()) {
        throw except::runtime_error("gpu::installer: cannot turn the GPU on, since no GPU has been chosen; choose one with \"ausaxs gpu install\"");
    }
    toggle.create();
}

bool installer::remove() {
    io::Folder folder = installer::folder();
    io::File record = record_file(folder);
    discard(staging_folder(folder));
    discard(toggle_file(folder));
    if (!record.exists()) {
        if (remove_if_empty(folder)) {return false;}
        throw except::runtime_error(
            "gpu::installer: refusing to remove \"" + folder.path() + "\": there is no record of what was installed there (\"" + record.path() +
            "\"). Remove it by hand."
        );
    }

    uninstall(folder, read_record(record));
    record.remove();
    if (!remove_if_empty(folder)) {
        console::print_warning("Kept \"" + folder.path() + "\", since it holds files that were not part of the GPU backend.");
    }
    return true;
}
