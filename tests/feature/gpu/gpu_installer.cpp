// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>
#include <settings/GeneralSettings.h>

#include <filesystem>
#include <fstream>

using namespace ausaxs;
namespace fs = std::filesystem;

namespace {
    struct TempFolder {
        TempFolder() : previous(settings::general::gpu_folder) {
            fs::remove_all(root);
            settings::general::gpu_folder = (root/"gpu").string() + "/";
        }
        ~TempFolder() {
            settings::general::gpu_folder = previous;
            fs::remove_all(root);
        }
        fs::path root = "temp/tests/gpu_installer";
        std::string previous;
    };

    void touch(const fs::path& file) {
        fs::create_directories(file.parent_path());
        std::ofstream(file) << "content";
    }
}

TEST_CASE("installer::folder") {
    TempFolder temp;
    CHECK(gpu::installer::folder() == temp.root/"gpu");

    settings::general::gpu_folder = (temp.root/"gpu").string();
    CHECK(gpu::installer::folder() == temp.root/"gpu");

    settings::general::gpu_folder = "";
    CHECK(gpu::installer::folder().empty());
    CHECK(gpu::installer::installed_library().empty());
}

TEST_CASE("installer::asset_url") {
    if (gpu::installer::asset_name().empty()) {SKIP("no prebuilt backend for this platform");}
    auto url = gpu::installer::asset_url("v1.2.3");
    CHECK(url.starts_with("https://"));
    CHECK(url.ends_with("/v1.2.3/" + gpu::installer::asset_name()));
}

TEST_CASE("installer::installed_library") {
    TempFolder temp;
    auto folder = gpu::installer::folder();
    auto library = folder/gpu::GPULoader::library_name();

    SECTION("nothing installed") {
        CHECK(gpu::installer::installed_library().empty());
        CHECK(gpu::installer::installed_source().empty());
    }

    SECTION("an unfinished installation is ignored") {
        touch(library);
        CHECK(gpu::installer::installed_library().empty());
    }

    SECTION("a finished installation is found") {
        touch(library);
        std::ofstream(folder/"installed.txt") << "https://example.org/backend.zip\n";
        CHECK(gpu::installer::installed_library() == library);
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
    }
}

TEST_CASE("installer::remove") {
    TempFolder temp;
    auto folder = gpu::installer::folder();

    SECTION("nothing to remove") {
        CHECK_FALSE(gpu::installer::remove());
    }

    SECTION("an installation is removed") {
        touch(folder/gpu::GPULoader::library_name());
        touch(folder/"hipSYCL"/"plugin");
        CHECK(gpu::installer::remove());
        CHECK_FALSE(fs::exists(folder));
    }

    SECTION("a folder that is not an installation is left alone") {
        touch(folder/"unrelated.txt");
        CHECK_THROWS(gpu::installer::remove());
        CHECK(fs::exists(folder/"unrelated.txt"));
    }
}

TEST_CASE("installer::install") {
    TempFolder temp;
    auto folder = gpu::installer::folder();
    if (!gpu::installer::supported()) {
        CHECK_THROWS(gpu::installer::install("backend.zip"));
        return;
    }

    SECTION("a folder that is not an installation is left alone") {
        touch(folder/"unrelated.txt");
        CHECK_THROWS(gpu::installer::install("backend.zip"));
        CHECK(fs::exists(folder/"unrelated.txt"));
    }

    SECTION("a file that is not a zip is rejected, and leaves nothing behind") {
        auto archive = temp.root/"not_a.zip";
        touch(archive);
        CHECK_THROWS(gpu::installer::install(archive.string()));
        CHECK_FALSE(fs::exists(folder));
        CHECK_FALSE(fs::exists(temp.root/"gpu.staging"));
    }
}
