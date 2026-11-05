// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>
#include <settings/GeneralSettings.h>

#include <miniz.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace ausaxs;
namespace fs = std::filesystem;

namespace {
    // Points the cache, and with it the install folder, at a fresh temporary folder.
    struct TempCache {
        TempCache() : previous(settings::general::cache) {
            fs::remove_all(root);
            settings::general::cache = root.string() + "/";
        }
        ~TempCache() {
            settings::general::cache = previous;
            fs::remove_all(root);
        }
        fs::path root = "temp/tests/gpu_installer";
        fs::path folder = root/"gpu";
        fs::path record = root/"gpu.record";
        std::string previous;
    };

    void touch(const fs::path& file) {
        fs::create_directories(file.parent_path());
        std::ofstream(file) << "content";
    }

    void write_record(const fs::path& record, const std::vector<std::string>& entries) {
        std::ofstream out(record);
        for (const auto& entry : entries) {out << entry << "\n";}
    }

    // A zip shaped like a real backend package. The library is not a real one, so it always fails the load test.
    fs::path make_backend_zip(const fs::path& zip) {
        fs::create_directories(zip.parent_path());
        mz_zip_archive archive{};
        REQUIRE(mz_zip_writer_init_file(&archive, zip.string().c_str(), 0));
        std::string content = "not a library";
        for (const std::string& name : {std::string(gpu::GPULoader::library_name()), std::string("hipSYCL/plugin.so")}) {
            REQUIRE(mz_zip_writer_add_mem(&archive, name.c_str(), content.data(), content.size(), MZ_DEFAULT_COMPRESSION));
        }
        REQUIRE(mz_zip_writer_finalize_archive(&archive));
        REQUIRE(mz_zip_writer_end(&archive));
        return zip;
    }
}

TEST_CASE("installer::folder follows the cache") {
    TempCache temp;
    CHECK(gpu::installer::folder().path() == temp.folder.string());

    settings::general::cache = (temp.root/"other").string();
    CHECK(gpu::installer::folder().path() == (temp.root/"other"/"gpu").string());
}

TEST_CASE("installer::asset_url") {
    if (gpu::installer::asset_name().empty()) {SKIP("no prebuilt backend for this platform");}
    auto url = gpu::installer::asset_url();
    CHECK(url.starts_with("https://"));
    CHECK(url.ends_with("/latest/download/" + gpu::installer::asset_name()));
}

TEST_CASE("installer::installed_library") {
    TempCache temp;
    auto library = temp.folder/gpu::GPULoader::library_name();

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
        std::ofstream(temp.folder/"installed.txt") << "https://example.org/backend.zip\n";
        CHECK(gpu::installer::installed_library().path() == library.string());
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
    }
}

TEST_CASE("installer::remove") {
    TempCache temp;
    std::string library(gpu::GPULoader::library_name());

    SECTION("nothing to remove") {
        CHECK_FALSE(gpu::installer::remove());
    }

    SECTION("an installation is removed, along with its record") {
        touch(temp.folder/library);
        touch(temp.folder/"hipSYCL"/"plugin.so");
        write_record(temp.record, {library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(fs::exists(temp.folder));
        CHECK_FALSE(fs::exists(temp.record));
    }

    SECTION("files that were not installed are kept, along with the folders holding them") {
        touch(temp.folder/library);
        touch(temp.folder/"hipSYCL"/"plugin.so");
        touch(temp.folder/"hipSYCL"/"mine.txt");
        touch(temp.folder/"notes.txt");
        write_record(temp.record, {library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(fs::exists(temp.folder/library));
        CHECK_FALSE(fs::exists(temp.folder/"hipSYCL"/"plugin.so"));
        CHECK(fs::exists(temp.folder/"hipSYCL"/"mine.txt"));
        CHECK(fs::exists(temp.folder/"notes.txt"));
    }

    SECTION("a folder without a record is left alone") {
        touch(temp.folder/library);
        CHECK_THROWS(gpu::installer::remove());
        CHECK(fs::exists(temp.folder/library));
    }

    SECTION("an empty folder without a record is tidied up") {
        fs::create_directories(temp.folder);
        CHECK_FALSE(gpu::installer::remove());
        CHECK_FALSE(fs::exists(temp.folder));
    }
}

TEST_CASE("installer::install") {
    TempCache temp;
    auto zip = make_backend_zip(temp.root/"downloads"/"backend.zip");

    SECTION("a file that is not a zip is rejected, and leaves nothing behind") {
        auto archive = temp.root/"not_a.zip";
        touch(archive);
        CHECK_THROWS(gpu::installer::install(archive.string()));
        CHECK_FALSE(fs::exists(temp.folder));
        CHECK_FALSE(fs::exists(temp.root/"gpu.staging"));
        CHECK_FALSE(fs::exists(temp.record));
    }

    SECTION("an existing path that is not from an earlier installation is never replaced") {
        touch(temp.folder/"hipSYCL"/"mine.txt");
        CHECK_THROWS(gpu::installer::install(zip.string()));
        CHECK(fs::exists(temp.folder/"hipSYCL"/"mine.txt"));
        CHECK_FALSE(fs::exists(temp.folder/gpu::GPULoader::library_name()));
        CHECK_FALSE(fs::exists(temp.root/"gpu.staging"));
    }

    SECTION("a backend that fails to load is removed again, and other files are kept") {
        touch(temp.folder/"notes.txt");
        CHECK_THROWS(gpu::installer::install(zip.string()));
        CHECK(fs::exists(temp.folder/"notes.txt"));
        CHECK_FALSE(fs::exists(temp.folder/gpu::GPULoader::library_name()));
        CHECK_FALSE(fs::exists(temp.folder/"hipSYCL"));
        CHECK_FALSE(fs::exists(temp.record));
        CHECK_FALSE(fs::exists(temp.root/"gpu.staging"));
    }

    SECTION("an earlier installation is replaced") {
        std::string library(gpu::GPULoader::library_name());
        touch(temp.folder/library);
        touch(temp.folder/"old.so");
        write_record(temp.record, {library, "old.so"});
        CHECK_THROWS(gpu::installer::install(zip.string())); // the fake library fails the load test
        CHECK_FALSE(fs::exists(temp.folder/"old.so"));
        CHECK_FALSE(fs::exists(temp.folder));
    }
}
