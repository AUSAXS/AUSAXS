// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <catch2/catch_test_macros.hpp>

#include <gpu/GPUInstaller.h>
#include <gpu/GPULoader.h>
#include <io/File.h>
#include <io/Folder.h>
#include <settings/GeneralSettings.h>

#include <string>
#include <vector>

#include <miniz.h>

using namespace ausaxs;

namespace {
    // Points the cache, and with it the install folder, at a fresh temporary folder.
    struct TempCache {
        TempCache() : previous(settings::general::cache) {
            root.remove();
            settings::general::cache = root.path() + "/";
        }
        ~TempCache() {
            settings::general::cache = previous;
            root.remove();
        }

        // a file or folder inside the install folder
        io::File file(std::string_view relative) const {return folder.path() + "/" + std::string(relative);}
        io::Folder dir(std::string_view relative) const {return folder.path() + "/" + std::string(relative);}

        io::Folder root = "temp/tests/gpu_installer";
        io::Folder folder = root.path() + "/gpu";
        io::Folder staging = root.path() + "/gpu.staging";
        io::File record = root.path() + "/gpu.record";
        std::string previous;
    };

    void write_record(const io::File& record, const std::vector<std::string>& entries) {
        std::string contents;
        for (const auto& entry : entries) {contents += entry + "\n";}
        record.create(contents);
    }

    // A zip shaped like a real backend package. The library is not a real one, so it always fails the load test.
    io::File make_backend_zip(const io::File& zip) {
        zip.directory().create();
        mz_zip_archive archive{};
        REQUIRE(mz_zip_writer_init_file(&archive, zip.path().c_str(), 0));
        std::string content = "not a library";
        for (const std::string& name : {std::string(gpu::GPULoader::library_name()), std::string("hipSYCL/plugin.so")}) {
            REQUIRE(mz_zip_writer_add_mem(&archive, name.c_str(), content.data(), content.size(), MZ_DEFAULT_COMPRESSION));
        }
        REQUIRE(mz_zip_writer_finalize_archive(&archive));
        REQUIRE(mz_zip_writer_end(&archive));
        return zip;
    }
}

TEST_CASE("installer::installed_library") {
    TempCache temp;
    auto library = temp.file(gpu::GPULoader::library_name());

    SECTION("nothing installed") {
        CHECK(gpu::installer::installed_library().empty());
        CHECK(gpu::installer::installed_source().empty());
    }

    SECTION("an unfinished installation is ignored") {
        library.create();
        CHECK(gpu::installer::installed_library().empty());
    }

    SECTION("a finished installation is found") {
        library.create();
        temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
        CHECK(gpu::installer::installed_library().path() == library.path());
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
        CHECK(gpu::installer::installed_device() == "cuda: NVIDIA GeForce RTX 3050");
    }

    SECTION("an installation that found no GPU, or predates choosing one, has no chosen device") {
        library.create();
        temp.file("installed.txt").create("https://example.org/backend.zip\n");
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
        CHECK(gpu::installer::installed_device().empty());
    }

    SECTION("an unfinished installation has no chosen device") {
        temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
        CHECK(gpu::installer::installed_device().empty());
    }
}

TEST_CASE("installer::remove") {
    TempCache temp;
    std::string library(gpu::GPULoader::library_name());

    SECTION("nothing to remove") {
        CHECK_FALSE(gpu::installer::remove());
    }

    SECTION("an installation is removed, along with its record") {
        temp.file(library).create();
        temp.file("hipSYCL/plugin.so").create();
        write_record(temp.record, {library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(temp.folder.exists());
        CHECK_FALSE(temp.record.exists());
    }

    SECTION("files that were not installed are kept, along with the folders holding them") {
        temp.file(library).create();
        temp.file("hipSYCL/plugin.so").create();
        temp.file("hipSYCL/mine.txt").create();
        temp.file("notes.txt").create();
        write_record(temp.record, {library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(temp.file(library).exists());
        CHECK_FALSE(temp.file("hipSYCL/plugin.so").exists());
        CHECK(temp.file("hipSYCL/mine.txt").exists());
        CHECK(temp.file("notes.txt").exists());
    }

    SECTION("a folder without a record is left alone") {
        temp.file(library).create();
        CHECK_THROWS(gpu::installer::remove());
        CHECK(temp.file(library).exists());
    }

    SECTION("an empty folder without a record is tidied up") {
        temp.folder.create();
        CHECK_FALSE(gpu::installer::remove());
        CHECK_FALSE(temp.folder.exists());
    }
}

TEST_CASE("installer::install") {
    TempCache temp;
    auto zip = make_backend_zip(temp.root.path() + "/downloads/backend.zip");

    SECTION("a file that is not a zip is rejected, and leaves nothing behind") {
        io::File archive = temp.root.path() + "/not_a.zip";
        archive.create();
        CHECK_THROWS(gpu::installer::install(archive.path()));
        CHECK_FALSE(temp.folder.exists());
        CHECK_FALSE(temp.staging.exists());
        CHECK_FALSE(temp.record.exists());
    }

    SECTION("an existing path that is not from an earlier installation is never replaced") {
        temp.file("hipSYCL/mine.txt").create();
        CHECK_THROWS(gpu::installer::install(zip.path()));
        CHECK(temp.file("hipSYCL/mine.txt").exists());
        CHECK_FALSE(temp.file(gpu::GPULoader::library_name()).exists());
        CHECK_FALSE(temp.staging.exists());
    }

    SECTION("a backend that fails to load is removed again, and other files are kept") {
        temp.file("notes.txt").create();
        CHECK_THROWS(gpu::installer::install(zip.path()));
        CHECK(temp.file("notes.txt").exists());
        CHECK_FALSE(temp.file(gpu::GPULoader::library_name()).exists());
        CHECK_FALSE(temp.dir("hipSYCL").exists());
        CHECK_FALSE(temp.record.exists());
        CHECK_FALSE(temp.staging.exists());
    }

    SECTION("an earlier installation is replaced") {
        std::string library(gpu::GPULoader::library_name());
        temp.file(library).create();
        temp.file("old.so").create();
        write_record(temp.record, {library, "old.so"});
        CHECK_THROWS(gpu::installer::install(zip.path())); // the fake library fails the load test
        CHECK_FALSE(temp.file("old.so").exists());
        CHECK_FALSE(temp.folder.exists());
    }
}
