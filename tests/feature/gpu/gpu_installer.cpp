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
    // Points the cache, and with it the install folder, at a fresh temporary folder. CTest runs each test case as its own process,
    // in parallel with the others, so every test case needs a folder of its own.
    struct TempCache {
        explicit TempCache(std::string_view name) : root("temp/tests/gpu_installer/" + std::string(name)), previous(settings::general::cache) {
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

        io::Folder root;
        io::Folder folder = root.path() + "/gpu";
        io::Folder staging = root.path() + "/gpu.staging";
        io::File record = root.path() + "/gpu.record";
        io::File toggle = root.path() + "/.autotoggle_gpu";
        std::string library{gpu::GPULoader::library_name()};
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
    TempCache temp("installed_library");

    SECTION("an incomplete installation is ignored") {
        CHECK(gpu::installer::installed_library().empty());
        CHECK(gpu::installer::installed_source().empty());

        // the library, but no marker: the installation did not finish
        temp.file(temp.library).create();
        CHECK(gpu::installer::installed_library().empty());
        CHECK(gpu::installer::installed_source().empty());

        // a marker, but no library
        temp.file(temp.library).remove();
        temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
        CHECK(gpu::installer::installed_library().empty());
        CHECK(gpu::installer::installed_device().empty());
    }

    SECTION("a finished installation is found") {
        temp.file(temp.library).create();
        temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
        CHECK(gpu::installer::installed_library().path() == temp.file(temp.library).path());
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
        CHECK(gpu::installer::installed_device() == "cuda: NVIDIA GeForce RTX 3050");

        // one that found no GPU, or predates choosing one, has no chosen device
        temp.file("installed.txt").create("https://example.org/backend.zip\n");
        CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
        CHECK(gpu::installer::installed_device().empty());
    }
}

TEST_CASE("installer::remove") {
    TempCache temp("remove");

    SECTION("without a record, only an empty folder is removed") {
        CHECK_FALSE(gpu::installer::remove());

        temp.folder.create();
        CHECK_FALSE(gpu::installer::remove());
        CHECK_FALSE(temp.folder.exists());

        temp.file(temp.library).create();
        CHECK_THROWS(gpu::installer::remove());
        CHECK(temp.file(temp.library).exists());
    }

    SECTION("an installation is removed, along with its record") {
        temp.file(temp.library).create();
        temp.file("hipSYCL/plugin.so").create();
        write_record(temp.record, {temp.library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(temp.folder.exists());
        CHECK_FALSE(temp.record.exists());
    }

    SECTION("files that were not installed are kept, along with the folders holding them") {
        temp.file(temp.library).create();
        temp.file("hipSYCL/plugin.so").create();
        temp.file("hipSYCL/mine.txt").create();
        temp.file("notes.txt").create();
        write_record(temp.record, {temp.library, "hipSYCL", "hipSYCL/plugin.so"});
        CHECK(gpu::installer::remove());
        CHECK_FALSE(temp.file(temp.library).exists());
        CHECK_FALSE(temp.file("hipSYCL/plugin.so").exists());
        CHECK(temp.file("hipSYCL/mine.txt").exists());
        CHECK(temp.file("notes.txt").exists());
        CHECK_FALSE(temp.record.exists());
    }
}

TEST_CASE("installer::enabled") {
    TempCache temp("enabled");

    SECTION("cannot be turned on without an installation, and a leftover toggle does nothing") {
        temp.toggle.create();
        CHECK_FALSE(gpu::installer::enabled());
        CHECK_THROWS(gpu::installer::set_enabled(true));
    }

    SECTION("turned on and off, leaving the installation alone") {
        temp.file(temp.library).create();
        temp.file("installed.txt").create("https://example.org/backend.zip\n\n");
        write_record(temp.record, {temp.library, "installed.txt"});
        CHECK_FALSE(gpu::installer::enabled());
        CHECK_THROWS(gpu::installer::set_enabled(true)); // no GPU has been chosen
        CHECK_FALSE(gpu::installer::enabled());

        temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
        gpu::installer::set_enabled(true);
        CHECK(gpu::installer::enabled());

        gpu::installer::set_enabled(false);
        CHECK_FALSE(gpu::installer::enabled());
        CHECK(temp.file(temp.library).exists());

        gpu::installer::set_enabled(true);
        CHECK(gpu::installer::remove());
        CHECK_FALSE(temp.toggle.exists());
    }
}

TEST_CASE("installer::select_device") {
    TempCache temp("select_device");
    CHECK_FALSE(gpu::installer::select_device()); // nothing installed

    // an installation that cannot be loaded is left as it was
    temp.file(temp.library).create();
    temp.file("installed.txt").create("https://example.org/backend.zip\ncuda: NVIDIA GeForce RTX 3050\n");
    write_record(temp.record, {temp.library, "installed.txt"});
    gpu::installer::set_enabled(true);
    CHECK_FALSE(gpu::installer::select_device());
    CHECK(temp.file(temp.library).exists());
    CHECK(gpu::installer::installed_source() == "https://example.org/backend.zip");
    CHECK(gpu::installer::installed_device() == "cuda: NVIDIA GeForce RTX 3050");
    CHECK(gpu::installer::enabled());
}

TEST_CASE("installer::install") {
    TempCache temp("install");

    SECTION("a file that is not a zip is rejected, and leaves nothing behind") {
        io::File archive = temp.root.path() + "/not_a.zip";
        archive.create();
        CHECK_THROWS(gpu::installer::install(archive.path()));
        CHECK_FALSE(temp.folder.exists());
        CHECK_FALSE(temp.staging.exists());
        CHECK_FALSE(temp.record.exists());
    }

    SECTION("an existing path that is not from an earlier installation is never replaced") {
        auto zip = make_backend_zip(temp.root.path() + "/downloads/backend.zip");
        temp.file("hipSYCL/mine.txt").create();
        CHECK_THROWS(gpu::installer::install(zip.path()));
        CHECK(temp.file("hipSYCL/mine.txt").exists());
        CHECK_FALSE(temp.file(temp.library).exists());
        CHECK_FALSE(temp.staging.exists());
    }

    SECTION("an earlier installation is replaced, and a backend that fails to load is removed again, keeping other files") {
        auto zip = make_backend_zip(temp.root.path() + "/downloads/backend.zip");
        temp.file(temp.library).create(); // also in the new backend, but recorded, so it may be replaced
        temp.file("old.so").create();
        temp.file("notes.txt").create();
        write_record(temp.record, {temp.library, "old.so"});
        CHECK_THROWS(gpu::installer::install(zip.path())); // the fake library fails the load test
        CHECK_FALSE(temp.file("old.so").exists());
        CHECK_FALSE(temp.file(temp.library).exists());
        CHECK_FALSE(temp.dir("hipSYCL").exists());
        CHECK(temp.file("notes.txt").exists());
        CHECK_FALSE(temp.record.exists());
        CHECK_FALSE(temp.staging.exists());
    }
}
