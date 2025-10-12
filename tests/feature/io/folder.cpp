#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/Folder.h>

#include <filesystem>

using namespace ausaxs;

TEST_CASE("Folder::Folder") {
    SECTION("empty") {
        io::Folder folder("");
        CHECK(folder.path() == ".");
    }

    SECTION("simple") {
        io::Folder folder("test");
        CHECK(folder.path() == "test");
    }

    SECTION("simple with trailing slash") {
        io::Folder folder("tests/");
        CHECK(folder.path() == "tests");
    }

    SECTION("nested") {
        io::Folder folder("tests/folder");
        CHECK(folder.path() == "tests/folder");
    }

    SECTION("nested with trailing slash") {
        io::Folder folder("tests/folder/");
        CHECK(folder.path() == "tests/folder");
    }
}

TEST_CASE("Folder::exists") {
    SECTION("false") {
        io::Folder folder("fake");
        CHECK(folder.exists() == false);

        io::Folder folder2("fake/folder");
        CHECK(folder2.exists() == false);
    }

    SECTION("true") {
        io::Folder folder("tests");
        CHECK(folder.exists());

        io::Folder folder2("tests/files");
        CHECK(folder2.exists());
    }
}

TEST_CASE("Folder::create") {
    SECTION("simple") {
        std::string path = "___dummy";
        io::Folder folder(path);
        folder.create();
        CHECK(folder.exists());
        std::filesystem::remove(path);
    }

    SECTION("nested") {
        std::string path = "dummy/nested";
        io::Folder folder(path);
        folder.create();
        CHECK(folder.exists());
        std::filesystem::remove(path);
    }
}