#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/File.h>

#include <string>
#include <filesystem>

TEST_CASE("file_constructor") {
    SECTION("simple") {
        io::File file("test");
        CHECK(file.path() == "./test");
    }

    SECTION("with folder") {
        io::File file("test/file.txt");
        CHECK(file.directory().path() == "test");
        CHECK(file.stem() == "file");
        CHECK(file.extension() == ".txt");
        CHECK(file.path() == "test/file.txt");
    }

    SECTION("nested folders") {
        io::File file("test/folder/file.txt");
        CHECK(file.directory().path() == "test/folder");
        CHECK(file.stem() == "file");
        CHECK(file.extension() == ".txt");
        CHECK(file.path() == "test/folder/file.txt");
    }
}

TEST_CASE("file_exists") {
    SECTION("false") {
        io::File file("fake");
        CHECK(file.exists() == false);

        io::File file2("fake/file.txt");
        CHECK(file2.exists() == false);
    }

    SECTION("true") {
        io::File file("test/files/2epe.dat");
        CHECK(file.exists());
    }
}

TEST_CASE("file_append") {
    io::File file("test/file.txt");
    file.append("_2");
    CHECK(file.path() == "test/file_2.txt");
}

TEST_CASE("file_stem") {
    io::File file("test/file.txt");
    CHECK(file.stem() == "file");
}

TEST_CASE("file_split") {
    auto [dir, name, ext] = io::File::split("test/file.txt");
    CHECK(dir == "test");
    CHECK(name == "file");
    CHECK(ext == ".txt");
}

TEST_CASE("file_create") {
    std::string path = "temp/dummy.txt";
    io::File file(path);
    file.create();
    CHECK(file.exists());
    std::filesystem::remove(path);
}