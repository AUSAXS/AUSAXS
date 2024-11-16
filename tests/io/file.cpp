#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/File.h>

#include <string>
#include <filesystem>
#include <fstream>

using namespace ausaxs;

TEST_CASE("File::File") {
    SECTION("simple") {
        io::File file("test");
        CHECK(file.path() == "./test");
    }

    SECTION("with folder") {
        io::File file("tests/file.txt");
        CHECK(file.directory().path() == "tests");
        CHECK(file.stem() == "file");
        CHECK(file.extension() == ".txt");
        CHECK(file.path() == "tests/file.txt");
    }

    SECTION("nested folders") {
        io::File file("tests/folder/file.txt");
        CHECK(file.directory().path() == "tests/folder");
        CHECK(file.stem() == "file");
        CHECK(file.extension() == ".txt");
        CHECK(file.path() == "tests/folder/file.txt");
    }
}

TEST_CASE("File::exists") {
    SECTION("false") {
        io::File file("fake");
        CHECK(file.exists() == false);

        io::File file2("fake/file.txt");
        CHECK(file2.exists() == false);
    }

    SECTION("true") {
        io::File file("tests/files/2epe.dat");
        CHECK(file.exists());
    }
}

TEST_CASE("File::append") {
    io::File file("tests/file.txt");
    file.append("_2");
    CHECK(file.path() == "tests/file_2.txt");
}

TEST_CASE("File::stem") {
    io::File file("tests/file.txt");
    CHECK(file.stem() == "file");
}

TEST_CASE("File::split") {
    auto [dir, name, ext] = io::File::split("tests/file.txt");
    CHECK(dir == "tests");
    CHECK(name == "file");
    CHECK(ext == ".txt");
}

TEST_CASE("File::create") {
    SECTION("empty") {
        std::string path = "temp/dummy.txt";
        io::File file(path);
        file.create();
        CHECK(file.exists());
        file.remove();
    }

    SECTION("with contents") {
        std::string path = "temp/dummy.txt";
        io::File file(path);
        file.create("test");
        CHECK(file.exists());
        std::ifstream f(path);
        std::string contents;
        f >> contents;
        CHECK(contents == "test");
        file.remove();
    }
}

TEST_CASE("File::remove") {
    std::string path = "temp/dummy.txt";
    io::File file(path);
    file.create();
    CHECK(file.exists() == true);
    file.remove();
    CHECK(file.exists() == false);
}