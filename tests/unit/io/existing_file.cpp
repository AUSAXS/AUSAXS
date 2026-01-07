#include <catch2/catch_test_macros.hpp>

#include <io/ExistingFile.h>

using namespace ausaxs;

TEST_CASE("ExistingFile::ExistingFile") {
    SECTION("non-existing") {
        CHECK_THROWS(io::ExistingFile("fake"));
        CHECK_THROWS(io::ExistingFile("fake/file.txt"));
    }

    SECTION("existing") {
        io::ExistingFile file("tests/files/2epe.dat");
        CHECK(file.path() == "tests/files/2epe.dat");
    }
}

TEST_CASE("ExistingFile::empty") {
    io::ExistingFile file;
    CHECK(file.empty());
}
