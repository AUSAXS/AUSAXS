#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <io/ExistingFile.h>

TEST_CASE("existing_file_constructor") {
    SECTION("non-existing") {
        CHECK_THROWS(io::ExistingFile("fake"));
        CHECK_THROWS(io::ExistingFile("fake/file.txt"));
    }

    SECTION("existing") {
        io::ExistingFile file("test/files/2epe.dat");
        CHECK(file.path() == "test/files/2epe.dat");
    }
}