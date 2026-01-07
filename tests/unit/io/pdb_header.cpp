#include <catch2/catch_test_macros.hpp>

#include <io/pdb/Header.h>

using namespace ausaxs;
using namespace io::pdb;

TEST_CASE("Header::get_type") {
    Header header;
    REQUIRE(header.get_type() == RecordType::HEADER);
}

TEST_CASE("Header::parse_pdb") {
    Header header;
    header.parse_pdb("KEYWDS    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(header.as_pdb() == "KEYWDS    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n");
}

TEST_CASE("Header::as_pdb") {
    Header header;
    REQUIRE(header.get() == header.as_pdb());
}

TEST_CASE("Header::add") {
    Header header;
    header.add("HELIX    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(header.get() == "HELIX    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n");
}

TEST_CASE("Header::remove") {
    Header header;
    header.add("JRNL    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    header.add("AUTHOR    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2");
    header.add("JRNL    3 3 3 3 3 3 3 3 3 3 3 3 3 3 3");
    header.remove("JRNL");
    REQUIRE(header.get() == "AUTHOR    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2\n");
}

TEST_CASE("Header::get") {
    Header h1;
    h1.add("HEADER test1");
    CHECK(h1.get() == "HEADER test1\n");

    h1.add("TITLE  test2");
    CHECK(h1.get() == "HEADER test1\nTITLE  test2\n");
}

TEST_CASE("Header::size") {
    Header header;
    header.add("SEQRES    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(header.size() == 1);

    header.add("SSBOND    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2");
    REQUIRE(header.size() == 2);
}
