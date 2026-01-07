#include <catch2/catch_test_macros.hpp>

#include <io/pdb/Terminate.h>

using namespace ausaxs;
using namespace io::pdb;

TEST_CASE("Terminate::Terminate") {
    SECTION("int, string&, string&, int, string&") {
        Terminate t1(1, "LYS", 'A', 2, "B");
        CHECK(t1.serial == 1);
        CHECK(t1.resName == "LYS");
        CHECK(t1.chainID == 'A');
        CHECK(t1.resSeq == 2);
        CHECK(t1.iCode == "B");
        CHECK(t1.get_type() == RecordType::TERMINATE);
    }
}

TEST_CASE("Terminate::get_type") {
    Terminate terminate;
    REQUIRE(terminate.get_type() == RecordType::TERMINATE);
}

TEST_CASE("Terminate::parse_pdb") {
    Terminate terminate;
    terminate.parse_pdb("TER       1      LYS A   2");
    CHECK(terminate.serial == 1);
    CHECK(terminate.resName == "LYS");
    CHECK(terminate.chainID == 'A');
    CHECK(terminate.resSeq == 2);
    CHECK(terminate.iCode == " ");
}

TEST_CASE("Terminate::as_pdb") {
    Terminate terminate;
    REQUIRE(terminate.get() == terminate.as_pdb());
}

TEST_CASE("Terminate::get") {
    Terminate terminate;
    std::string res = "TER       1      LYS A   2";
    terminate.parse_pdb(res);
    REQUIRE(terminate.get().substr(0, res.size()) == res);
}

TEST_CASE("Terminate::set_serial") {
    Terminate terminate;
    terminate.set_serial(5);
    CHECK(terminate.serial == 5);

    terminate.parse_pdb("TER       1      LYS A   2");
    terminate.set_serial(2);
    CHECK(terminate.serial == 2);
}
