#include <catch2/catch_test_macros.hpp>

#include <io/pdb/Record.h>

using namespace ausaxs;
using namespace io::pdb;

TEST_CASE("Record::get_type") {
    CHECK(Record::get_type("ATOM") == RecordType::ATOM);
    CHECK(Record::get_type("HETATM") == RecordType::ATOM);
    CHECK(Record::get_type("TER") == RecordType::TERMINATE);
    CHECK(Record::get_type("HEADER") == RecordType::HEADER);
    CHECK(Record::get_type("CONECT") == RecordType::FOOTER);
}

TEST_CASE("Record::find_type") {
    CHECK(Record::find_type("ATOM  ") == RecordType::ATOM);
    CHECK(Record::find_type("ANISOU") == RecordType::NOTYPE);
    CHECK_FALSE(Record::find_type("FOOBAR").has_value());
    CHECK_THROWS(Record::get_type("FOOBAR"));
}
