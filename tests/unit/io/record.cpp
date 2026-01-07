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
