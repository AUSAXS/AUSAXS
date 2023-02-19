#include <catch2/catch_test_macros.hpp>

#include <data/Record.h>
#include <data/Terminate.h>
#include <data/Footer.h>
#include <data/Header.h>

TEST_CASE("record") {
    CHECK(Record::get_type("ATOM") == Record::RecordType::ATOM);
    CHECK(Record::get_type("HETATM") == Record::RecordType::ATOM);
    CHECK(Record::get_type("TER") == Record::RecordType::TERMINATE);
    CHECK(Record::get_type("HEADER") == Record::RecordType::HEADER);
    CHECK(Record::get_type("CONECT") == Record::RecordType::FOOTER);
}


TEST_CASE("terminate") {
    Terminate t1(1, "LYS", "A", 2, "B");
    CHECK(t1.serial == 1);
    CHECK(t1.resName == "LYS");
    CHECK(t1.chainID == "A");
    CHECK(t1.resSeq == 2);
    CHECK(t1.iCode == "B");
    CHECK(t1.get_type() == Record::RecordType::TERMINATE);

    Terminate t2(3, "ARG", "C", 4, "D");
    std::string t2_str = t2.as_pdb();

    t1.parse_pdb(t2_str);
    CHECK(t1.serial == 3);
    CHECK(t1.resName == "ARG");
    CHECK(t1.chainID == "C");
    CHECK(t1.resSeq == 4);
    CHECK(t1.iCode == "D");

    t1.set_serial(5);
    CHECK(t1.serial == 5);
}

TEST_CASE("footer") {
    Footer f1;
    CHECK(f1.get_type() == Record::RecordType::FOOTER);

    f1.add("CONECT test1");
    CHECK(f1.get() == "CONECT test1\n");

    f1.add("CONECT test2");
    CHECK(f1.get() == "CONECT test1\nCONECT test2\n");

    f1.add("END    test3");
    CHECK(f1.size() == 3);
    f1.remove("CONECT");
    CHECK(f1.size() == 1);
    CHECK(f1.get() == "END    test3\n");
}

TEST_CASE("header") {
    Header f1;
    CHECK(f1.get_type() == Record::RecordType::HEADER);

    f1.add("HEADER test1");
    CHECK(f1.get() == "HEADER test1\n");

    f1.add("TITLE  test2");
    CHECK(f1.get() == "HEADER test1\nTITLE  test2\n");

    f1.add("TITLE  test3");
    CHECK(f1.size() == 3);
    f1.remove("TITLE");
    CHECK(f1.size() == 1);
    CHECK(f1.get() == "HEADER test1\n");
}