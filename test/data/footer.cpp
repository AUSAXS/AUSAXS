#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Footer.h>
#include <settings/All.h>

TEST_CASE("Footer::get_type") {
    Footer footer;
    REQUIRE(footer.get_type() == RecordType::FOOTER);
}

TEST_CASE("Footer::parse_pdb") {
    Footer footer;
    footer.parse_pdb("CONECT    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(footer.get() == "CONECT    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n");
}

TEST_CASE("Footer::as_pdb") {
    Footer footer;
    REQUIRE(footer.as_pdb() == footer.get());
}

TEST_CASE("Footer::add") {
    Footer footer;
    footer.add("CONECT    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(footer.get() == "CONECT    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n");
}

TEST_CASE("Footer::remove") {
    Footer footer;
    footer.add("ENDMDL    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    footer.add("CONECT    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2");
    footer.add("ENDMDL    3 3 3 3 3 3 3 3 3 3 3 3 3 3 3");
    footer.remove("ENDMDL");
    REQUIRE(footer.get() == "CONECT    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2\n");
}

TEST_CASE("Footer::get") {
    Footer f1;
    f1.add("CONECT test1");
    CHECK(f1.get() == "CONECT test1\n");

    f1.add("CONECT test2");
    CHECK(f1.get() == "CONECT test1\nCONECT test2\n");
}

TEST_CASE("Footer::size") {
    Footer footer;
    footer.add("CONECT    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1");
    REQUIRE(footer.size() == 1);

    footer.add("ENDMDL    2 2 2 2 2 2 2 2 2 2 2 2 2 2 2");
    REQUIRE(footer.size() == 2);
}