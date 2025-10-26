#include <catch2/catch_test_macros.hpp>

#include <grid/detail/GridMember.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <form_factor/FormFactorType.h>

using namespace ausaxs;
using namespace ausaxs::grid;
using namespace ausaxs::data;

TEST_CASE("GridMember<AtomFF>::constructor") {
    SECTION("default") {
        GridMember<AtomFF> member;
        REQUIRE_FALSE(member.is_expanded());
    }

    SECTION("with atom and location") {
        AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
        Vector3<int> loc(5, 6, 7);
        GridMember<AtomFF> member(atom, loc);
        
        REQUIRE(member.get_bin_loc() == loc);
        REQUIRE(member.get_absolute_loc() == atom.coordinates());
        REQUIRE_FALSE(member.is_expanded());
    }

    SECTION("copy constructor") {
        AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
        Vector3<int> loc(5, 6, 7);
        GridMember<AtomFF> member1(atom, loc);
        member1.set_expanded(true);
        
        GridMember<AtomFF> member2(member1);
        REQUIRE(member2.get_bin_loc() == member1.get_bin_loc());
        REQUIRE(member2.get_absolute_loc() == member1.get_absolute_loc());
        REQUIRE(member2.is_expanded() == member1.is_expanded());
    }
}

TEST_CASE("GridMember<Water>::constructor") {
    SECTION("default") {
        GridMember<Water> member;
        REQUIRE_FALSE(member.is_expanded());
    }

    SECTION("with water and location") {
        Water water({1.0, 2.0, 3.0});
        Vector3<int> loc(5, 6, 7);
        GridMember<Water> member(water, loc);
        
        REQUIRE(member.get_bin_loc() == loc);
        REQUIRE(member.get_absolute_loc() == water.coordinates());
        REQUIRE_FALSE(member.is_expanded());
    }
}

TEST_CASE("GridMember::get_bin_loc") {
    AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    Vector3<int> loc(5, 6, 7);
    GridMember<AtomFF> member(atom, loc);

    SECTION("const accessor") {
        const GridMember<AtomFF>& const_member = member;
        REQUIRE(const_member.get_bin_loc() == loc);
    }

    SECTION("non-const accessor") {
        member.get_bin_loc().x() = 10;
        REQUIRE(member.get_bin_loc().x() == 10);
    }
}

TEST_CASE("GridMember::get_absolute_loc") {
    AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    Vector3<int> loc(5, 6, 7);
    GridMember<AtomFF> member(atom, loc);

    SECTION("const accessor") {
        const GridMember<AtomFF>& const_member = member;
        REQUIRE(const_member.get_absolute_loc() == Vector3<double>(1.0, 2.0, 3.0));
    }

    SECTION("non-const accessor") {
        member.get_absolute_loc().x() = 10.0;
        REQUIRE(member.get_absolute_loc().x() == 10.0);
    }
}

TEST_CASE("GridMember::expansion state") {
    AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    GridMember<AtomFF> member(atom, Vector3<int>(5, 6, 7));

    SECTION("initially not expanded") {
        REQUIRE_FALSE(member.is_expanded());
    }

    SECTION("set expanded") {
        member.set_expanded(true);
        REQUIRE(member.is_expanded());
        member.set_expanded(false);
        REQUIRE_FALSE(member.is_expanded());
    }
}

TEST_CASE("GridMember::get_atom") {
    AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    GridMember<AtomFF> member(atom, Vector3<int>(5, 6, 7));

    SECTION("const accessor") {
        const GridMember<AtomFF>& const_member = member;
        REQUIRE(const_member.get_atom() == atom);
    }

    SECTION("non-const accessor") {
        member.get_atom().coordinates().x() = 10.0;
        REQUIRE(member.get_atom().coordinates().x() == 10.0);
    }
}

TEST_CASE("GridMember::get_atom_type") {
    SECTION("AtomFF") {
        AtomFF atom({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
        GridMember<AtomFF> member(atom, Vector3<int>(5, 6, 7));
        REQUIRE(member.get_atom_type() == form_factor::form_factor_t::C);
    }

    SECTION("Water") {
        Water water({1.0, 2.0, 3.0});
        GridMember<Water> member(water, Vector3<int>(5, 6, 7));
        REQUIRE(member.get_atom_type() == form_factor::form_factor_t::OH);
    }
}

TEST_CASE("GridMember::operator==") {
    AtomFF atom1({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    AtomFF atom2({4.0, 5.0, 6.0}, form_factor::form_factor_t::N);
    Vector3<int> loc(5, 6, 7);

    SECTION("compare with GridMember") {
        GridMember<AtomFF> member1(atom1, loc);
        GridMember<AtomFF> member2(atom1, loc);
        GridMember<AtomFF> member3(atom2, loc);
        
        REQUIRE(member1 == member2);
        REQUIRE_FALSE(member1 == member3);
    }

    SECTION("compare with atom") {
        GridMember<AtomFF> member(atom1, loc);
        REQUIRE(member == atom1);
        REQUIRE_FALSE(member == atom2);
    }
}

TEST_CASE("GridMember::operator=") {
    AtomFF atom1({1.0, 2.0, 3.0}, form_factor::form_factor_t::C);
    AtomFF atom2({4.0, 5.0, 6.0}, form_factor::form_factor_t::N);

    SECTION("copy assignment") {
        GridMember<AtomFF> member1(atom1, Vector3<int>(5, 6, 7));
        member1.set_expanded(true);
        GridMember<AtomFF> member2(atom2, Vector3<int>(1, 2, 3));
        
        member2 = member1;
        REQUIRE(member2.get_bin_loc() == member1.get_bin_loc());
        REQUIRE(member2.get_atom() == member1.get_atom());
        REQUIRE(member2.is_expanded() == member1.is_expanded());
    }

    SECTION("move assignment") {
        GridMember<AtomFF> member1(atom1, Vector3<int>(5, 6, 7));
        member1.set_expanded(true);
        GridMember<AtomFF> member2(atom2, Vector3<int>(1, 2, 3));
        
        Vector3<int> original_loc = member1.get_bin_loc();
        member2 = std::move(member1);
        REQUIRE(member2.get_bin_loc() == original_loc);
        REQUIRE(member2.is_expanded());
    }
}
