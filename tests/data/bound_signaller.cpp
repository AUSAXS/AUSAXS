#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/state/BoundSignaller.h>
#include <data/state/StateManager.h>

using namespace ausaxs;

struct fixture {
    fixture() : sm(5) {
        sm.reset_to_false();
    }
    state::StateManager sm;
};

TEST_CASE("BoundSignaller::BoundSignaller") {
    SECTION("uint, StateManager*") {
        signaller::BoundSignaller bs(0, nullptr);
        REQUIRE(bs.get_id() == 0);

        signaller::BoundSignaller bs2(1, nullptr);
        REQUIRE(bs2.get_id() == 1);
    }
}

TEST_CASE_METHOD(fixture, "BoundSignaller::external_change") {
    SECTION("calls StateManager::externally_modified") {
        auto bs0 = std::make_shared<signaller::BoundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::BoundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_externally_modified(0) == false);
        bs0->external_change();
        CHECK(sm.is_externally_modified(0) == true);

        CHECK(sm.is_externally_modified(2) == false);
        bs2->external_change();
        CHECK(sm.is_externally_modified(2) == true);
    }
}

TEST_CASE_METHOD(fixture, "BoundSignaller::internal_change") {
    SECTION("calls StateManager::internally_modified") {
        auto bs0 = std::make_shared<signaller::BoundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::BoundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_internally_modified(0) == false);
        bs0->internal_change();
        CHECK(sm.is_internally_modified(0) == true);

        CHECK(sm.is_internally_modified(2) == false);
        bs2->internal_change();
        CHECK(sm.is_internally_modified(2) == true);
    }
}