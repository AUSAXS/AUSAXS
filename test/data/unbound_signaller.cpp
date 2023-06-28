#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/state/UnboundSignaller.h>
#include <data/state/StateManager.h>

struct fixture {
    fixture() : sm(5) {
        sm.reset();
    }
    StateManager sm;
};

TEST_CASE_METHOD(fixture, "UnboundSignaller::external_change") {
    SECTION("calls StateManager::externally_modified") {
        auto bs0 = std::make_shared<signaller::UnboundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::UnboundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_externally_modified(0) == false);
        bs0->external_change();
        CHECK(sm.is_externally_modified(0) == false);

        CHECK(sm.is_externally_modified(2) == false);
        bs2->external_change();
        CHECK(sm.is_externally_modified(2) == false);
    }
}

TEST_CASE_METHOD(fixture, "UnboundSignaller::internal_change") {
    SECTION("calls StateManager::internally_modified") {
        auto bs0 = std::make_shared<signaller::UnboundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::UnboundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_internally_modified(0) == false);
        bs0->internal_change();
        CHECK(sm.is_internally_modified(0) == false);

        CHECK(sm.is_internally_modified(2) == false);
        bs2->internal_change();
        CHECK(sm.is_internally_modified(2) == false);
    }
}