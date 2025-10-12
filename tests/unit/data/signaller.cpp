#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>
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

        signaller::BoundSignaller bs3(10, nullptr);
        REQUIRE(bs3.get_id() == 10);
    }
}

TEST_CASE("BoundSignaller::get_id") {
    signaller::BoundSignaller bs0(0, nullptr);
    signaller::BoundSignaller bs5(5, nullptr);
    signaller::BoundSignaller bs100(100, nullptr);
    
    CHECK(bs0.get_id() == 0);
    CHECK(bs5.get_id() == 5);
    CHECK(bs100.get_id() == 100);
}

TEST_CASE_METHOD(fixture, "BoundSignaller::modified_external") {
    SECTION("calls StateManager::externally_modified") {
        auto bs0 = std::make_shared<signaller::BoundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::BoundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_externally_modified(0) == false);
        bs0->modified_external();
        CHECK(sm.is_externally_modified(0) == true);

        CHECK(sm.is_externally_modified(2) == false);
        bs2->modified_external();
        CHECK(sm.is_externally_modified(2) == true);
    }

    SECTION("multiple calls") {
        auto bs = std::make_shared<signaller::BoundSignaller>(1, &sm);
        sm.set_probe(1, bs);

        bs->modified_external();
        CHECK(sm.is_externally_modified(1) == true);
        
        sm.reset_to_false();
        CHECK(sm.is_externally_modified(1) == false);
        
        bs->modified_external();
        CHECK(sm.is_externally_modified(1) == true);
    }
}

TEST_CASE_METHOD(fixture, "BoundSignaller::modified_internal") {
    SECTION("calls StateManager::internally_modified") {
        auto bs0 = std::make_shared<signaller::BoundSignaller>(0, &sm);
        auto bs2 = std::make_shared<signaller::BoundSignaller>(2, &sm);
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_internally_modified(0) == false);
        bs0->modified_internal();
        CHECK(sm.is_internally_modified(0) == true);

        CHECK(sm.is_internally_modified(2) == false);
        bs2->modified_internal();
        CHECK(sm.is_internally_modified(2) == true);
    }

    SECTION("multiple calls") {
        auto bs = std::make_shared<signaller::BoundSignaller>(3, &sm);
        sm.set_probe(3, bs);

        bs->modified_internal();
        CHECK(sm.is_internally_modified(3) == true);
        
        sm.reset_to_false();
        CHECK(sm.is_internally_modified(3) == false);
        
        bs->modified_internal();
        CHECK(sm.is_internally_modified(3) == true);
    }
}

TEST_CASE("UnboundSignaller::UnboundSignaller") {
    SECTION("default constructor") {
        signaller::UnboundSignaller us;
    }
}

TEST_CASE_METHOD(fixture, "UnboundSignaller::modified_external") {
    SECTION("does not call StateManager") {
        auto bs0 = std::make_shared<signaller::UnboundSignaller>();
        auto bs2 = std::make_shared<signaller::UnboundSignaller>();
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_externally_modified(0) == false);
        bs0->modified_external();
        CHECK(sm.is_externally_modified(0) == false);

        CHECK(sm.is_externally_modified(2) == false);
        bs2->modified_external();
        CHECK(sm.is_externally_modified(2) == false);
    }

    SECTION("multiple calls have no effect") {
        auto bs = std::make_shared<signaller::UnboundSignaller>();
        sm.set_probe(1, bs);

        bs->modified_external();
        bs->modified_external();
        bs->modified_external();
        CHECK(sm.is_externally_modified(1) == false);
    }
}

TEST_CASE_METHOD(fixture, "UnboundSignaller::modified_internal") {
    SECTION("does not call StateManager") {
        auto bs0 = std::make_shared<signaller::UnboundSignaller>();
        auto bs2 = std::make_shared<signaller::UnboundSignaller>();
        sm.set_probe(0, bs0);
        sm.set_probe(2, bs2);

        CHECK(sm.is_internally_modified(0) == false);
        bs0->modified_internal();
        CHECK(sm.is_internally_modified(0) == false);

        CHECK(sm.is_internally_modified(2) == false);
        bs2->modified_internal();
        CHECK(sm.is_internally_modified(2) == false);
    }

    SECTION("multiple calls have no effect") {
        auto bs = std::make_shared<signaller::UnboundSignaller>();
        sm.set_probe(3, bs);

        bs->modified_internal();
        bs->modified_internal();
        bs->modified_internal();
        CHECK(sm.is_internally_modified(3) == false);
    }
}
