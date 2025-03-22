#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/state/StateManager.h>
#include <data/state/Signaller.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>

#include <memory>

using namespace ausaxs;
using namespace state;

struct fixture {
    fixture() : manager(5) {
        manager.reset_to_false();
    }
    StateManager manager;
};

TEST_CASE("StateManager::StateManager") {
    SECTION("uint") {
        unsigned int size = 5;
        StateManager manager(size);
        REQUIRE(manager.get_internally_modified_bodies().size() == size);
        CHECK(manager.get_probes().size() == size);
        CHECK(manager.get_externally_modified_bodies() == std::vector{true, true, true, true, true});
        CHECK(manager.get_internally_modified_bodies() == std::vector{true, true, true, true, true});
        CHECK(manager.is_modified_hydration() == true);
    }
}

TEST_CASE_METHOD(fixture, "StateManager::externally_modified_all") {
    manager.externally_modified_all();
    CHECK(manager.get_externally_modified_bodies() == std::vector{true, true, true, true, true});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.is_modified_hydration() == false);
}

TEST_CASE_METHOD(fixture, "StateManager::internally_modified_all") {
    manager.internally_modified_all();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{true, true, true, true, true});
    CHECK(manager.is_modified_hydration() == false);
}

TEST_CASE_METHOD(fixture, "StateManager::externally_modified") {
    manager.externally_modified(2);
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, true, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.is_modified_hydration() == false);
    CHECK(manager.is_externally_modified(2) == true);
}

TEST_CASE_METHOD(fixture, "StateManager::internally_modified") {
    manager.internally_modified(2);
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, true, false, false});
    CHECK(manager.is_modified_hydration() == false);
    CHECK(manager.is_internally_modified(2) == true);
}

TEST_CASE_METHOD(fixture, "StateManager::modified_hydration_layer") {
    manager.modified_hydration_layer();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.is_modified_hydration() == true);
}

TEST_CASE("StateManager::reset") {
    unsigned int size = 5;
    StateManager manager(size);
    manager.reset_to_false();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.is_modified_hydration() == false);
}

TEST_CASE("StateManager::get_probe") {
    unsigned int size = 5;
    StateManager manager(size);

    for (unsigned int i = 0; i < size; i++) {
        auto probe = std::dynamic_pointer_cast<signaller::BoundSignaller>(manager.get_probe(i));
        CHECK(probe != nullptr);
        CHECK(probe->get_id() == i);
    }
}

TEST_CASE("StateManager::set_probe") {
    unsigned int size = 5;
    StateManager manager(size);

    SECTION("UnboundSignaller") {
        for (unsigned int i = 0; i < size; i++) {
            auto probe = std::make_shared<signaller::UnboundSignaller>();
            manager.set_probe(i, probe);
            CHECK(manager.get_probe(i) == probe);
        }
    }

    SECTION("BoundSignaller") {
        for (unsigned int i = 0; i < size; i++) {
            auto probe = std::make_shared<signaller::BoundSignaller>(i, &manager);
            manager.set_probe(i, probe);
            CHECK(manager.get_probe(i) == probe);
        }
    }
}

TEST_CASE_METHOD(fixture, "StateManager probe signalling") {
    manager.get_probe(2)->modified_external();
    manager.get_probe(4)->modified_external();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, true, false, true});
    CHECK(!manager.is_externally_modified(0));
    CHECK(!manager.is_externally_modified(1));
    CHECK(manager.is_externally_modified(2));
    CHECK(!manager.is_externally_modified(3));
    CHECK(manager.is_externally_modified(4));
}