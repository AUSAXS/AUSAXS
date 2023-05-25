#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Protein.h>
#include <data/Body.h>
#include <data/state/StateManager.h>
#include <data/state/Signaller.h>
#include <data/state/BoundSignaller.h>
#include <data/state/UnboundSignaller.h>

#include <memory>

TEST_CASE("StateManager_constructor") {
    unsigned int size = 5;
    StateManager manager(size);
    REQUIRE(manager.get_internally_modified_bodies().size() == size);
    CHECK(manager.get_probes().size() == size);
    CHECK(manager.get_externally_modified_bodies() == std::vector{true, true, true, true, true});
    CHECK(manager.get_internally_modified_bodies() == std::vector{true, true, true, true, true});
    CHECK(manager.get_modified_hydration() == true);
}

TEST_CASE("StateManager_reset") {
    unsigned int size = 5;
    StateManager manager(size);
    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_modified_hydration() == false);
}

TEST_CASE("StateManager_get_probe") {
    unsigned int size = 5;
    StateManager manager(size);

    for (unsigned int i = 0; i < size; i++) {
        auto probe = std::dynamic_pointer_cast<signaller::BoundSignaller>(manager.get_probe(i));
        CHECK(probe != nullptr);
        CHECK(probe->get_id() == i);
    }
}

TEST_CASE("StateManager_signalling") {
    unsigned int size = 5;
    StateManager manager(size);

    manager.get_probe(2)->external_change();
    manager.get_probe(4)->external_change();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, true, false, true});
    CHECK(!manager.is_externally_modified(0));
    CHECK(!manager.is_externally_modified(1));
    CHECK(manager.is_externally_modified(2));
    CHECK(!manager.is_externally_modified(3));
    CHECK(manager.is_externally_modified(4));
}