#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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

    SECTION("different sizes") {
        StateManager manager1(1);
        CHECK(manager1.get_probes().size() == 1);
        
        StateManager manager3(3);
        CHECK(manager3.get_probes().size() == 3);
        
        StateManager manager10(10);
        CHECK(manager10.get_probes().size() == 10);
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
    SECTION("single index") {
        manager.externally_modified(2);
        CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, true, false, false});
        CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.is_modified_hydration() == false);
        CHECK(manager.is_externally_modified(2) == true);
    }

    SECTION("multiple indices") {
        manager.externally_modified(0);
        manager.externally_modified(4);
        CHECK(manager.get_externally_modified_bodies() == std::vector{true, false, false, false, true});
        CHECK(manager.is_externally_modified(0) == true);
        CHECK(manager.is_externally_modified(4) == true);
    }
}

TEST_CASE_METHOD(fixture, "StateManager::internally_modified") {
    SECTION("single index") {
        manager.internally_modified(2);
        CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, true, false, false});
        CHECK(manager.is_modified_hydration() == false);
        CHECK(manager.is_internally_modified(2) == true);
    }

    SECTION("multiple indices") {
        manager.internally_modified(1);
        manager.internally_modified(3);
        CHECK(manager.get_internally_modified_bodies() == std::vector{false, true, false, true, false});
        CHECK(manager.is_internally_modified(1) == true);
        CHECK(manager.is_internally_modified(3) == true);
    }
}

TEST_CASE_METHOD(fixture, "StateManager::modified_hydration_layer") {
    manager.modified_hydration_layer();
    CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
    CHECK(manager.is_modified_hydration() == true);
}

TEST_CASE("StateManager::reset_to_false") {
    unsigned int size = 5;
    StateManager manager(size);
    
    SECTION("reset after modifications") {
        manager.externally_modified_all();
        manager.internally_modified_all();
        manager.modified_hydration_layer();
        
        manager.reset_to_false();
        CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.is_modified_hydration() == false);
    }

    SECTION("reset from initial state") {
        manager.reset_to_false();
        CHECK(manager.get_externally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.get_internally_modified_bodies() == std::vector{false, false, false, false, false});
        CHECK(manager.is_modified_hydration() == false);
    }
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

TEST_CASE("StateManager::get_probes") {
    unsigned int size = 3;
    StateManager manager(size);
    
    auto probes = manager.get_probes();
    CHECK(probes.size() == size);
    for (unsigned int i = 0; i < size; i++) {
        CHECK(probes[i] != nullptr);
    }
}

TEST_CASE_METHOD(fixture, "StateManager::is_externally_modified") {
    SECTION("initially false") {
        for (int i = 0; i < 5; i++) {
            CHECK(manager.is_externally_modified(i) == false);
        }
    }

    SECTION("after modification") {
        manager.externally_modified(2);
        CHECK(manager.is_externally_modified(2) == true);
        CHECK(manager.is_externally_modified(0) == false);
        CHECK(manager.is_externally_modified(4) == false);
    }
}

TEST_CASE_METHOD(fixture, "StateManager::is_internally_modified") {
    SECTION("initially false") {
        for (int i = 0; i < 5; i++) {
            CHECK(manager.is_internally_modified(i) == false);
        }
    }

    SECTION("after modification") {
        manager.internally_modified(3);
        CHECK(manager.is_internally_modified(3) == true);
        CHECK(manager.is_internally_modified(0) == false);
        CHECK(manager.is_internally_modified(4) == false);
    }
}
