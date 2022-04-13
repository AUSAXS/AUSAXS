#include <catch2/catch.hpp>

#include <data/StateManager.h>
#include <data/Protein.h>

using std::vector, std::shared_ptr;

TEST_CASE("state_manager", "[managers]") {
    unsigned int size = 5;
    StateManager manager(size);

    CHECK(manager.get_externally_modified_bodies() == vector{true, true, true, true, true});
    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, false, false, false});
    manager.get_probe(2)->external_change();
    manager.get_probe(4)->external_change();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, true, false, true});
    CHECK(!manager.is_externally_modified(0));
    CHECK(!manager.is_externally_modified(1));
    CHECK(manager.is_externally_modified(2));
    CHECK(!manager.is_externally_modified(3));
    CHECK(manager.is_externally_modified(4));
}

TEST_CASE("partial_histogram_manager_works", "[managers]") {
    vector<Body> bodies(5);
    Protein protein(bodies);
    shared_ptr<PartialHistogramManager> phm = protein.get_histogram_manager();
    StateManager& manager = phm->get_state_manager();

    manager.reset();
    phm->get_probe(0)->external_change();
    phm->get_probe(2)->external_change();
    CHECK(manager.get_externally_modified_bodies() == vector{true, false, true, false, false});
}

TEST_CASE("protein_manager", "[managers]") {
    vector<Body> bodies(5);
    Protein protein(bodies);

    shared_ptr<PartialHistogramManager> phm = protein.get_histogram_manager();
    StateManager& manager = phm->get_state_manager();

    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, false, false, false});
    shared_ptr<StateManager::BoundSignaller> probe0 = phm->get_probe(0);
    shared_ptr<StateManager::BoundSignaller> probe2 = phm->get_probe(2);
    probe0->external_change();
    probe2->external_change();
    CHECK(manager.get_externally_modified_bodies() == vector{true, false, true, false, false});

    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, false, false, false});
    Body body;
    body.register_probe(probe0);
    body.changed_external_state();
    CHECK(manager.get_externally_modified_bodies() == vector{true, false, false, false, false});

    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, false, false, false});
    protein.bind_body_signallers();
    protein.bodies[0].changed_external_state();
    protein.bodies[2].changed_external_state();
    CHECK(manager.get_externally_modified_bodies() == vector{true, false, true, false, false});

    manager.reset();
    CHECK(manager.get_externally_modified_bodies() == vector{false, false, false, false, false});
    protein.bodies[0] = Body();
    protein.bodies[4] = Body();
    protein.bodies[0].changed_external_state();
    protein.bodies[4].changed_external_state();
    CHECK(manager.get_externally_modified_bodies() == vector{true, false, false, false, true});
}