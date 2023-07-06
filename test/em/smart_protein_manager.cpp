#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <em/manager/SmartProteinManager.h>
#include <em/ImageStack.h>
#include <data/Protein.h>

struct fixture {
    em::ImageStack stack = em::ImageStack("test/file/A2M_2020_Q4.ccp4");
    em::managers::SmartProteinManager manager = em::managers::SmartProteinManager(&stack);
};

TEST_CASE_METHOD(fixture, "SmartProteinManager::ProteinManager") {
    SECTION("ImageStackBase*") {
        CHECK(!manager.get_charge_levels().empty());
    }
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::set_charge_levels") {
    std::vector<double> charge_levels = {1, 2, 3};
    manager.set_charge_levels(charge_levels);
    CHECK(manager.get_charge_levels() == charge_levels);
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::get_protein") {
    // we just check the size of the returned protein
    unsigned int size = manager.get_protein(2)->atom_size();
    CHECK(size != 0);
    CHECK(manager.get_protein(1)->atom_size() > size);
    CHECK(manager.get_protein(3)->atom_size() < size);
    CHECK(manager.get_protein(2)->atom_size() == size);
}

TEST_CASE_METHOD(fixture, "SmartProteinManager::get_histogram") {
    CHECK(manager.get_histogram(1) == manager.get_protein(1)->get_histogram());
}