#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/RandomSelect.h>
#include <rigidbody/selection/SequentialSelect.h>
#include <data/record/Atom.h>
#include <data/Body.h>
#include <rigidbody/RigidBody.h>
#include <settings/MoleculeSettings.h>

using namespace data;
using namespace data::record;
using namespace rigidbody;

TEST_CASE("BodySelectStrategy::next") {
    settings::molecule::use_effective_charge = false;
    settings::molecule::implicit_hydrogens = false;
    std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1, -1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>(-1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, constants::atom_t::C, "C", 1), Atom(Vector3<double>( 1, 1,  1), 1, constants::atom_t::C, "C", 1)};
    std::vector<Body> atoms = {Body(b1), Body(b2), Body(b3), Body(b4)};
    RigidBody rigidbody(Molecule{atoms});
    auto manager = rigidbody.get_constraint_manager();

    // add a varying number of constraints to each body
    for (unsigned int i = 0; i < rigidbody.body_size(); i++) {
        for (unsigned int j = i+1; j < rigidbody.body_size(); j++) {
            for (unsigned int k = j; k < 5; k++) {
                manager->add_constraint(rigidbody::constraints::DistanceConstraint(&rigidbody, i, j, 0, 0));
            }
        }
    }

    SECTION("RandomSelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::RandomSelect>(&rigidbody);
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> count;

        // count how many times each body and constraint is selected
        unsigned int iterations = 1000;
        for (unsigned int i = 0; i < iterations; i++) {
            auto[ibody, iconstraint] = strat->next();
            if (ibody >= rigidbody.body_size()) {
                std::cout << "Strategy selected a body outside the allowed range. Number: " << ibody << std::endl;
                REQUIRE(false);
            } 
            if (iconstraint >= rigidbody.get_constraint_manager()->distance_constraints_map.at(ibody).size()) {
                std::cout << "Strategy selected a constraint outside the allowed range. Number: " << iconstraint << std::endl;
                REQUIRE(false);
            }
            count[ibody][iconstraint]++;
        }

        for (unsigned int i = 0; i < rigidbody.body_size(); i++) {
            // calculate how many times each body was selected
            double sum = 0;
            for (unsigned int j = 0; j < rigidbody.get_constraint_manager()->distance_constraints_map.at(i).size(); j++) {
                sum += count[i][j];
            }

            // check that each body was selected at least 20% of the time
            REQUIRE(sum > iterations*0.2);

            // check that the constraints were randomly selected
            for (unsigned int j = i; j < rigidbody.get_constraint_manager()->distance_constraints_map.at(i).size(); j++) {
                REQUIRE(count[i][j] > 0.7*sum/rigidbody.get_constraint_manager()->distance_constraints_map.at(i).size());
            }
        }
    }

    SECTION("SequentialSelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::SequentialSelect>(&rigidbody);
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> count;

        // check that the constraints are selected sequentially
        for (unsigned int i = 0; i < rigidbody.body_size(); i++) {
            for (unsigned int j = 0; j < rigidbody.get_constraint_manager()->distance_constraints_map.at(i).size(); j++) {
                auto[ibody, iconstraint] = strat->next();
                REQUIRE(ibody == i);
                REQUIRE(iconstraint == j);
            }
        }

        // check that the strategy loops back to the beginning
        auto[ibody, iconstraint] = strat->next();
        REQUIRE(ibody == 0);
        REQUIRE(iconstraint == 0);
    }

    SECTION("RandomConstraintSelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::RandomConstraintSelect>(&rigidbody);
        std::unordered_map<unsigned int, unsigned int> count;

        // count how many times each constraint is selected
        unsigned int iterations = 1000;
        for (unsigned int i = 0; i < iterations; i++) {
            auto[ibody, iconstraint] = strat->next();
            if (ibody >= rigidbody.body_size()) {
                std::cout << "Strategy selected a body outside the allowed range. Number: " << ibody << std::endl;
                REQUIRE(false);
            } 
            if (iconstraint >= rigidbody.get_constraint_manager()->distance_constraints_map.at(ibody).size()) {
                std::cout << "Strategy selected a constraint outside the allowed range. Number: " << iconstraint << std::endl;
                REQUIRE(false);
            }
            count[iconstraint]++;
        }

        for (unsigned int i = 0; i < rigidbody.get_constraint_manager()->distance_constraints.size(); i++) {
            REQUIRE(count[i] > 0.8*iterations/rigidbody.body_size());
        }
    }
}