#include <catch2/catch_test_macros.hpp>

#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <data/Body.h>
#include <rigidbody/Rigidbody.h>
#include <settings/MoleculeSettings.h>
#include <settings/GeneralSettings.h>

using namespace ausaxs;
using namespace data;
using namespace rigidbody;

TEST_CASE("BodySelectStrategy::next") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    std::vector<AtomFF> b1 = {AtomFF({-1, -1, -1}, form_factor::form_factor_t::C), AtomFF({-1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b2 = {AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C), AtomFF({ 1, 1, -1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b3 = {AtomFF({-1, -1,  1}, form_factor::form_factor_t::C), AtomFF({-1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<AtomFF> b4 = {AtomFF({ 1, -1,  1}, form_factor::form_factor_t::C), AtomFF({ 1, 1,  1}, form_factor::form_factor_t::C)};
    std::vector<Body> atoms = {Body(b1), Body(b2), Body(b3), Body(b4)};
    Rigidbody rigidbody(Molecule{atoms});
    auto& manager = rigidbody.constraints;

    // add a varying number of constraints to each body
    for (unsigned int i = 0; i < rigidbody.molecule.size_body(); i++) {
        for (unsigned int j = i+1; j < rigidbody.molecule.size_body(); j++) {
            for (unsigned int k = j; k < 5; k++) {
                manager->add_constraint(std::make_unique<rigidbody::constraints::DistanceConstraintBond>(&rigidbody.molecule, i, j));
            }
        }
    }

    SECTION("SequentialConstraintSelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::SequentialConstraintSelect>(&rigidbody);
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> count;

        // check that the constraints are selected sequentially
        for (unsigned int i = 0; i < rigidbody.molecule.size_body(); i++) {
            for (unsigned int j = 0; j < rigidbody.constraints->get_body_constraints(i).size(); j++) {
                auto[ibody, iconstraint] = strat->next();
                REQUIRE(ibody == i);
                REQUIRE(iconstraint == int(j));
            }
        }

        // check that the strategy loops back to the beginning
        auto[ibody, iconstraint] = strat->next();
        REQUIRE(ibody == 0);
        REQUIRE(iconstraint == 0);
    }

    SECTION("SequentialBodySelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::SequentialBodySelect>(&rigidbody);
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> count;

        // check that the bodies are selected sequentially
        for (unsigned int i = 0; i < rigidbody.molecule.size_body(); i++) {
            auto[ibody, _] = strat->next();
            REQUIRE(ibody == i);
        }

        // check that the strategy loops back to the beginning
        auto[ibody, iconstraint] = strat->next();
        REQUIRE(ibody == 0);
    }

    SECTION("RandomConstraintSelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::RandomConstraintSelect>(&rigidbody);
        std::unordered_map<unsigned int, unsigned int> count;

        // count how many times each constraint is selected
        unsigned int iterations = 10000;
        for (unsigned int i = 0; i < iterations; i++) {
            auto[ibody, iconstraint] = strat->next();
            if (ibody >= rigidbody.molecule.size_body()) {
                std::cout << "Strategy selected a body outside the allowed range. Number: " << ibody << std::endl;
                REQUIRE(false);
            } 
            if (iconstraint >= int(rigidbody.constraints->get_body_constraints(ibody).size())) {
                std::cout << "Strategy selected a constraint outside the allowed range. Number: " << iconstraint << std::endl;
                REQUIRE(false);
            }

            auto constraint = rigidbody.constraints->get_body_constraints(ibody).at(iconstraint);
            for (unsigned int j = 0; j < rigidbody.constraints->discoverable_constraints.size(); j++) {
                if (rigidbody.constraints->discoverable_constraints[j].get() == constraint) {
                    count[j]++;
                }
            }
        }

        for (unsigned int i = 0; i < rigidbody.constraints->discoverable_constraints.size(); i++) {
            REQUIRE(count[i] > 0.8*iterations/rigidbody.constraints->discoverable_constraints.size());
        }
    }

    SECTION("RandomBodySelect::next") {
        std::unique_ptr<rigidbody::selection::BodySelectStrategy> strat = std::make_unique<rigidbody::selection::RandomBodySelect>(&rigidbody);
        std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> count;

        // count how many times each body and constraint is selected
        unsigned int iterations = 10000;
        for (unsigned int i = 0; i < iterations; i++) {
            auto[ibody, iconstraint] = strat->next();
            if (ibody >= rigidbody.molecule.size_body()) {
                std::cout << "Strategy selected a body outside the allowed range. Number: " << ibody << std::endl;
                REQUIRE(false);
            } 
            if (iconstraint >= int(rigidbody.constraints->get_body_constraints(ibody).size())) {
                std::cout << "Strategy selected a constraint outside the allowed range. Number: " << iconstraint << std::endl;
                REQUIRE(false);
            }
            count[ibody][iconstraint]++;
        }

        for (unsigned int i = 0; i < rigidbody.molecule.size_body(); i++) {
            // calculate how many times each body was selected
            double sum = 0;
            for (unsigned int j = 0; j < rigidbody.constraints->get_body_constraints(i).size(); j++) {
                sum += count[i][j];
            }

            // check that each body was selected at least 20% of the time
            REQUIRE(sum > iterations*0.2);

            // check that the constraints were randomly selected
            for (unsigned int j = i; j < rigidbody.constraints->get_body_constraints(i).size(); j++) {
                REQUIRE(count[i][j] > 0.7*sum/rigidbody.constraints->get_body_constraints(i).size());
            }
        }
    }
}