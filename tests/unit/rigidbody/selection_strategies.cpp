#include <catch2/catch_test_macros.hpp>

#include <rigidbody/selection/RandomBodySelect.h>
#include <rigidbody/selection/SequentialBodySelect.h>
#include <rigidbody/selection/RandomConstraintSelect.h>
#include <rigidbody/selection/SequentialConstraintSelect.h>
#include <rigidbody/selection/ManualSelect.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::selection;

struct SelectionStrategiesFixture {
    SelectionStrategiesFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
        
        // Create test bodies
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({10, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a4({15, 0, 0}, form_factor::form_factor_t::C);
        
        Body b1(std::vector{a1});
        Body b2(std::vector{a2});
        Body b3(std::vector{a3});
        Body b4(std::vector{a4});
        
        rb = std::make_unique<Rigidbody>(Molecule{std::vector<Body>{b1, b2, b3, b4}});
        
        // Add some constraints for constraint-based selection
        rb->get_constraint_manager()->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(rb->get_molecule(), 0, 1)
        );
        rb->get_constraint_manager()->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(rb->get_molecule(), 1, 2)
        );
        rb->get_constraint_manager()->add_constraint(
            std::make_unique<constraints::DistanceConstraintCM>(rb->get_molecule(), 2, 3)
        );
    }
    
    std::unique_ptr<Rigidbody> rb;
};

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::RandomBodySelect") {
    RandomBodySelect selector(rb.get());
    
    SECTION("Select returns valid body indices") {
        for (int i = 0; i < 10; ++i) {
            auto selected = selector.select();
            CHECK(selected.size() == 1);
            CHECK(selected[0] >= 0);
            CHECK(selected[0] < rb->get_bodies().size());
        }
    }
    
    SECTION("Multiple selections can return different bodies") {
        std::set<unsigned int> seen_bodies;
        for (int i = 0; i < 20; ++i) {
            auto selected = selector.select();
            seen_bodies.insert(selected[0]);
        }
        // With 4 bodies and 20 selections, we should see more than 1 unique body
        CHECK(seen_bodies.size() > 1);
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::SequentialBodySelect") {
    SequentialBodySelect selector(rb.get());
    
    SECTION("Select cycles through bodies in order") {
        unsigned int num_bodies = rb->get_bodies().size();
        
        for (unsigned int i = 0; i < num_bodies * 2; ++i) {
            auto selected = selector.select();
            CHECK(selected.size() == 1);
            CHECK(selected[0] == i % num_bodies);
        }
    }
    
    SECTION("Select wraps around after last body") {
        unsigned int num_bodies = rb->get_bodies().size();
        
        // Advance to last body
        for (unsigned int i = 0; i < num_bodies - 1; ++i) {
            selector.select();
        }
        
        auto last = selector.select();
        CHECK(last[0] == num_bodies - 1);
        
        auto first_again = selector.select();
        CHECK(first_again[0] == 0);
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::RandomConstraintSelect") {
    RandomConstraintSelect selector(rb.get());
    
    SECTION("Select returns valid body indices based on constraints") {
        for (int i = 0; i < 10; ++i) {
            auto selected = selector.select();
            CHECK(selected.size() == 2);  // Should select two bodies involved in a constraint
            CHECK(selected[0] >= 0);
            CHECK(selected[0] < rb->get_bodies().size());
            CHECK(selected[1] >= 0);
            CHECK(selected[1] < rb->get_bodies().size());
            CHECK(selected[0] != selected[1]);  // Should be different bodies
        }
    }
    
    SECTION("Multiple selections can return different constraint pairs") {
        std::set<std::pair<unsigned int, unsigned int>> seen_pairs;
        for (int i = 0; i < 20; ++i) {
            auto selected = selector.select();
            auto pair = std::make_pair(
                std::min(selected[0], selected[1]),
                std::max(selected[0], selected[1])
            );
            seen_pairs.insert(pair);
        }
        // With 3 constraints and 20 selections, we should see more than 1 unique pair
        CHECK(seen_pairs.size() > 1);
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::SequentialConstraintSelect") {
    SequentialConstraintSelect selector(rb.get());
    
    SECTION("Select cycles through constraints in order") {
        unsigned int num_constraints = rb->get_constraint_manager()->get_constraints().size();
        
        for (unsigned int i = 0; i < num_constraints * 2; ++i) {
            auto selected = selector.select();
            CHECK(selected.size() == 2);
            CHECK(selected[0] != selected[1]);
        }
    }
    
    SECTION("Select returns pairs of bodies involved in constraints") {
        // First constraint is between bodies 0 and 1
        auto selected1 = selector.select();
        CHECK(((selected1[0] == 0 && selected1[1] == 1) || (selected1[0] == 1 && selected1[1] == 0)));
        
        // Second constraint is between bodies 1 and 2
        auto selected2 = selector.select();
        CHECK(((selected2[0] == 1 && selected2[1] == 2) || (selected2[0] == 2 && selected2[1] == 1)));
        
        // Third constraint is between bodies 2 and 3
        auto selected3 = selector.select();
        CHECK(((selected3[0] == 2 && selected3[1] == 3) || (selected3[0] == 3 && selected3[1] == 2)));
        
        // Should wrap around to first constraint
        auto selected4 = selector.select();
        CHECK(((selected4[0] == 0 && selected4[1] == 1) || (selected4[0] == 1 && selected4[1] == 0)));
    }
}

TEST_CASE_METHOD(SelectionStrategiesFixture, "SelectionStrategies::ManualSelect") {
    SECTION("Select specific bodies") {
        std::vector<unsigned int> bodies_to_select = {1, 3};
        ManualSelect selector(rb.get(), bodies_to_select);
        
        auto selected = selector.select();
        CHECK(selected == bodies_to_select);
    }
    
    SECTION("Select single body") {
        std::vector<unsigned int> bodies_to_select = {2};
        ManualSelect selector(rb.get(), bodies_to_select);
        
        auto selected = selector.select();
        CHECK(selected.size() == 1);
        CHECK(selected[0] == 2);
    }
    
    SECTION("Select all bodies") {
        std::vector<unsigned int> bodies_to_select = {0, 1, 2, 3};
        ManualSelect selector(rb.get(), bodies_to_select);
        
        auto selected = selector.select();
        CHECK(selected == bodies_to_select);
    }
    
    SECTION("Same selection every time") {
        std::vector<unsigned int> bodies_to_select = {0, 2};
        ManualSelect selector(rb.get(), bodies_to_select);
        
        auto selected1 = selector.select();
        auto selected2 = selector.select();
        auto selected3 = selector.select();
        
        CHECK(selected1 == selected2);
        CHECK(selected2 == selected3);
        CHECK(selected3 == bodies_to_select);
    }
}
