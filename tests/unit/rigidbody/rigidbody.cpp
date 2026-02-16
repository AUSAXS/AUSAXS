#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <grid/Grid.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

struct fixture {
    fixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::molecule::center = false;
        settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    }

    AtomFF a1 = AtomFF({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2 = AtomFF({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3 = AtomFF({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a4 = AtomFF({ 1,  1, -1}, form_factor::form_factor_t::C);
    Body b1 = Body(std::vector{a1, a2});
    Body b2 = Body(std::vector{a3, a4});
};

TEST_CASE_METHOD(fixture, "Rigidbody::construction") {
    SECTION("move construction from Molecule") {
        Molecule mol(std::vector<Body>{b1, b2});
        Rigidbody rb(std::move(mol));
        
        CHECK(rb.molecule.size_body() == 2);
        CHECK(rb.constraints != nullptr);
        CHECK(rb.conformation != nullptr);
        CHECK(rb.body_selector != nullptr);
        CHECK(rb.transformer != nullptr);
        CHECK(rb.parameter_generator != nullptr);
    }

    SECTION("move assignment") {
        Molecule mol1(std::vector<Body>{b1});
        Molecule mol2(std::vector<Body>{b2});
        Rigidbody rb1(std::move(mol1));
        Rigidbody rb2(std::move(mol2));
        
        rb1 = std::move(rb2);
        CHECK(rb1.molecule.size_body() == 1);
    }
}

TEST_CASE_METHOD(fixture, "Rigidbody::refresh_grid") {
    settings::grid::scaling = 2;
    Molecule mol(std::vector<Body>{b1, b2});
    Rigidbody rb(std::move(mol));
    
    auto initial_grid = rb.molecule.get_grid();
    auto initial_axes = initial_grid->get_axes();
    
    rb.molecule.get_body(0).translate(Vector3<double>(10, 0, 0));
    
    rb.refresh_grid();
    
    auto new_grid = rb.molecule.get_grid();
    auto new_axes = new_grid->get_axes();
    
    CHECK(new_axes.x.span() >= initial_axes.x.span());
}

TEST_CASE_METHOD(fixture, "Rigidbody::conformation initialization") {
    Molecule mol(std::vector<Body>{b1, b2});
    Rigidbody rb(std::move(mol));
    
    SECTION("initial conformation is stored") {
        CHECK(rb.conformation->initial_conformation.size() == 2);
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            CHECK(rb.conformation->initial_conformation[i].size_atom() == rb.molecule.get_body(i).size_atom());
        }
    }

    SECTION("absolute parameters initialized") {
        CHECK(rb.conformation->absolute_parameters.parameters.size() == 2);
        for (const auto& param : rb.conformation->absolute_parameters.parameters) {
            CHECK(param.rotation == Vector3<double>(0, 0, 0));
        }
    }
}

TEST_CASE("Rigidbody::component creation from factories") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));

    SECTION("constraints created") {
        CHECK(rb.constraints != nullptr);
        CHECK(rb.constraints->discoverable_constraints.size() > 0);
    }

    SECTION("body selector created") {
        CHECK(rb.body_selector != nullptr);
        auto [ibody, iconstraint] = rb.body_selector->next();
        CHECK(ibody >= 0);
        CHECK(ibody < static_cast<int>(rb.molecule.size_body()));
    }

    SECTION("transformer created") {
        CHECK(rb.transformer != nullptr);
    }

    SECTION("parameter generator created") {
        CHECK(rb.parameter_generator != nullptr);
        auto params = rb.parameter_generator->next(0);
        bool has_params = params.rotation.has_value() || params.translation.has_value();
        CHECK(has_params);
    }
}
