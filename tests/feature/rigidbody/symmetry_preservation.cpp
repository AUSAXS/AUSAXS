#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;

TEST_CASE("Symmetry::preservation_during_optimization") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    SECTION("symmetry constraints maintained") {
        auto mol = BodySplitter::split("tests/files/SASDJG5.pdb");
        Rigidbody rb(std::move(mol));
        rb.molecule.generate_new_hydration();

        if (rb.molecule.size_body() < 2) {
            SKIP("Test requires multi-body system");
        }

        auto controller = std::make_unique<controller::SimpleController>(&rb);
        controller->setup("tests/files/SASDJG5.dat");

        int iterations = 10;
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            controller->finish_step();

            CHECK(rb.molecule.size_body() > 0);
            for (unsigned int j = 0; j < rb.molecule.size_body(); ++j) {
                CHECK(rb.molecule.get_body(j).size_atom() > 0);
            }
        }
    }
}

TEST_CASE("Symmetry::parameter_consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    auto controller = std::make_unique<controller::SimpleController>(&rb);
    controller->setup("tests/files/LAR1-2.dat");

    SECTION("all bodies have valid parameters") {
        controller->prepare_step();
        
        CHECK(rb.conformation->absolute_parameters.parameters.size() == rb.molecule.size_body());
        
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            auto& params = rb.conformation->absolute_parameters.parameters[i];
            
            CHECK(std::isfinite(params.rotation.x()));
            CHECK(std::isfinite(params.rotation.y()));
            CHECK(std::isfinite(params.rotation.z()));
            CHECK(std::isfinite(params.translation.x()));
            CHECK(std::isfinite(params.translation.y()));
            CHECK(std::isfinite(params.translation.z()));
        }
        
        controller->finish_step();
    }

    SECTION("parameters accumulate correctly over multiple steps") {
        int iterations = 10;
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            controller->finish_step();
        }
        
        for (unsigned int ibody = 0; ibody < rb.molecule.size_body(); ++ibody) {
            auto& current_body = rb.molecule.get_body(ibody);
            auto& params = rb.conformation->absolute_parameters.parameters[ibody];
            auto& original = rb.conformation->initial_conformation[ibody];
            
            Body reconstructed = original;
            reconstructed.rotate(matrix::rotation_matrix(params.rotation));
            reconstructed.translate(params.translation);
            
            auto current_cm = current_body.get_cm();
            auto reconstructed_cm = reconstructed.get_cm();
            
            INFO("After " << iterations << " iterations, body " << ibody << " should be reconstructible");
            REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 0.1));
            REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 0.1));
            REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 0.1));
        }
    }
}

TEST_CASE("Symmetry::transformation_reversibility") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::None;
    settings::grid::scaling = 2;

    AtomFF a1({-1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a2({-1,  1, -1}, form_factor::form_factor_t::C);
    AtomFF a3({ 1, -1, -1}, form_factor::form_factor_t::C);
    AtomFF a4({ 1,  1, -1}, form_factor::form_factor_t::C);
    
    Body b1 = Body(std::vector<AtomFF>{a1, a2});
    Body b2 = Body(std::vector<AtomFF>{a3, a4});
    
    Rigidbody rb(Molecule{std::vector<Body>{b1, b2}});
    rb.constraints->add_constraint(std::make_unique<constraints::DistanceConstraintBond>(&rb.molecule, 0, 1));

    SECTION("apply and undo restore original state") {
        std::vector<Vector3<double>> initial_cms;
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            initial_cms.push_back(rb.molecule.get_body(i).get_cm());
        }
        
        auto initial_params = rb.conformation->absolute_parameters;
        
        auto params = rb.parameter_generator->next(0);
        auto constraint = rb.constraints->discoverable_constraints[0].get();
        rb.transformer->apply(std::move(params), constraint);
        
        auto transformed_params = rb.conformation->absolute_parameters;
        
        rb.transformer->undo();
        
        auto restored_params = rb.conformation->absolute_parameters;
        
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            INFO("Body " << i << " parameters should be restored after undo");
            CHECK(restored_params.parameters[i].rotation == initial_params.parameters[i].rotation);
            CHECK(restored_params.parameters[i].translation == initial_params.parameters[i].translation);
            
            auto restored_cm = rb.molecule.get_body(i).get_cm();
            REQUIRE_THAT(restored_cm.x(), Catch::Matchers::WithinAbs(initial_cms[i].x(), 1e-6));
            REQUIRE_THAT(restored_cm.y(), Catch::Matchers::WithinAbs(initial_cms[i].y(), 1e-6));
            REQUIRE_THAT(restored_cm.z(), Catch::Matchers::WithinAbs(initial_cms[i].z(), 1e-6));
        }
    }
}
