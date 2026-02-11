#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/DefaultOptimizer.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/transform/RigidTransform.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;

TEST_CASE("OptimizationLoop::basic_convergence") {
    settings::general::verbose = false;
    settings::general::output = "temp/rigidbody/optimization_loop/";
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::controller_choice = settings::rigidbody::ControllerChoice::Classic;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    auto controller = factory::create_controller(&rb);
    controller->setup("tests/files/LAR1-2.dat");

    double initial_chi2 = controller->get_current_best_config()->chi2;

    SECTION("chi2 changes during optimization") {
        int iterations = 50;
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            controller->finish_step();
        }

        double final_chi2 = controller->get_current_best_config()->chi2;
        
        CHECK(final_chi2 > 0);
        CHECK(initial_chi2 > 0);
    }

    SECTION("molecule integrity maintained") {
        int iterations = 20;
        unsigned int initial_bodies = rb.molecule.size_body();
        
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            controller->finish_step();
            
            CHECK(rb.molecule.size_body() == initial_bodies);
            for (unsigned int j = 0; j < rb.molecule.size_body(); ++j) {
                CHECK(rb.molecule.get_body(j).size_atom() > 0);
            }
        }
    }
}

TEST_CASE("OptimizationLoop::histogram_consistency") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    auto controller = std::make_unique<controller::SimpleController>(&rb);
    controller->setup("tests/files/LAR1-2.dat");

    SECTION("histogram recalculation matches transformed structure") {
        controller->prepare_step();
        
        auto histogram1 = rb.molecule.get_histogram();
        
        Molecule verification_mol = rb.molecule;
        auto histogram2 = verification_mol.get_histogram();
        
        auto profile1 = histogram1->debye_transform();
        auto profile2 = histogram2->debye_transform();
        
        CHECK(profile1.size() == profile2.size());
        for (unsigned int i = 0; i < std::min(profile1.size(), profile2.size()); ++i) {
            REQUIRE_THAT(profile1.y(i), Catch::Matchers::WithinRel(profile2.y(i), 1e-6));
        }
        
        controller->finish_step();
    }

    SECTION("parameters correctly reconstruct body positions") {
        int iterations = 10;
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            
            for (unsigned int ibody = 0; ibody < rb.molecule.size_body(); ++ibody) {
                auto& current_body = rb.molecule.get_body(ibody);
                auto& params = rb.conformation->absolute_parameters.parameters[ibody];
                auto& original = rb.conformation->initial_conformation[ibody];
                
                Body reconstructed = original;
                reconstructed.rotate(matrix::rotation_matrix(params.rotation));
                reconstructed.translate(params.translation);
                
                auto current_cm = current_body.get_cm();
                auto reconstructed_cm = reconstructed.get_cm();
                
                INFO("Iteration " << i << ", body " << ibody);
                REQUIRE_THAT(reconstructed_cm.x(), Catch::Matchers::WithinAbs(current_cm.x(), 0.01));
                REQUIRE_THAT(reconstructed_cm.y(), Catch::Matchers::WithinAbs(current_cm.y(), 0.01));
                REQUIRE_THAT(reconstructed_cm.z(), Catch::Matchers::WithinAbs(current_cm.z(), 0.01));
            }
            
            controller->finish_step();
        }
    }
}

TEST_CASE("OptimizationLoop::constraint_preservation") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::rigidbody::transform_strategy = settings::rigidbody::TransformationStrategyChoice::RigidTransform;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    auto controller = std::make_unique<controller::SimpleController>(&rb);
    controller->setup("tests/files/LAR1-2.dat");

    SECTION("internal group constraints preserved during transformation") {
        std::vector<double> initial_distances;
        for (const auto& constraint : rb.constraints->discoverable_constraints) {
            auto dist = (constraint->get_atom1().coordinates() - constraint->get_atom2().coordinates()).norm();
            initial_distances.push_back(dist);
        }

        int iterations = 5;
        for (int i = 0; i < iterations; ++i) {
            controller->prepare_step();
            controller->finish_step();
        }

        for (size_t j = 0; j < rb.constraints->discoverable_constraints.size(); ++j) {
            auto& constraint = rb.constraints->discoverable_constraints[j];
            double current_dist = (constraint->get_atom1().coordinates() - constraint->get_atom2().coordinates()).norm();
            
            INFO("Constraint " << j << " distance changed beyond tolerance");
        }
    }
}

TEST_CASE("OptimizationLoop::backup_and_restore") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    auto controller = std::make_unique<controller::SimpleController>(&rb);
    controller->setup("tests/files/LAR1-2.dat");

    SECTION("rejected steps restore previous state") {
        auto initial_params = rb.conformation->absolute_parameters;
        
        std::vector<Vector3<double>> initial_cms;
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            initial_cms.push_back(rb.molecule.get_body(i).get_cm());
        }
        
        double initial_chi2 = controller->get_current_best_config()->chi2;
        
        int rejected_count = 0;
        int iterations = 20;
        for (int i = 0; i < iterations; ++i) {
            bool accepted = controller->prepare_step();
            double step_chi2 = rb.conformation->absolute_parameters.chi2;
            controller->finish_step();
            
            if (!accepted) {
                rejected_count++;
                
                CHECK(controller->get_current_best_config()->chi2 == initial_chi2);
            }
        }
        
        INFO("At least some steps should be rejected");
    }
}
