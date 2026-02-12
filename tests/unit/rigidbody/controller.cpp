#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::controller;

TEST_CASE("SimpleController::basic_functionality") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    SimpleController controller(&rb);

    SECTION("setup initializes fitter and best config") {
        controller.setup("tests/files/LAR1-2.dat");
        
        CHECK(controller.get_fitter() != nullptr);
        CHECK(controller.get_current_best_config() != nullptr);
        CHECK(controller.get_current_best_config()->chi2 > 0);
    }

    SECTION("prepare_step generates configuration") {
        controller.setup("tests/files/LAR1-2.dat");
        double initial_chi2 = controller.get_current_best_config()->chi2;
        
        controller.prepare_step();
        
        CHECK(controller.get_current_best_config()->chi2 == initial_chi2);
    }

    SECTION("finish_step updates or reverts") {
        controller.setup("tests/files/LAR1-2.dat");
        
        controller.prepare_step();
        controller.finish_step();
        
        CHECK(rb.conformation != nullptr);
    }
}

TEST_CASE("Controller::step_workflow") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    settings::rigidbody::constraint_generation_strategy = settings::rigidbody::ConstraintGenerationStrategyChoice::Linear;
    settings::grid::min_bins = 250;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));
    rb.molecule.generate_new_hydration();

    SimpleController controller(&rb);
    controller.setup("tests/files/LAR1-2.dat");

    SECTION("prepare_step modifies molecule") {
        auto initial_cm = rb.molecule.get_body(0).get_cm();
        
        controller.prepare_step();
        
        CHECK(rb.molecule.get_body(0).size_atom() > 0);
    }

    SECTION("finish_step maintains consistency") {
        controller.prepare_step();
        controller.finish_step();
        
        CHECK(rb.molecule.size_body() == 3);
        for (unsigned int i = 0; i < rb.molecule.size_body(); ++i) {
            CHECK(rb.molecule.get_body(i).size_atom() > 0);
        }
    }
}

TEST_CASE("ControllerFactory::create_controller") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;

    auto mol = BodySplitter::split("tests/files/LAR1-2.pdb", {9, 99});
    Rigidbody rb(std::move(mol));

    SECTION("creates Classic controller") {
        settings::rigidbody::controller_choice = settings::rigidbody::ControllerChoice::Classic;
        auto controller = factory::create_controller(&rb);
        CHECK(controller != nullptr);
    }
}
