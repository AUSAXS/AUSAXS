#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/Body.h>
#include <data/Molecule.h>
#include <io/ExistingFile.h>
#include <rigidbody/BodySplitter.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/detail/MoleculeTransformParametersAbsolute.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <settings/All.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::controller;

struct ControllerFixture {
    ControllerFixture() {
        settings::general::verbose = false;
        settings::molecule::implicit_hydrogens = false;
        settings::grid::min_bins = 250;

        // Create a rigidbody for testing
        auto bodies = BodySplitter::split("tests/files/SASDJG5.pdb");
        rb = std::make_unique<Rigidbody>(std::move(bodies));
        rb->constraints->generate_constraints(settings::rigidbody::ConstraintGenerationStrategyChoice::Backbone);
    }
    
    std::unique_ptr<Rigidbody> rb;
};

TEST_CASE_METHOD(ControllerFixture, "Controllers::SimpleController basic functionality") {
    SimpleController ctrl(rb.get());
    
    SECTION("Setup initializes controller") {
        REQUIRE_NOTHROW(ctrl.setup(io::ExistingFile("tests/files/SASDJG5.dat")));
        CHECK(ctrl.get_fitter() != nullptr);
        CHECK(ctrl.get_current_best_config() != nullptr);
    }
    
    SECTION("Prepare and finish step") {
        ctrl.setup(io::ExistingFile("tests/files/SASDJG5.dat"));
        
        // Run a few optimization steps
        for (int i = 0; i < 5; ++i) {
            bool accepted = ctrl.prepare_step();
            ctrl.finish_step();
            
            // Step acceptance is deterministic for SimpleController (better chi2 = accepted)
            // We just check it doesn't crash
            REQUIRE((accepted == true || accepted == false));
        }
    }
    
    SECTION("Best configuration is updated on improvement") {
        ctrl.setup(io::ExistingFile("tests/files/SASDJG5.dat"));
        
        double initial_chi2 = ctrl.get_current_best_config()->chi2;
        
        // Run several steps
        bool improvement_found = false;
        for (int i = 0; i < 100; ++i) {
            bool accepted = ctrl.prepare_step();
            ctrl.finish_step();
            
            if (accepted) {
                improvement_found = true;
                double new_chi2 = ctrl.get_current_best_config()->chi2;
                CHECK(new_chi2 <= initial_chi2);
                initial_chi2 = new_chi2;
                continue;
            }
        }
        
        // With 100 steps, at least one should be accepted
        CHECK(improvement_found);
    }

    SECTION("A rejected step leaves the controller describing the restored conformation") {
        ctrl.setup(io::ExistingFile("tests/files/SASDJG5.dat"));

        // run until a step is rejected
        bool rejected = false;
        for (int i = 0; i < 100 && !rejected; ++i) {
            rejected = !ctrl.prepare_step();
            ctrl.finish_step();
        }
        REQUIRE(rejected);

        // the poses were restored to the best configuration, and so must its chi2 be
        CHECK(rb->conformation->absolute_parameters.chi2 == ctrl.get_current_best_config()->chi2);

        // and after an update, the fitter must evaluate the restored molecule rather than the rejected candidate
        ctrl.update_fitter();
        fitter::ConstrainedFitter reference(rb->constraints.get(), io::ExistingFile("tests/files/SASDJG5.dat"), rb->molecule.get_histogram());
        CHECK_THAT(ctrl.get_fitter()->fit_chi2_only(), Catch::Matchers::WithinRel(reference.fit_chi2_only(), 1e-12));
    }
}

TEST_CASE("Controllers::ControllerFactory") {
    settings::general::verbose = false;
    settings::molecule::implicit_hydrogens = false;
    
    AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
    AtomFF a2({5, 0, 0}, form_factor::form_factor_t::C);
    Rigidbody rb(Molecule{std::vector<Body>{Body(std::vector{a1}), Body(std::vector{a2})}});
    
    SECTION("Create SimpleController") {
        auto ctrl = factory::create_controller(&rb, settings::rigidbody::ControllerChoice::Classic);
        REQUIRE(ctrl != nullptr);
        CHECK(dynamic_cast<SimpleController*>(ctrl.get()) != nullptr);
    }
}
