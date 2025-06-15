#include <rigidbody/controller/SimpleController.h>
#include <rigidbody/parameters/ParameterGenerationStrategy.h>
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/selection/BodySelectStrategy.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/Rigidbody.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <grid/Grid.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::parameter;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::rigidbody::controller;

SimpleController::~SimpleController() = default;

void SimpleController::setup(const io::ExistingFile& measurement_path) {
    if (!calibration) {
        this->fitter = std::make_unique<fitter::ConstrainedFitter>(
            rigidbody->constraints.get(), 
            measurement_path, 
            rigidbody->molecule.get_histogram()
        );
    } else {
        auto histogram = rigidbody->molecule.get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        this->fitter = std::make_unique<fitter::ConstrainedFitter>(
            rigidbody->constraints.get(), 
            measurement_path, 
            std::move(histogram)
        );
    }
}

void SimpleController::update_fitter() {
    assert(fitter.get() != nullptr && "RigidBody::update_fitter: Fitter not initialized.");
    if (!calibration) {
        fitter->set_model(rigidbody->molecule.get_histogram());
    } else {
        auto histogram = rigidbody->molecule.get_histogram();
        histogram->apply_water_scaling_factor(calibration->get_parameter("c"));
        fitter->set_model(std::move(histogram));
    }
}

bool SimpleController::run_step() {
    auto& molecule = rigidbody->molecule;
    auto grid = molecule.get_grid();

    // select a body to be modified this iteration
    auto [ibody, iconstraint] = rigidbody->body_selector->next();
    if (iconstraint == -1) {    // transform free body
        Parameter param = rigidbody->parameter_generator->next(ibody);
        rigidbody->transformer->apply(std::move(param), ibody);
    } else {                    // transform constrained body
        DistanceConstraint& constraint = rigidbody->constraints->distance_constraints_map.at(ibody).at(iconstraint).get();
        Parameter param = rigidbody->parameter_generator->next(ibody);
        rigidbody->transformer->apply(std::move(param), constraint);
    }
    molecule.generate_new_hydration(); 

    // update the body location in the fitter
    update_fitter();
    double new_chi2 = fitter->fit_chi2_only();

    // if the old configuration was better
    if (new_chi2 >= best->chi2) {
        rigidbody->transformer->undo();         // undo the body transforms
        *grid = *best->grid;                    // restore the old grid
        molecule.get_waters() = best->waters;   // restore the old waters
        molecule.signal_modified_hydration_layer();
        return false;
    } else {
        // accept the changes
        best->grid = std::make_shared<grid::Grid>(*grid);
        best->waters = molecule.get_waters();
        best->chi2 = new_chi2;
        return true;
    }
}