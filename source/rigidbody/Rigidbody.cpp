#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <rigidbody/controller/ControllerFactory.h>
#include <fitter/SmartFitter.h>
#include <fitter/FitResult.h>
#include <data/Molecule.h>
#include <utility/Logging.h>

using namespace ausaxs::rigidbody;

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>

Rigidbody::~Rigidbody() = default;

Rigidbody::Rigidbody(data::Molecule&& _molecule) : molecule(std::move(_molecule)) {
    logging::log("setting histogram manager to PartialHistogramManagerMT.");
    molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
    controller = factory::create_controller(this);
    body_selector = factory::create_selection_strategy(this);
    transformer = factory::create_transform_strategy(this);
    parameter_generator = factory::create_parameter_strategy(
        this, 
        settings::rigidbody::iterations,
        5,
        std::numbers::pi/3
    );
    constraints = std::make_unique<constraints::ConstraintManager>(this);
}

Rigidbody::Rigidbody(Rigidbody&& other) = default;
Rigidbody& Rigidbody::operator=(Rigidbody&& other) = default;