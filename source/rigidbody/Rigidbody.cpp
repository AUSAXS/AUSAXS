#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/selection/BodySelectFactory.h>
#include <rigidbody/transform/TransformFactory.h>
#include <rigidbody/parameters/ParameterGenerationFactory.h>
#include <fitter/SmartFitter.h>
#include <fitter/FitResult.h>
#include <data/Molecule.h>

using namespace ausaxs::rigidbody;

Rigidbody::~Rigidbody() = default;

Rigidbody::Rigidbody(data::Molecule&& molecule)
    : molecule(std::move(molecule)),
      constraints(std::make_unique<constraints::ConstraintManager>(&this->molecule)),
      body_selector(factory::create_selection_strategy(this)),
      transformer(factory::create_transform_strategy(this)),
      parameter_generator(factory::create_parameter_strategy(this, settings::rigidbody::iterations, 5, std::numbers::pi/3)) 
{
        molecule.set_histogram_manager(settings::hist::HistogramManagerChoice::PartialHistogramManagerMT);
}