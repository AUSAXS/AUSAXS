#include <rigidbody/controller/IController.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/BestConf.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::controller;

std::unique_ptr<rigidbody::detail::BestConf> init_best(observer_ptr<Rigidbody> rigidbody) {
    return std::make_unique<rigidbody::detail::BestConf>(
        std::make_shared<grid::Grid>(*rigidbody->molecule.get_grid()), 
        rigidbody->molecule.get_waters(), 
        std::numeric_limits<double>::max()
    );
}

IController::IController(observer_ptr<Rigidbody> rigidbody) : rigidbody(rigidbody), best(init_best(rigidbody)) {
    assert(rigidbody != nullptr && "IController: RigidBody must not be null.");
}

IController::IController(observer_ptr<Rigidbody> rigidbody, std::unique_ptr<fitter::FitResult> calibration) 
    : rigidbody(rigidbody), calibration(std::move(calibration)), best(init_best(rigidbody))
{
    assert(rigidbody != nullptr && "IController: RigidBody must not be null.");
}

IController::~IController() = default;

observer_ptr<rigidbody::detail::BestConf> IController::current_best() const {
    assert(best != nullptr && "IController::current_best: Best configuration not set.");
    return best.get();
}

observer_ptr<fitter::ConstrainedFitter> IController::get_fitter() const {
    assert(fitter != nullptr && "IController::fitter: Fitter not set.");
    return fitter.get();
}

observer_ptr<const fitter::FitResult> IController::get_calibration() const {
    assert(calibration != nullptr && "IController::get_calibration: Calibration not set.");
    return calibration.get();
}