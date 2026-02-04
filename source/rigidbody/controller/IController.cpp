// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/controller/IController.h>
#include <rigidbody/constraints/ConstrainedFitter.h>
#include <rigidbody/detail/Configuration.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>

#include <memory>
#include <cassert>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::controller;

std::unique_ptr<rigidbody::detail::Configuration> init_config() {
    return std::make_unique<rigidbody::detail::Configuration>();
}

IController::IController(observer_ptr<Rigidbody> rigidbody) : rigidbody(rigidbody), current_config(init_config()), current_best_config(init_config()) {
    assert(rigidbody != nullptr && "IController: RigidBody must not be null.");
}

IController::IController(observer_ptr<Rigidbody> rigidbody, std::unique_ptr<fitter::FitResult> calibration) 
    : rigidbody(rigidbody), calibration(std::move(calibration)), current_config(init_config())
{
    assert(rigidbody != nullptr && "IController: RigidBody must not be null.");
}

IController::~IController() = default;

observer_ptr<rigidbody::detail::Configuration> IController::get_current_best_config() const {
    assert(current_best_config != nullptr && "IController::current_best: Best configuration not set.");
    return current_best_config.get();
}

observer_ptr<rigidbody::detail::Configuration> IController::get_current_config() const {
    assert(current_config != nullptr && "IController::current_best: Current configuration not set.");
    return current_config.get();
}

observer_ptr<fitter::ConstrainedFitter> IController::get_fitter() const {
    assert(fitter != nullptr && "IController::fitter: Fitter not set.");
    return fitter.get();
}

observer_ptr<const fitter::FitResult> IController::get_calibration() const {
    assert(calibration != nullptr && "IController::get_calibration: Calibration not set.");
    return calibration.get();
}