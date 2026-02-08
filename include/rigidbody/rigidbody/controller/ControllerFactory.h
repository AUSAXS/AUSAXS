// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/controller/IController.h>
#include <settings/RigidBodySettings.h>

#include <memory>

namespace ausaxs::rigidbody::factory {
    std::unique_ptr<controller::IController> create_controller(observer_ptr<Rigidbody> molecule);
    std::unique_ptr<controller::IController> create_controller(observer_ptr<Rigidbody> molecule, settings::rigidbody::ControllerChoice choice);
}