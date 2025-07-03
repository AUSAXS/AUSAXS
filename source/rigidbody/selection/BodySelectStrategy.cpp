// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectStrategy.h>

#include <rigidbody/RigidBody.h>

using namespace ausaxs::rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(observer_ptr<const RigidBody> rigidbody) : rigidbody(rigidbody), N(rigidbody->size_body()) {}