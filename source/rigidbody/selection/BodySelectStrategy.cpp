// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/selection/BodySelectStrategy.h>

#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::selection;

BodySelectStrategy::BodySelectStrategy(observer_ptr<const Rigidbody> rigidbody) : rigidbody(rigidbody), N(rigidbody->molecule.size_body()) {}