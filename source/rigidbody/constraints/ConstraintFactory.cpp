// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/ConstraintFactory.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::factory;

std::unique_ptr<constraints::Constraint> factory::create_constraint_cm(
    observer_ptr<const data::Molecule> owner, unsigned int body1, unsigned int body2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1, body2, true);
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_closest(
    observer_ptr<const data::Molecule> owner, unsigned int body1, unsigned int body2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1, body2, false);
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_attractor(
    observer_ptr<const data::Molecule> owner, unsigned int body1, unsigned int body2, double target_distance
) {
    return std::make_unique<constraints::AttractorConstraint>(owner, body1, body2, target_distance);
}

std::unique_ptr<constraints::Constraint> factory::create_constraint(
    observer_ptr<const data::Molecule> owner, unsigned int body1, unsigned int body2, unsigned int iatom1, unsigned int iatom2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1, body2, iatom1, iatom2);
}