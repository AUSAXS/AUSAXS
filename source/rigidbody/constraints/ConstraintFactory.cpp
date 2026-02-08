// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/ConstraintFactory.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::factory;
using namespace ausaxs::rigidbody::sequencer::detail;

std::unique_ptr<constraints::Constraint> factory::create_constraint_cm(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1.body, body2.body, true);
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_closest(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1.body, body2.body, false);
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_attractor(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2, double target_distance
) {
    return std::make_unique<constraints::AttractorConstraint>(
        owner, target_distance, 
        body1.body, body2.body, std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica}
    );
}

std::unique_ptr<constraints::Constraint> factory::create_constraint(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2, unsigned int iatom1, unsigned int iatom2
) {
    return std::make_unique<constraints::DistanceConstraint>(owner, body1.body, body2.body, iatom1, iatom2);
}