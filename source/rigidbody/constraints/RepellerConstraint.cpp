// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/RepellerConstraint.h>
#include <rigidbody/constraints/DistanceConstraintFunctions.h>
#include <data/Molecule.h>

using namespace ausaxs::rigidbody::constraints;

RepellerConstraint::RepellerConstraint(
    observer_ptr<const data::Molecule> molecule, double target_distance, int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : DistanceConstraintCM(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {
    d_target = target_distance;
}

double RepellerConstraint::evaluate() const {
    return transform(evaluate_distance(), d_target);
}

double RepellerConstraint::transform(double distance, double r_base) {
    if (r_base < distance) {return 0;}
    double offset = distance - r_base;
    return functions::attractor_repulsor(offset);
}