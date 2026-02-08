// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/AttractorConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

AttractorConstraint::AttractorConstraint(
    observer_ptr<const data::Molecule> molecule, double target_distance, int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : DistanceConstraintCM(molecule, ibody1, ibody2, std::move(isym1), std::move(isym2)) {
    d_target = target_distance;
}

double AttractorConstraint::evaluate() const {
    return transform(evaluate_distance(), d_target);
}

double AttractorConstraint::transform(double distance, double r_base) {
    if (distance < r_base) {return 0;}
    double offset = distance - r_base;
    return 10*offset*offset;
}