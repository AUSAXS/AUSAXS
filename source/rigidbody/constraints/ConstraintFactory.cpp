// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/ConstraintFactory.h>
#include <rigidbody/constraints/DistanceConstraintCM.h>
#include <rigidbody/constraints/DistanceConstraintAtom.h>
#include <rigidbody/constraints/DistanceConstraintBond.h>
#include <rigidbody/constraints/AttractorConstraint.h>
#include <rigidbody/constraints/RepellerConstraint.h>
#include <data/Molecule.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody;
using namespace ausaxs::rigidbody::factory;
using namespace ausaxs::rigidbody::sequencer::detail;

std::unique_ptr<constraints::Constraint> factory::create_constraint_cm(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2
) {
    return std::make_unique<constraints::DistanceConstraintCM>(owner, body1.body, body2.body, std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica});
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_bond(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2
) {
    return std::make_unique<constraints::DistanceConstraintBond>(owner, body1.body, body2.body, std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica});
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_attractor(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2, double target_distance
) {
    return std::make_unique<constraints::AttractorConstraint>(
        owner, target_distance, 
        body1.body, body2.body, std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica}
    );
}

std::unique_ptr<constraints::Constraint> factory::create_constraint_repeller(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2, double target_distance
) {
    return std::make_unique<constraints::RepellerConstraint>(
        owner, target_distance, 
        body1.body, body2.body, std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica}
    );
}

std::unique_ptr<constraints::Constraint> factory::create_constraint(
    observer_ptr<const data::Molecule> owner, BodySymmetrySelector body1, BodySymmetrySelector body2, unsigned int iatom1, unsigned int iatom2
) {
    return std::make_unique<constraints::DistanceConstraintAtom>(
        owner, body1.body, iatom1, body2.body, iatom2, 
        std::pair{body1.symmetry, body1.replica}, std::pair{body2.symmetry, body2.replica}
    );
}