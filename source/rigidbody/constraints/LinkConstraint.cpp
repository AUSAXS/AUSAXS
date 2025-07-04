// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/LinkConstraint.h>
#include <settings/RigidBodySettings.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

LinkConstraint::LinkConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2) 
    : protein(protein), ibody1(ibody1), ibody2(ibody2) {
    const Body& body1 = protein->get_body(ibody1);
    const Body& body2 = protein->get_body(ibody2);

    // constraints within the same body doesn't make sense
    if (body1.get_uid() == body2.get_uid()) {
        throw except::invalid_argument("LinkConstraint::LinkConstraint: Cannot create a constraint between atoms in the same body!");
    }
}

double LinkConstraint::evaluate() const {
    return 0;
}

const Body& LinkConstraint::get_body1() const {
    return protein->get_body(ibody1);
}

const Body& LinkConstraint::get_body2() const {
    return protein->get_body(ibody2);
}

Body& LinkConstraint::get_body1() {
    return protein->get_body(ibody1);
}

Body& LinkConstraint::get_body2() {
    return protein->get_body(ibody2);
}