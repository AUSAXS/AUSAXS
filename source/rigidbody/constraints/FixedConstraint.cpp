/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/FixedConstraint.h>
#include <settings/RigidBodySettings.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

FixedConstraint::FixedConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2) 
    : protein(protein), ibody1(ibody1), ibody2(ibody2) {
    const Body& body1 = protein->get_body(ibody1);
    const Body& body2 = protein->get_body(ibody2);

    // constraints within the same body doesn't make sense
    if (body1.get_uid() == body2.get_uid()) {
        throw except::invalid_argument("FixedConstraint::FixedConstraint: Cannot create a constraint between atoms in the same body!");
    }
}

double FixedConstraint::evaluate() const {
    return 0;
}

const Body& FixedConstraint::get_body1() const {
    return protein->get_body(ibody1);
}

const Body& FixedConstraint::get_body2() const {
    return protein->get_body(ibody2);
}

Body& FixedConstraint::get_body1() {
    return protein->get_body(ibody1);
}

Body& FixedConstraint::get_body2() {
    return protein->get_body(ibody2);
}