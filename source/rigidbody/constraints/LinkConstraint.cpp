/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/LinkConstraint.h>
#include <settings/RigidBodySettings.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>

using namespace rigidbody::constraints;
using namespace data;
using namespace data::record;

LinkConstraint::LinkConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2) 
    : protein(protein), ibody1(ibody1), ibody2(ibody2) {
    const Body& body1 = protein->get_body(ibody1);
    const Body& body2 = protein->get_body(ibody2);

    // constraints within the same body doesn't make sense
    if (body1.get_id() == body2.get_id()) {
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