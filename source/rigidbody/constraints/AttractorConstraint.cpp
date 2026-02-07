// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/AttractorConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

AttractorConstraint::AttractorConstraint(observer_ptr<const data::Molecule> molecule, unsigned int ibody1, unsigned int ibody2, double target) 
    : DistanceConstraint(molecule, ibody1, ibody2, true) 
{
    r_base = target;
}

double AttractorConstraint::evaluate() const {
    const AtomFF& atom1 = molecule->get_body(ibody1).get_atom(iatom1);
    const AtomFF& atom2 = molecule->get_body(ibody2).get_atom(iatom2);
    return transform(atom1.coordinates().distance(atom2.coordinates()), r_base);
}

bool AttractorConstraint::operator==(const AttractorConstraint& constraint) const {
    return DistanceConstraint::operator==(constraint);
}

double AttractorConstraint::transform(double distance, double r_base) {
    if (distance < r_base) {return 0;}
    double offset = distance - r_base;
    return 10*offset*offset;
}

std::string AttractorConstraint::to_string() const {
    std::stringstream ss;
    ss << 
        "AttractorConstraint between (" << form_factor::to_string(molecule->get_body(ibody1).get_atom(iatom1).form_factor_type()) << ") and "
        "(" << form_factor::to_string(molecule->get_body(ibody2).get_atom(iatom2).form_factor_type()) << ")"
    ;
    return ss.str();
}

std::ostream& operator<<(std::ostream& os, const AttractorConstraint& constraint) {os << constraint.to_string(); return os;}