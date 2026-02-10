// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/IDistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody::constraints;

IDistanceConstraint::IDistanceConstraint(observer_ptr<const Molecule> molecule, int ibody1, int ibody2, std::pair<int, int> isym1, std::pair<int, int> isym2)
    : molecule(molecule), ibody1(ibody1), ibody2(ibody2), isym1(std::move(isym1)), isym2(std::move(isym2))
{}

IDistanceConstraint::IDistanceConstraint(
    observer_ptr<const data::Molecule> molecule, int ibody1, int iatom1, int ibody2, int iatom2, std::pair<int, int> isym1, std::pair<int, int> isym2
) : molecule(molecule), ibody1(ibody1), ibody2(ibody2), iatom1(iatom1), iatom2(iatom2), isym1(std::move(isym1)), isym2(std::move(isym2))
{}

const AtomFF& IDistanceConstraint::get_atom1() const {
    return molecule->get_body(ibody1).get_atom(iatom1);
}

const AtomFF& IDistanceConstraint::get_atom2() const {
    return molecule->get_body(ibody2).get_atom(iatom2);
}

const Body& IDistanceConstraint::get_body1() const {
    return molecule->get_body(ibody1);
}

const Body& IDistanceConstraint::get_body2() const {
    return molecule->get_body(ibody2);
}