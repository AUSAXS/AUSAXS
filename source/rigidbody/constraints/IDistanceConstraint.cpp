// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/constraints/IDistanceConstraint.h>
#include <data/Molecule.h>
#include <data/Body.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::data;
using namespace ausaxs::rigidbody::constraints;

namespace {
    // Rotate a body-frame offset into the body's current orientation. An empty orientation means the body has not been reoriented at all.
    Vector3<double> to_current_frame(const Body& body, const Vector3<double>& v) {
        const auto& F = body.symmetry().get_obj()->orientation;
        return F.has_value() ? F.value()*v : v;
    }
}

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

void IDistanceConstraint::cache_cm_offsets() {
    assert(0 <= iatom1 && 0 <= iatom2 && "IDistanceConstraint::cache_cm_offsets: the atom indices must be resolved first.");
    const auto& body1 = get_body1();
    const auto& body2 = get_body2();

    // constraints are built during setup, so the bodies are still in the frame their symmetries were defined in and the offsets need no conversion
    assert(
        !body1.symmetry().get_obj()->orientation.has_value() && !body2.symmetry().get_obj()->orientation.has_value()
        && "IDistanceConstraint::cache_cm_offsets: cannot build a constraint against an already-reoriented body."
    );
    cm_offset1 = body1.get_cm() - body1.get_atom(iatom1).coordinates();
    cm_offset2 = body2.get_cm() - body2.get_atom(iatom2).coordinates();
}

Vector3<double> IDistanceConstraint::cm1() const {
    const auto& body = get_body1();
    return body.get_atom(iatom1).coordinates() + to_current_frame(body, cm_offset1);
}

Vector3<double> IDistanceConstraint::cm2() const {
    const auto& body = get_body2();
    return body.get_atom(iatom2).coordinates() + to_current_frame(body, cm_offset2);
}