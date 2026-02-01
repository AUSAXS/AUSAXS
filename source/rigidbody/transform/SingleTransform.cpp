// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/Rigidbody.h>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::BodyTransformParameters&& par, constraints::DistanceConstraint& constraint) {
    TransformGroup group({&rigidbody->molecule.get_body(constraint.ibody1)}, {constraint.ibody1}, constraint, constraint.get_atom2().coordinates());
    backup(group);
    rotate(matrix::rotation_matrix(par.rotation), group);
    translate(par.translation, group);
    symmetry(std::move(par.symmetry_pars), *group.bodies.front());
}