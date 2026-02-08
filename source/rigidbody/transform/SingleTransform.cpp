// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/constraints/DistanceConstraint.h>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<RigidBody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::Parameter&& par, constraints::DistanceConstraint& constraint) {
    // Check if this is effectively a no-op to avoid floating point drift
    bool is_zero_rotation = par.rotation.magnitude() < 1e-10;
    bool is_zero_translation = par.translation.magnitude() < 1e-10;
    
    if (is_zero_rotation && is_zero_translation) {
        // Still need to backup in case undo is called
        TransformGroup group({&constraint.get_body1()}, {constraint.ibody1}, constraint, constraint.get_atom2().coordinates());
        backup(group);
        return;
    }
    
    TransformGroup group({&constraint.get_body1()}, {constraint.ibody1}, constraint, constraint.get_atom2().coordinates());
    backup(group);
    rotate(matrix::rotation_matrix(par.rotation), group);
    translate(par.translation, group);
    symmetry(std::move(par.symmetry_pars), *group.bodies.front());
}