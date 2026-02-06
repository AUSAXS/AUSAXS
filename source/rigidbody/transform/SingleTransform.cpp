// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::BodyTransformParametersRelative&& par, constraints::DistanceConstraint& constraint) {
    auto& body = rigidbody->molecule.get_body(constraint.ibody1);
    TransformGroup group({&body}, {constraint.ibody1}, constraint, constraint.get_atom2().coordinates());
    backup(group);
    body = rigidbody->conformation->initial_conformation[constraint.ibody1];

    // remove body from grid since it does not track transforms
    auto grid = rigidbody->molecule.get_grid();
    grid->remove(body);

    // compute new absolute transform parameters for the body
    auto& body_params = rigidbody->conformation->absolute_parameters.parameters[constraint.ibody1];
    if (par.rotation.has_value()) {body_params.rotation += par.rotation.value();}
    if (par.translation.has_value()) {body_params.translation += par.translation.value();}

    // apply transformations
    if (par.rotation.has_value() || par.translation.has_value()) {
        rotate_and_translate(matrix::rotation_matrix(body_params.rotation), body_params.translation, group);
    }

    // update and apply symmetry parameters
    if (par.symmetry_pars.has_value()) {
        body_params.symmetry_pars = add_symmetries(body_params.symmetry_pars, par.symmetry_pars.value());
        apply_symmetry(body_params.symmetry_pars, *group.bodies.front());
    }

    // re-add body and refresh grid
    grid->add(body);
    rigidbody->refresh_grid();
}