// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>

#include <data/Body.h>
#include <grid/Grid.h>
#include <math/MatrixUtils.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/transform/BackupBody.h>  // IWYU pragma: keep

#include <cassert>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(
    const parameter::BodyTransformParametersRelative& par, observer_ptr<const constraints::IDistanceConstraint> constraint, int ibody
) {
    // the selected body is the one that moves, pivoting about the constraint's atom on the *other* body. The constraint was drawn from ibody's own list,
    // so ibody is one of its two ends; which end it is depends only on the order the constraint was declared in, so it must not decide what moves
    assert((constraint->ibody1 == ibody || constraint->ibody2 == ibody) && "SingleTransform::apply: the constraint does not involve the selected body.");
    Vector3<double> pivot = (constraint->ibody1 == ibody ? constraint->get_atom2() : constraint->get_atom1()).coordinates();

    // remove body from grid since it does not track transforms
    auto* grid = rigidbody->molecule.get_grid();
    {   // backup body and parameters for undo
        auto& body = rigidbody->molecule.get_body(ibody);
        grid->remove(body);
        bodybackup.clear();
        bodybackup.emplace_back(body, ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);
    }

    // compute new absolute transform parameters for the body
    auto& body_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
    if (par.rotation.has_value()) {
        auto dR = matrix::rotation_matrix(par.rotation.value());
        body_params.transform(pivot, dR);
    } if (par.translation.has_value()) {
        body_params.transform(par.translation.value());
    }

    // reconstruct body from initial conformation using absolute parameters
    auto& body = rigidbody->molecule.get_body(ibody);
    if (par.rotation.has_value() || par.translation.has_value()) {
        body = rigidbody->conformation->initial_conformation[ibody];
        rotate_and_translate(matrix::rotation_matrix(body_params.rotation), body_params.translation, body.get_cm(), body);
        restore_symmetry(ibody); // rebuilding from the initial conformation also reset the symmetries, so put the accumulated values back
    } // no transformation: nothing has touched the body, so it is already in the right state and is left alone

    // the symmetry deltas were generated from ibody's own symmetry list, and ibody is the body being moved
    if (par.symmetry_pars.has_value()) {apply_symmetry_delta(ibody, par.symmetry_pars.value());}

    // re-add body and refresh grid
    rigidbody->refresh_grid();
    rigidbody->molecule.get_grid()->add(body); // refresh_grid may reallocate the grid, so re-fetch the pointer
}