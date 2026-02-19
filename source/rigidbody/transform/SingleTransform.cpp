// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/constraints/IDistanceConstraint.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::BodyTransformParametersRelative&& par, observer_ptr<const constraints::IDistanceConstraint> constraint) {
    // remove body from grid since it does not track transforms
    int ibody = constraint->ibody1;
    auto grid = rigidbody->molecule.get_grid();
    {   // backup body and parameters for undo
        auto& body = rigidbody->molecule.get_body(ibody);
        grid->remove(body);
        bodybackup.clear();
        bodybackup.emplace_back(std::move(body), ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);
    }

    // compute new absolute transform parameters for the body
    Vector3<double> pivot = constraint->get_atom2().coordinates();
    auto& body_params = rigidbody->conformation->absolute_parameters.parameters[constraint->ibody1];
    if (par.rotation.has_value()) {
        auto dR = matrix::rotation_matrix(par.rotation.value());
        body_params.transform(pivot, dR);
    } if (par.translation.has_value()) {
        body_params.transform(par.translation.value());
    }

    // reconstruct body from initial conformation using absolute parameters
    auto& body = rigidbody->molecule.get_body(ibody);
    if (par.rotation.has_value() || par.translation.has_value()) {
        body = rigidbody->conformation->initial_conformation[constraint->ibody1];
        rotate_and_translate(matrix::rotation_matrix(body_params.rotation), body_params.translation, body.get_cm(), body);
    } else { // no transformation, so just restore the original conformation
        body = std::move(bodybackup.front().body.value());
        bodybackup.front().body.reset();
    }

    // apply symmetry parameters
    if (par.symmetry_pars.has_value()) {
        body_params.symmetry_pars = add_symmetries(body_params.symmetry_pars, par.symmetry_pars.value());
        apply_symmetry(body_params.symmetry_pars, body);
    }

    // re-add body and refresh grid
    grid->add(body);
    rigidbody->refresh_grid();
}