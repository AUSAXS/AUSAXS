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

void SingleTransform::apply(
    parameter::BodyTransformParametersRelative&& par, observer_ptr<const constraints::IDistanceConstraint> constraint, unsigned int isymmetry_body
) {
    // remove body from grid since it does not track transforms
    int ibody = constraint->ibody1;
    auto grid = rigidbody->molecule.get_grid();
    {   // backup body and parameters for undo
        auto& body = rigidbody->molecule.get_body(ibody);
        grid->remove(body);
        bodybackup.clear();
        bodybackup.emplace_back(std::move(body), ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);
    }

    // the symmetry deltas were generated from isymmetry_body's own symmetry list, so that is the only body they can be applied to; it needs its own backup
    // entry and grid round-trip whenever it is not the transformed body itself (the grid tracks the symmetry copies of every body)
    bool symmetry_other_body = par.symmetry_pars.has_value() && isymmetry_body != static_cast<unsigned int>(ibody);
    if (symmetry_other_body) {
        bodybackup.emplace_back(
            rigidbody->molecule.get_body(isymmetry_body), isymmetry_body, rigidbody->conformation->absolute_parameters.parameters[isymmetry_body]
        );
        grid->remove(rigidbody->molecule.get_body(isymmetry_body));
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
        restore_symmetry(ibody); // rebuilding from the initial conformation also reset the symmetries, so put the accumulated values back
    } else { // no transformation, so just restore the original conformation
        body = std::move(bodybackup.front().body.value());
        bodybackup.front().body.reset();
    }

    // apply the symmetry deltas to the body they were generated for
    if (par.symmetry_pars.has_value()) {apply_symmetry_delta(isymmetry_body, par.symmetry_pars.value());}

    // re-add body and refresh grid
    rigidbody->refresh_grid();
    rigidbody->molecule.get_grid()->add(body); // refresh_grid may reallocate the grid, so re-fetch the pointer
    if (symmetry_other_body) {rigidbody->molecule.get_grid()->add(rigidbody->molecule.get_body(isymmetry_body));}
}