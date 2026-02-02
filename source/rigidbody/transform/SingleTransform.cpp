// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/BodyTransformParameters.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <data/Body.h>

using namespace ausaxs::rigidbody::transform;

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::BodyTransformParameters&& par, constraints::DistanceConstraint& constraint) {
    unsigned int ibody = constraint.ibody1;
    auto& body_ref = rigidbody->molecule.get_body(ibody);
    auto grid = rigidbody->molecule.get_grid();

    // Backup current state
    TransformGroup group({&body_ref}, {ibody}, constraint, constraint.get_atom2().coordinates());
    backup(group);

    // Remove old body from grid
    grid->remove(body_ref);

    // Get fresh body from original_conformation and apply absolute transformation
    auto body = rigidbody->conformation->original_conformation[ibody];

    body.rotate(matrix::rotation_matrix(par.rotation));
    body.translate(par.translation);
    symmetry(std::move(par.symmetry_pars), body);

    // Update molecule with transformed body
    rigidbody->molecule.get_body(ibody) = std::move(body);

    // Update configuration with new absolute parameters
    rigidbody->conformation->configuration.parameters[ibody] = std::move(par);

    // Ensure grid has space
    rigidbody->refresh_grid();

    // Add body back to grid
    grid->add(rigidbody->molecule.get_body(ibody));
}