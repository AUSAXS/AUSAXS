// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/SingleTransform.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/RelativeTransformParameters.h>
#include <rigidbody/constraints/DistanceConstraint.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <grid/Grid.h>
#include <data/Body.h>
#include <math/MatrixUtils.h>

using namespace ausaxs::rigidbody::transform;

namespace {
    // Helper to compose rotation matrices and extract Euler angles (XYZ extrinsic convention)
    ausaxs::Vector3<double> compose_rotation(const ausaxs::Matrix<double>& R_delta, const ausaxs::Vector3<double>& current_rotation) {
        auto R_current = ausaxs::matrix::rotation_matrix(current_rotation);
        auto R_new = R_delta * R_current;
        
        // Extract Euler angles (alpha, beta, gamma) from the rotation matrix
        double beta = std::asin(std::clamp(-R_new(2, 0), -1.0, 1.0));
        double cos_beta = std::cos(beta);
        
        double alpha, gamma;
        if (std::abs(cos_beta) > 1e-10) {
            alpha = std::atan2(R_new(2, 1), R_new(2, 2));
            gamma = std::atan2(R_new(1, 0), R_new(0, 0));
        } else {
            // Gimbal lock
            alpha = 0.0;
            gamma = std::atan2(-R_new(0, 1), R_new(1, 1));
        }
        
        return {alpha, beta, gamma};
    }
}

SingleTransform::SingleTransform(observer_ptr<Rigidbody> rigidbody) : TransformStrategy(rigidbody) {}

SingleTransform::~SingleTransform() = default;

void SingleTransform::apply(parameter::RelativeTransformParameters&& par, constraints::DistanceConstraint& constraint) {
    unsigned int ibody = constraint.ibody1;
    auto& body_ref = rigidbody->molecule.get_body(ibody);
    auto grid = rigidbody->molecule.get_grid();

    // Backup current state
    // Pivot is the position of the constraining atom (atom2 is in the non-moving body)
    Vector3<double> pivot = constraint.get_atom2().coordinates();
    TransformGroup group({&body_ref}, {ibody}, constraint, pivot);
    backup(group);

    // Remove old body from grid
    grid->remove(body_ref);

    // Compute new absolute parameters from current + delta
    // Rotation happens around the pivot point
    auto& current_params = rigidbody->conformation->configuration.parameters[ibody];
    auto R_delta = matrix::rotation_matrix(par.rotation);
    
    auto new_rotation = compose_rotation(R_delta, current_params.rotation);
    // t_new = R_delta * (t_old - pivot) + pivot + t_delta
    auto new_translation = R_delta * (current_params.translation - pivot) + pivot + par.translation;

    // Get fresh body from original_conformation and apply new absolute transformation
    auto body = rigidbody->conformation->original_conformation[ibody];

    body.rotate(matrix::rotation_matrix(new_rotation));
    body.translate(new_translation);
    symmetry(std::move(par.symmetry_pars), body);

    // Update molecule with transformed body
    rigidbody->molecule.get_body(ibody) = std::move(body);

    // Update configuration with new absolute parameters
    current_params.rotation = new_rotation;
    current_params.translation = new_translation;

    // Ensure grid has space
    rigidbody->refresh_grid();

    // Add body back to grid
    grid->add(rigidbody->molecule.get_body(ibody));
}