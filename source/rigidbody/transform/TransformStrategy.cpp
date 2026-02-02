// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include "data/symmetry/Symmetry.h"
#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/RelativeTransformParameters.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <math/MatrixUtils.h>

#include <vector>

using namespace ausaxs::rigidbody::transform;

TransformStrategy::TransformStrategy(observer_ptr<Rigidbody> rigidbody) : rigidbody(rigidbody) {}

TransformStrategy::~TransformStrategy() = default;

void TransformStrategy::rotate(const Matrix<double>& M, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (observer_ptr<data::Body> body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M]     (observer_ptr<data::Body> body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (observer_ptr<data::Body> body) {body->translate(group.pivot);});
}

void TransformStrategy::translate(const Vector3<double>& t, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&t] (data::Body* body) {body->translate(t);});
}

void TransformStrategy::rotate_and_translate(const Matrix<double>& M, const Vector3<double>& t, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group]     (observer_ptr<data::Body> body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M]         (observer_ptr<data::Body> body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group, &t] (observer_ptr<data::Body> body) {body->translate(group.pivot+t);});
}

void TransformStrategy::symmetry(std::vector<symmetry::Symmetry>&& symmetry_pars, data::Body& body) {
    assert(symmetry_pars.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        body.symmetry().get(i) = symmetry_pars[i];
    }
}

namespace {
    // Helper to compose rotation matrices and extract Euler angles (XYZ extrinsic convention)
    // The rotation matrix uses the convention: R = Rz(gamma) * Ry(beta) * Rx(alpha)
    // Matrix layout:
    //   [0,0] = cos(beta)*cos(gamma)
    //   [1,0] = cos(beta)*sin(gamma)
    //   [2,0] = -sin(beta)
    //   [2,1] = sin(alpha)*cos(beta)
    //   [2,2] = cos(alpha)*cos(beta)
    ausaxs::Vector3<double> compose_rotation(const ausaxs::Vector3<double>& current_rotation, const ausaxs::Vector3<double>& delta_rotation) {
        auto R_current = ausaxs::matrix::rotation_matrix(current_rotation);
        auto R_delta = ausaxs::matrix::rotation_matrix(delta_rotation);
        auto R_new = R_delta * R_current;
        
        // Extract Euler angles (alpha, beta, gamma) from the rotation matrix
        // Using XYZ extrinsic convention matching rotation_matrix(alpha, beta, gamma)
        double beta = std::asin(std::clamp(-R_new(2, 0), -1.0, 1.0));
        double cos_beta = std::cos(beta);
        
        double alpha, gamma;
        if (std::abs(cos_beta) > 1e-10) {
            // Normal case: cos(beta) != 0
            alpha = std::atan2(R_new(2, 1), R_new(2, 2));
            gamma = std::atan2(R_new(1, 0), R_new(0, 0));
        } else {
            // Gimbal lock: beta = +/- pi/2
            // In this case, alpha and gamma are not uniquely determined
            // Convention: set alpha = 0 and solve for gamma
            alpha = 0.0;
            gamma = std::atan2(-R_new(0, 1), R_new(1, 1));
        }
        
        return {alpha, beta, gamma};
    }
}

void TransformStrategy::apply(parameter::RelativeTransformParameters&& par, unsigned int ibody) {
    auto grid = rigidbody->molecule.get_grid();

    {   // remove old body and backup
        auto& body = rigidbody->molecule.get_body(ibody);

        bodybackup.clear();
        bodybackup.emplace_back(body, ibody, rigidbody->conformation->configuration.parameters[ibody]);

        grid->remove(body);
    }

    {   // Compute new absolute parameters from current + delta
        auto& current_params = rigidbody->conformation->configuration.parameters[ibody];
        auto R_delta = matrix::rotation_matrix(par.rotation);
        
        // R_new = R_delta * R_current, t_new = R_delta * t_current + t_delta
        auto new_rotation = compose_rotation(current_params.rotation, par.rotation);
        auto new_translation = R_delta * current_params.translation + par.translation;

        // Get fresh body and apply the new absolute transformation
        assert(ibody < rigidbody->conformation->original_conformation.size() && "ibody out of bounds");
        auto body = rigidbody->conformation->original_conformation[ibody];

        body.rotate(matrix::rotation_matrix(new_rotation));
        body.translate(new_translation);

        // Update symmetry parameters
        symmetry(std::move(par.symmetry_pars), body);

        // Update the molecule with the transformed body
        rigidbody->molecule.get_body(ibody) = std::move(body);

        // Update configuration with new absolute parameters
        current_params.rotation = new_rotation;
        current_params.translation = new_translation;

        // Ensure there is space for the new conformation in the grid
        rigidbody->refresh_grid();
    }

    // Re-add the body to the grid
    grid->add(rigidbody->molecule.get_body(ibody));
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        rigidbody->molecule.get_body(body.index) = std::move(body.body);
        rigidbody->conformation->configuration.parameters[body.index] = std::move(body.params);
    }
    bodybackup.clear();
}

void TransformStrategy::backup(TransformGroup& group) {
    bodybackup.clear();
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i], rigidbody->conformation->configuration.parameters[group.indices[i]]);
    }
}
