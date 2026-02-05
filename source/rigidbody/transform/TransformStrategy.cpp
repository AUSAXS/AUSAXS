// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/BodyTransformParametersRelative.h>
#include <rigidbody/detail/SystemSpecification.h>
#include <rigidbody/Rigidbody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>
#include <data/state/BoundSignaller.h>
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

void TransformStrategy::apply_symmetry(std::vector<symmetry::Symmetry>&& symmetry, data::Body& body) {
    assert(symmetry.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        auto& current_sym = body.symmetry().get(i);
        current_sym.initial_relation.translation = symmetry[i].initial_relation.translation;
        current_sym.initial_relation.orientation = symmetry[i].initial_relation.orientation;
        body.get_signaller()->modified_symmetry(i);
    }
}

void TransformStrategy::update_body(unsigned int ibody, parameter::BodyTransformParametersAbsolute&& pars) {
    auto& current_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
    auto& body = rigidbody->molecule.get_body(ibody);

    // Check what actually changed to avoid unnecessary recalculations
    bool rotation_changed = (pars.rotation != current_params.rotation);
    bool translation_changed = (pars.translation != current_params.translation);
    bool symmetry_changed = !pars.symmetry_pars.empty();

    // If rotation or translation changed, we need to reconstruct from initial_conformation
    if (rotation_changed || translation_changed) {
        auto fresh_body = rigidbody->conformation->initial_conformation[ibody];
        fresh_body.rotate(matrix::rotation_matrix(pars.rotation));
        fresh_body.translate(pars.translation);
        
        // Apply current symmetry state (will be updated below if symmetry_changed)
        if (!pars.symmetry_pars.empty()) {
            apply_symmetry(std::move(pars.symmetry_pars), fresh_body);
        } else {
            for (unsigned int i = 0; i < current_params.symmetry_pars.size() && i < fresh_body.size_symmetry(); ++i) {
                fresh_body.symmetry().get(i) = current_params.symmetry_pars[i];
                fresh_body.get_signaller()->modified_symmetry(i);
            }
        }
        
        body = std::move(fresh_body);
        current_params.rotation = std::move(pars.rotation);
        current_params.translation = std::move(pars.translation);
    } else if (symmetry_changed) {
        // Only symmetry changed - update in place without reconstruction
        apply_symmetry(std::move(pars.symmetry_pars), body);
    }

    // Extract final symmetry parameters
    if (symmetry_changed) {
        current_params.symmetry_pars.clear();
        for (unsigned int i = 0; i < body.size_symmetry(); ++i) {
            current_params.symmetry_pars.push_back(body.symmetry().get(i));
        }
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

void TransformStrategy::apply(parameter::BodyTransformParametersRelative&& par, unsigned int ibody) {
    auto grid = rigidbody->molecule.get_grid();
    auto& body = rigidbody->molecule.get_body(ibody);

    // Backup and remove old body
    bodybackup.clear();
    bodybackup.emplace_back(body, ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);
    grid->remove(body);

    // Compute new absolute parameters from current + delta
    auto& current_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
    auto R_delta = matrix::rotation_matrix(par.rotation);
    auto new_rotation = compose_rotation(current_params.rotation, par.rotation);
    auto new_translation = R_delta * current_params.translation + par.translation;
    auto new_symmetry = current_params.symmetry_pars;
    for (unsigned int i = 0; i < new_symmetry.size(); ++i) {
        auto& sym = new_symmetry[i];
        sym.initial_relation.orientation += par.symmetry_pars[i].initial_relation.orientation;
        sym.initial_relation.translation += par.symmetry_pars[i].initial_relation.translation;
    }

    // Update body with new absolute parameters
    update_body(ibody, {new_rotation, new_translation, std::move(new_symmetry)});

    // Refresh grid and re-add body
    rigidbody->refresh_grid();
    grid->add(rigidbody->molecule.get_body(ibody));
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        rigidbody->molecule.get_body(body.index) = std::move(body.body);
        rigidbody->conformation->absolute_parameters.parameters[body.index] = std::move(body.params);
    }
    bodybackup.clear();
}

void TransformStrategy::backup(TransformGroup& group) {
    bodybackup.clear();
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i], rigidbody->conformation->absolute_parameters.parameters[group.indices[i]]);
    }
}
