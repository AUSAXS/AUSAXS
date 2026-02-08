// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/Parameter.h>
#include <rigidbody/RigidBody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

#include <vector>

using namespace ausaxs::rigidbody::transform;

TransformStrategy::TransformStrategy(observer_ptr<RigidBody> rigidbody) : rigidbody(rigidbody) {}

TransformStrategy::~TransformStrategy() = default;

void TransformStrategy::rotate(const Matrix<double>& M, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (data::Body* body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M]     (data::Body* body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group] (data::Body* body) {body->translate(group.pivot);});
}

void TransformStrategy::translate(const Vector3<double>& t, TransformGroup& group) {
    std::for_each(group.bodies.begin(), group.bodies.end(), [&t] (data::Body* body) {body->translate(t);});
}

void TransformStrategy::rotate_and_translate(const Matrix<double>& M, const Vector3<double>& t, TransformGroup& group) {
    // Check if this is effectively a no-op to avoid floating point drift
    bool is_identity_rotation = true;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(M.index(i, j) - expected) > 1e-10) {
                is_identity_rotation = false;
                break;
            }
        }
        if (!is_identity_rotation) break;
    }
    
    bool is_zero_translation = t.magnitude() < 1e-10;
    
    // If both rotation and translation are effectively zero, skip the transformation
    if (is_identity_rotation && is_zero_translation) {
        return;
    }
    
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group]     (data::Body* body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M]         (data::Body* body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group, &t] (data::Body* body) {body->translate(group.pivot+t);});
}

void TransformStrategy::symmetry(std::vector<parameter::Parameter::SymmetryParameter>&& symmetry_pars, data::Body& body) {
    assert(symmetry_pars.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        body.symmetry().get(i).initial_relation.translation += symmetry_pars[i].rotation_cm;
        body.symmetry().get(i).initial_relation.orientation += symmetry_pars[i].rotation_angle;
        body.symmetry().get(i).repeat_relation.translate += symmetry_pars[i].translation;
    }
}

void TransformStrategy::apply(parameter::Parameter&& par, unsigned int ibody) {
    auto& body = rigidbody->get_body(ibody);

    bodybackup.clear();
    bodybackup.emplace_back(body, ibody);

    // Check if this is effectively a no-op to avoid floating point drift
    bool is_zero_rotation = par.rotation.magnitude() < 1e-10;
    bool is_zero_translation = par.translation.magnitude() < 1e-10;
    
    if (is_zero_rotation && is_zero_translation) {
        return;
    }

    auto grid = rigidbody->get_grid();
    grid->remove(body);

    // translate & rotate
    auto cm = body.get_cm();
    body.translate(-cm);
    body.rotate(matrix::rotation_matrix(par.rotation));
    body.translate(cm + par.translation);

    // update symmetry parameters
    symmetry(std::move(par.symmetry_pars), body);

    grid->add(body);
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        rigidbody->get_body(body.index) = std::move(body.body);
    }
    bodybackup.clear();
}

void TransformStrategy::backup(TransformGroup& group) {
    bodybackup.clear();
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i]);
    }
}
