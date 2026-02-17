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

void TransformStrategy::rotate_and_translate(const Matrix<double>& M, const Vector3<double>& t, const Vector3<double>& pivot, data::Body& body) {
    body.translate(-pivot);
    body.rotate(M);
    body.translate(pivot+t);
}

void TransformStrategy::apply_symmetry(const std::vector<symmetry::Symmetry>& symmetry, data::Body& body) {
    assert(symmetry.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        auto& current_sym = body.symmetry().get(i);
        current_sym.initial_relation.translation = symmetry[i].initial_relation.translation;
        current_sym.initial_relation.rotation = symmetry[i].initial_relation.rotation;
        body.get_signaller()->modified_symmetry(i);
    }
}

std::vector<ausaxs::symmetry::Symmetry> TransformStrategy::add_symmetries(const std::vector<symmetry::Symmetry>& current, const std::vector<symmetry::Symmetry>& delta) {
    assert(current.size() == delta.size());
    std::vector<symmetry::Symmetry> result = current;
    for (int i = 0; i < static_cast<int>(current.size()); ++i) {
        result[i].initial_relation.translation += delta[i].initial_relation.translation;
        result[i].initial_relation.rotation += delta[i].initial_relation.rotation;
    }
    return result;
}

void TransformStrategy::apply(parameter::BodyTransformParametersRelative&& par, unsigned int ibody) {
    auto& body = rigidbody->molecule.get_body(ibody);
    bodybackup.clear();
    bodybackup.emplace_back(body, ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);

    // remove body from grid since it does not track transforms
    auto grid = rigidbody->molecule.get_grid();
    grid->remove(body);

    // compute new absolute transform parameters for the body
    parameter::BodyTransformParametersAbsolute& body_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
    if (par.rotation.has_value()) {body_params.rotation = matrix::euler_angles(matrix::rotation_matrix(body_params.rotation)*matrix::rotation_matrix(par.rotation.value()));}
    if (par.translation.has_value()) {body_params.translation += par.translation.value();}

    // apply transformations
    body = rigidbody->conformation->initial_conformation[ibody];
    if (par.rotation.has_value() || par.translation.has_value()) {
        rotate_and_translate(matrix::rotation_matrix(body_params.rotation), body_params.translation, body.get_cm(), body);
    }

    // update and apply symmetry parameters
    if (par.symmetry_pars.has_value()) {
        body_params.symmetry_pars = add_symmetries(body_params.symmetry_pars, par.symmetry_pars.value());
        apply_symmetry(body_params.symmetry_pars, body);
    }

    // re-add body and refresh grid
    grid->add(body);
    rigidbody->refresh_grid();
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
    for (int i = 0; i < static_cast<int>(group.bodies.size()); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i], rigidbody->conformation->absolute_parameters.parameters[group.indices[i]]);
    }
}