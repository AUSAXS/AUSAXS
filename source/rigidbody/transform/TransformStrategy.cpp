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
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/PointSymmetry.h>
#include <math/MatrixUtils.h>

#include <vector>

using namespace ausaxs;
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

void TransformStrategy::apply_symmetry(const std::vector<std::unique_ptr<symmetry::ISymmetry>>& symmetry, data::Body& body) {
    assert(symmetry.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        auto current_sym = body.symmetry().get(i);
        auto src_t = symmetry[i]->span_translation();
        auto src_r = symmetry[i]->span_rotation();
        auto dst_t = current_sym->span_translation();
        auto dst_r = current_sym->span_rotation();

        assert(src_t.size() == dst_t.size() && src_r.size() == dst_r.size() && "TransformStrategy::apply_symmetry: Symmetry parameter size mismatch.");
        std::copy(src_t.begin(), src_t.end(), dst_t.begin());
        std::copy(src_r.begin(), src_r.end(), dst_r.begin());
        body.get_signaller()->modified_symmetry(i);
    }
}

void TransformStrategy::add_symmetries(
    std::vector<std::unique_ptr<symmetry::ISymmetry>>& current, const std::vector<std::unique_ptr<symmetry::ISymmetry>>& delta
) {
    // sanity checks to ensure compatible symmetries are being added
    assert(current.size() == delta.size() && "TransformStrategy::add_symmetries: Symmetry parameter size mismatch.");
    assert([&]() -> bool {
        for (unsigned int i = 0; i < current.size(); ++i) {
            assert(current[i] != nullptr && "TransformStrategy::add_symmetries: Current symmetry parameter cannot be null.");
            assert(delta[i] != nullptr && "TransformStrategy::add_symmetries: Delta symmetry parameter cannot be null.");
            if (       dynamic_cast<symmetry::CyclicSymmetry*>(current[i].get())) {
                assert(dynamic_cast<symmetry::CyclicSymmetry*>(delta[i].get()) && "TransformStrategy::add_symmetries: Symmetry type mismatch.");
            } else if (dynamic_cast<symmetry::PointSymmetry*>(current[i].get())) {
                assert(dynamic_cast<symmetry::PointSymmetry*>(delta[i].get()) && "TransformStrategy::add_symmetries: Symmetry type mismatch.");
            } else {
                assert(false && "TransformStrategy::add_symmetries: Unchecked symmetry type.");
            }
        }
        return true;
    }());

    for (unsigned int i = 0; i < current.size(); ++i) {current[i]->add(delta[i].get());}
}

void TransformStrategy::apply(parameter::BodyTransformParametersRelative&& par, unsigned int ibody) {
    assert(ibody < rigidbody->molecule.size_body() && "TransformStrategy::apply: Body index out of range.");
    assert(
        (par.rotation.has_value() || par.translation.has_value() || par.symmetry_pars.has_value()) 
        && "TransformStrategy::apply: No transformation specified."
    );

    // remove body from grid since it does not track transforms
    auto grid = rigidbody->molecule.get_grid();
    {   // backup body and parameters for undo
        auto& body = rigidbody->molecule.get_body(ibody);
        grid->remove(body);

        bodybackup.clear();
        bodybackup.emplace_back(std::move(body), ibody, rigidbody->conformation->absolute_parameters.parameters[ibody]);
    }

    // compute new absolute transform parameters for the body
    parameter::BodyTransformParametersAbsolute& body_params = rigidbody->conformation->absolute_parameters.parameters[ibody];
    if (par.rotation.has_value()) {body_params.rotation = matrix::euler_angles(matrix::rotation_matrix(body_params.rotation)*matrix::rotation_matrix(par.rotation.value()));}
    if (par.translation.has_value()) {body_params.translation += par.translation.value();}

    // apply transformations
    auto& body = rigidbody->molecule.get_body(ibody);
    if (par.rotation.has_value() || par.translation.has_value()) {
        body = rigidbody->conformation->initial_conformation[ibody];
        rotate_and_translate(matrix::rotation_matrix(body_params.rotation), body_params.translation, body.get_cm(), body);
    } else { // no transformation, so just restore the original conformation
        body = std::move(bodybackup.front().body.value());
        bodybackup.front().body.reset();
    }

    // update and apply symmetry parameters
    if (par.symmetry_pars.has_value()) {
        add_symmetries(body_params.symmetry_pars, par.symmetry_pars.value());
        apply_symmetry(body_params.symmetry_pars, body);
    }

    // re-add body and refresh grid
    rigidbody->refresh_grid();
    grid->add(body);
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        if (body.body.has_value()) {rigidbody->molecule.get_body(body.index) = std::move(body.body.value());}
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