// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/parameters/Parameter.h>
#include <rigidbody/detail/Conformation.h>
#include <rigidbody/Rigidbody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

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

void TransformStrategy::symmetry(std::vector<parameter::Parameter::SymmetryParameter>&& symmetry_pars, data::Body& body) {
    assert(symmetry_pars.size() == body.size_symmetry());
    for (int i = 0; i < static_cast<int>(body.size_symmetry()); ++i) {
        body.symmetry().get(i).initial_relation.translation += symmetry_pars[i].rotation_cm;
        body.symmetry().get(i).initial_relation.orientation += symmetry_pars[i].rotation_angle;
        body.symmetry().get(i).repeat_relation.translate += symmetry_pars[i].translation;
    }
}

void TransformStrategy::apply(parameter::Parameter&& par, unsigned int ibody) {
    auto grid = rigidbody->molecule.get_grid();

    {   // remove old body
        auto& body = rigidbody->molecule.get_body(ibody);

        bodybackup.clear();
        bodybackup.emplace_back(body, ibody); //! std::move

        grid->remove(body);
    }

    {   // get fresh body and apply absolute transformation
        assert(ibody < rigidbody->conformation->original_conformation.size() && "ibody out of bounds");
        auto body = rigidbody->conformation->original_conformation[ibody];

        // translate & rotate
        //! the first two ops can be optimized away by moving all original bodies to origo, and using their original cm as starting points
        auto cm = body.get_cm();
        body.translate(-cm);
        body.rotate(matrix::rotation_matrix(par.rotation));
        body.translate(cm + par.translation);

        // update symmetry parameters
        symmetry(std::move(par.symmetry_pars), body);

        // update the conformation
        rigidbody->molecule.get_body(ibody) = std::move(body);

        // ensure there is space for the new conformation in the grid
        rigidbody->refresh_grid();
    }

    // finally, re-add the body to the grid
    grid->add(rigidbody->molecule.get_body(ibody));
}

void TransformStrategy::undo() {
    for (auto& body : bodybackup) {
        rigidbody->molecule.get_body(body.index) = std::move(body.body);
    }
    bodybackup.clear();
}

void TransformStrategy::backup(TransformGroup& group) {
    bodybackup.clear();
    for (unsigned int i = 0; i < group.bodies.size(); i++) {
        bodybackup.emplace_back(*group.bodies[i], group.indices[i]);
    }
}
