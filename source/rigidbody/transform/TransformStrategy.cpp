/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/transform/TransformStrategy.h>
#include <rigidbody/transform/TransformGroup.h>
#include <rigidbody/transform/BackupBody.h>
#include <rigidbody/RigidBody.h>
#include <grid/detail/GridMember.h>
#include <grid/Grid.h>

#include <vector>

using namespace ausaxs::rigidbody::transform;

TransformStrategy::TransformStrategy(RigidBody* rigidbody) : rigidbody(rigidbody) {}

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
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group]     (data::Body* body) {body->translate(-group.pivot);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&M]         (data::Body* body) {body->rotate(M);});
    std::for_each(group.bodies.begin(), group.bodies.end(), [&group, &t] (data::Body* body) {body->translate(group.pivot+t);});
}

void TransformStrategy::apply(const Matrix<double>& M, const Vector3<double>& t, unsigned int ibody) {
    auto& body = rigidbody->get_body(ibody);

    bodybackup.clear();
    bodybackup.emplace_back(body, ibody);

    auto grid = rigidbody->get_grid();
    grid->remove(body);

    auto cm = body.get_cm();
    body.translate(-cm);
    body.rotate(M);
    body.translate(cm + t);

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
