#pragma once

#include <array>
#include <algorithm>

#include "data/Protein.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/TransformationStrategy.h"
#include "rigidbody/RigidTransform.h"

void RigidTransform::rotate(const Vector3& axis, const double rad, Constraint& constraint) {
    std::for_each(protein->constraints.begin(), protein->constraints.end(), [] (const Constraint& constraint) {});
}

void RigidTransform::translate(const Vector3& v, Constraint& constraint) {}