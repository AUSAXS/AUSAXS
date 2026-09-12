// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/transform/TransformGroup.h>

#include <data/Body.h>
#include <rigidbody/constraints/IDistanceConstraint.h>

using namespace ausaxs::rigidbody::transform;

TransformGroup::TransformGroup(
    std::vector<observer_ptr<data::Body>> bodies, std::vector<int> indices, 
    observer_ptr<const constraints::IDistanceConstraint> target, Vector3<double> pivot
)
    : bodies(std::move(bodies)), indices(std::move(indices)), target(target), pivot(std::move(pivot)) {}

TransformGroup::~TransformGroup() = default;