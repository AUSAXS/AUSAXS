// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/BodyTransformParametersRelative.h>

using namespace ausaxs::rigidbody::parameter;

BodyTransformParametersRelative::BodyTransformParametersRelative() = default;
BodyTransformParametersRelative::BodyTransformParametersRelative(
    const Vector3<double>& translation, 
    const Vector3<double>& rotation, 
    std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetry_pars
)
    : translation(translation), rotation(rotation), symmetry_pars(std::move(symmetry_pars))
{}
bool BodyTransformParametersRelative::operator==(const BodyTransformParametersRelative&) const = default;