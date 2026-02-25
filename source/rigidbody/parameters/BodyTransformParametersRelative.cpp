// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/BodyTransformParametersRelative.h>

using namespace ausaxs::rigidbody::parameter;

BodyTransformParametersRelative::BodyTransformParametersRelative() = default;
BodyTransformParametersRelative::BodyTransformParametersRelative(const BodyTransformParametersRelative& other) 
    : translation(other.translation), rotation(other.rotation) 
{
    if (!other.symmetry_pars.has_value()) {return;}
    symmetry_pars.emplace();
    symmetry_pars->reserve(other.symmetry_pars->size());
    for (const auto& sym : other.symmetry_pars.value()) {
        symmetry_pars->emplace_back(sym->clone());
    }
}
BodyTransformParametersRelative::BodyTransformParametersRelative(BodyTransformParametersRelative&&) noexcept = default;
BodyTransformParametersRelative::~BodyTransformParametersRelative() = default;

BodyTransformParametersRelative& BodyTransformParametersRelative::operator=(BodyTransformParametersRelative&&) noexcept = default;
BodyTransformParametersRelative& BodyTransformParametersRelative::operator=(const BodyTransformParametersRelative& other) {
    return *this = BodyTransformParametersRelative(other);
}

BodyTransformParametersRelative::BodyTransformParametersRelative(
    const Vector3<double>& translation, 
    const Vector3<double>& rotation, 
    std::vector<std::unique_ptr<symmetry::ISymmetry>> symmetry_pars
)
    : translation(translation), rotation(rotation), symmetry_pars(std::move(symmetry_pars))
{}
bool BodyTransformParametersRelative::operator==(const BodyTransformParametersRelative&) const = default;