// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/parameters/BodyTransformParametersAbsolute.h>

using namespace ausaxs::rigidbody::parameter;

BodyTransformParametersAbsolute::BodyTransformParametersAbsolute() : translation{0, 0, 0}, rotation{0, 0, 0} {}
BodyTransformParametersAbsolute::BodyTransformParametersAbsolute(const BodyTransformParametersAbsolute& other) 
    : translation(other.translation), rotation(other.rotation) 
{
    symmetry_pars.reserve(other.symmetry_pars.size());
    for (const auto& sym : other.symmetry_pars) {
        symmetry_pars.push_back(sym->clone());
    }
}
BodyTransformParametersAbsolute::BodyTransformParametersAbsolute(BodyTransformParametersAbsolute&&) noexcept = default;
BodyTransformParametersAbsolute::~BodyTransformParametersAbsolute() = default;

BodyTransformParametersAbsolute& BodyTransformParametersAbsolute::operator=(BodyTransformParametersAbsolute&&) noexcept = default;
BodyTransformParametersAbsolute& BodyTransformParametersAbsolute::operator=(const BodyTransformParametersAbsolute& other) {
    return *this = BodyTransformParametersAbsolute(other);
}

BodyTransformParametersAbsolute::BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles) : translation(std::move(dr)), rotation(std::move(euler_angles)) {}

BodyTransformParametersAbsolute::BodyTransformParametersAbsolute(Vector3<double> dr, Vector3<double> euler_angles, std::vector<std::unique_ptr<symmetry::ISymmetry>>&& symmetry_pars)
    : translation(std::move(dr)), rotation(std::move(euler_angles)), symmetry_pars(std::move(symmetry_pars))
{}

void BodyTransformParametersAbsolute::transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation) {
    translation = relative_rotation*(translation - pivot) + pivot;
    rotation = matrix::euler_angles(relative_rotation*matrix::rotation_matrix(rotation));
}

void BodyTransformParametersAbsolute::transform(const Vector3<double>& relative_translation) {
    translation += relative_translation;
}

void BodyTransformParametersAbsolute::transform(const Vector3<double>& pivot, const Matrix<double>& relative_rotation, const Vector3<double>& relative_translation) {
    transform(pivot, relative_rotation);
    transform(relative_translation);
}

bool BodyTransformParametersAbsolute::operator==(const BodyTransformParametersAbsolute&) const = default;