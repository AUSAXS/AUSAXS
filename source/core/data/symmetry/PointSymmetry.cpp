// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PointSymmetry.h>

#include <math/MatrixUtils.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

PointSymmetry::PointSymmetry() = default;

PointSymmetry::PointSymmetry(const Vector3<double>& translation, const Vector3<double>& rotation) : translation(translation), rotation(rotation) {}

bool PointSymmetry::is_closed() const { return false; }
std::string PointSymmetry::type_name() const { return "p2"; }

std::unique_ptr<ISymmetry> PointSymmetry::clone() const {
    return std::make_unique<PointSymmetry>(*this);
}

AffineTransform PointSymmetry::_make_transform(const Vector3<double>& anchor, int rep) const {
    assert(rep <= 1 && "PointSymmetry always generates exactly one copy (rep must be 1).");
    if (rep == 0) {return {};} // identity

    // final transform is v' = R*(v - cm) + cm + d
    //                       = R*v + (cm + d - R*cm)
    auto R = matrix::rotation_matrix<double>(rotation);
    auto T = anchor + translation - R*anchor;
    return {.rotation=std::move(R), .translation=T};
}

int PointSymmetry::repetitions() const {return 1;}
std::span<double> PointSymmetry::span_translation() {return {translation.begin(), translation.end()};}
std::span<double> PointSymmetry::span_rotation() {return {rotation.begin(), rotation.end()};}

ISymmetry& PointSymmetry::add(observer_ptr<const ISymmetry> other) {
    const auto* cast = dynamic_cast<const PointSymmetry*>(other);
    assert(cast != nullptr && "Can only add PointSymmetry with another PointSymmetry.");
    this->translation += cast->translation;
    this->rotation += cast->rotation;
    return *this;
}