// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PointSymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

PointSymmetry::PointSymmetry() = default;

PointSymmetry::PointSymmetry(const Vector3<double>& offset, const Vector3<double>& orient_axis) {
    translation = offset;
    rotation  = orient_axis;
}

bool PointSymmetry::is_closed() const { return false; }

std::unique_ptr<ISymmetry> PointSymmetry::clone() const {
    return std::make_unique<PointSymmetry>(*this);
}

std::function<ausaxs::Vector3<double>(ausaxs::Vector3<double>)> PointSymmetry::get_transform(const Vector3<double>& cm, int rep) const {
    assert(rep <= 1 && "PointSymmetry always generates exactly one copy (rep must be 1).");

    auto normed_axis = rotation.magnitude() > 1e-9 ? rotation / rotation.magnitude() : Vector3<double>{0, 0, 1};
    auto R = matrix::rotation_matrix<double>(normed_axis);

    // final transform is p' = R*(p - cm) + cm + d
    //                       = R*p + (cm + d - R*cm)
    auto T = cm + translation - R*cm;
    return [R=std::move(R), T=std::move(T)](Vector3<double> v) {
        return R * v + T;
    };
}

unsigned int PointSymmetry::repetitions() const {return 1;}
std::span<double> PointSymmetry::span_translation() {return std::span<double>(translation.begin(), translation.end());}
std::span<double> PointSymmetry::span_rotation() {return std::span<double>(rotation.begin(), rotation.end());}

ISymmetry& PointSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const PointSymmetry*>(other);
    assert(cast != nullptr && "Can only add PointSymmetry with another PointSymmetry.");
    this->translation += cast->translation;
    this->rotation += cast->rotation;
    return *this;
}