// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/PointSymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::symmetry;

PointSymmetry::PointSymmetry() = default;

PointSymmetry::PointSymmetry(const Vector3<double>& offset, const Vector3<double>& orient_axis, double orient_angle) {
    initial_relation.translation = offset;
    repeat_relation.axis  = orient_axis;
    repeat_relation.angle = orient_angle;
    repeat_relation.translation = {0, 0, 0};
}

bool PointSymmetry::is_closed() const { return false; }

std::unique_ptr<ISymmetry> PointSymmetry::clone() const {
    return std::make_unique<PointSymmetry>(*this);
}

std::function<ausaxs::Vector3<double>(ausaxs::Vector3<double>)> PointSymmetry::get_transform(const Vector3<double>& cm, int rep) const {
    assert(rep == 1 && "PointSymmetry always generates exactly one copy (rep must be 1).");

    const auto& d    = initial_relation.translation;
    const auto& axis = repeat_relation.axis;
    const double theta = repeat_relation.angle;

    auto normed_axis = axis.magnitude() > 1e-9 ? axis / axis.magnitude() : Vector3<double>{0, 0, 1};
    auto R = matrix::rotation_matrix<double>(normed_axis, theta);

    // T = (I - R)*cm + d
    auto T = cm - R * cm + d;

    return [R=std::move(R), T=std::move(T)](Vector3<double> v) {
        return R * v + T;
    };
}

ISymmetry& PointSymmetry::add(observer_ptr<const ISymmetry> other) {
    auto cast = dynamic_cast<const PointSymmetry*>(other);
    assert(cast != nullptr && "Can only add PointSymmetry with another PointSymmetry.");
    this->initial_relation.translation += cast->initial_relation.translation;
    this->repeat_relation.axis += cast->repeat_relation.axis;
    this->repeat_relation.angle += cast->repeat_relation.angle;
    return *this;
}