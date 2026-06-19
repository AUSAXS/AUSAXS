// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/TetrahedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <numbers>

using namespace ausaxs::symmetry;

// 3-fold body-diagonal rotation + 2-fold face rotation generate the rotation group A4
std::unique_ptr<ISymmetry> TetrahedralSymmetry::clone() const {return std::make_unique<TetrahedralSymmetry>(*this);}
const IPolyhedralSymmetry::GroupData& TetrahedralSymmetry::group() const {
    static const GroupData data = build({
        matrix::rotation_matrix<double>({1, 1, 1}, 2*std::numbers::pi/3), 
        matrix::rotation_matrix<double>({0, 0, 1}, std::numbers::pi)}
    , 12);
    return data;
}