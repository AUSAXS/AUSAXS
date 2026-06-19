// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/OctahedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <numbers>

using namespace ausaxs::symmetry;

// 4-fold face rotation + 3-fold body-diagonal rotation generate the rotation group S4
std::unique_ptr<ISymmetry> OctahedralSymmetry::clone() const {return std::make_unique<OctahedralSymmetry>(*this);}
const IPolyhedralSymmetry::GroupData& OctahedralSymmetry::group() const {
    static const GroupData data = build({
        matrix::rotation_matrix<double>({0, 0, 1}, std::numbers::pi/2), 
        matrix::rotation_matrix<double>({1, 1, 1}, 2*std::numbers::pi/3)}
    , 24);
    return data;
}