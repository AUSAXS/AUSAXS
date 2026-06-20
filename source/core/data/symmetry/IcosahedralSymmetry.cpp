// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/IcosahedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <numbers>

using namespace ausaxs::symmetry;

// 5-fold vertex + 2-fold edge + 3-fold face rotations of an icosahedron with vertices at the
// cyclic permutations of (0, +-1, +-phi) generate the group A5
std::unique_ptr<ISymmetry> IcosahedralSymmetry::clone() const {return std::make_unique<IcosahedralSymmetry>(*this);}
const IPolyhedralSymmetry::GroupData& IcosahedralSymmetry::group() const {
    static const double phi = (1 + std::sqrt(5.0))/2;
    static const GroupData data = build({
        matrix::rotation_matrix<double>({0, 1, phi}, 2*std::numbers::pi/5), 
        matrix::rotation_matrix<double>({1, 0, 0}, std::numbers::pi), 
        matrix::rotation_matrix<double>({1, 1, 1}, 2*std::numbers::pi/3)}
    , 60);
    return data;
}