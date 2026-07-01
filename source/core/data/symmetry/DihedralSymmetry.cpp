// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/symmetry/DihedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <cassert>
#include <numbers>

using namespace ausaxs::symmetry;

// C_n about the principal axis (z) + a C_2 about a perpendicular axis (x) generate the order-2n
// dihedral rotation group. The perpendicular flip is essential: about a shared axis the two
// rotations would only close into the cyclic group C_2n.
DihedralSymmetry::DihedralSymmetry(int n)
    : data(build({
        matrix::rotation_matrix<double>({0, 0, 1}, 2*std::numbers::pi/n),
        matrix::rotation_matrix<double>({1, 0, 0}, std::numbers::pi)
    }, 2*n))
{
    assert(2 <= n && "DihedralSymmetry: order n must be at least 2.");
}

std::unique_ptr<ISymmetry> DihedralSymmetry::clone() const {return std::make_unique<DihedralSymmetry>(*this);}
const IPolyhedralSymmetry::GroupData& DihedralSymmetry::group() const {return data;}
