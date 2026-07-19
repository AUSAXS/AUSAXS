// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SymmetryFit.h>
#include <rigidbody/sequencer/detail/Superposition.h>
#include <rigidbody/sequencer/detail/SymmetricEigen.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/IPolyhedralSymmetry.h>
#include <math/MatrixUtils.h>

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace ausaxs;
namespace rbsd = ausaxs::rigidbody::sequencer::detail;

namespace {
    double det3(const Matrix<double>& A) {
        return A(0,0)*(A(1,1)*A(2,2) - A(1,2)*A(2,1))
             - A(0,1)*(A(1,0)*A(2,2) - A(1,2)*A(2,0))
             + A(0,2)*(A(1,0)*A(2,1) - A(1,1)*A(2,0));
    }

    // Inverse of a 3x3 matrix via the adjugate. Assumes non-singular input.
    Matrix<double> inverse3(const Matrix<double>& A) {
        double a = A(0,0), b = A(0,1), c = A(0,2);
        double d = A(1,0), e = A(1,1), f = A(1,2);
        double g = A(2,0), h = A(2,1), i = A(2,2);
        double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
        Matrix<double> inv(3, 3);
        inv(0,0) =  (e*i - f*h)/det; inv(0,1) = -(b*i - c*h)/det; inv(0,2) =  (b*f - c*e)/det;
        inv(1,0) = -(d*i - f*g)/det; inv(1,1) =  (a*i - c*g)/det; inv(1,2) = -(a*f - c*d)/det;
        inv(2,0) =  (d*h - e*g)/det; inv(2,1) = -(a*h - b*g)/det; inv(2,2) =  (a*e - b*d)/det;
        return inv;
    }

    // Recover the (unit) rotation axis and angle of a rotation matrix.
    std::pair<Vector3<double>, double> axis_angle(const Matrix<double>& R) {
        double trace = R(0,0) + R(1,1) + R(2,2);
        double cos_t = std::clamp((trace - 1)/2, -1.0, 1.0);
        double angle = std::acos(cos_t);

        Vector3<double> axis{R(2,1) - R(1,2), R(0,2) - R(2,0), R(1,0) - R(0,1)};
        double s = axis.magnitude();
        if (1e-8 < s) {return {axis/s, angle};}       // generic case

        if (0 < cos_t) {return {{0, 0, 1}, 0.0};}      // near-identity: axis undefined

        // near-pi rotation: R is symmetric, so recover the axis from the diagonal of (R+I)/2 = a a^T
        Vector3<double> a{
            std::sqrt(std::max(0.0, (R(0,0) + 1)/2)),
            std::sqrt(std::max(0.0, (R(1,1) + 1)/2)),
            std::sqrt(std::max(0.0, (R(2,2) + 1)/2))
        };
        // fix relative signs against the off-diagonal entries, anchoring on the largest component
        if (a.x() >= a.y() && a.x() >= a.z()) {
            a = {a.x(), std::copysign(a.y(), R(0,1)), std::copysign(a.z(), R(0,2))};
        } else if (a.y() >= a.z()) {
            a = {std::copysign(a.x(), R(0,1)), a.y(), std::copysign(a.z(), R(1,2))};
        } else {
            a = {std::copysign(a.x(), R(0,2)), std::copysign(a.y(), R(1,2)), a.z()};
        }
        return {a.normalize(), angle};
    }

    // Solve the shared rotation centre c from the per-copy transforms: every point-group copy k
    // obeys t_k = (I - R_k) c. Stacking all copies and solving in the least-squares sense averages
    // out small asymmetries. A tiny regularisation towards cm pins the otherwise-free component along
    // a rotation axis (which lies in the null space of every I - R_k for a cyclic group).
    Vector3<double> solve_common_centre(
        const std::vector<Matrix<double>>& R, const std::vector<Vector3<double>>& t, const Vector3<double>& cm
    ) {
        Matrix<double> normal(3, 3);
        Vector3<double> rhs{0, 0, 0};
        for (std::size_t k = 0; k < R.size(); ++k) {
            Matrix<double> M = Matrix<double>::identity(3) - R[k];
            normal += M.transpose()*M;
            rhs += M.transpose()*t[k];
        }
        double lambda = 1e-9*(normal(0,0) + normal(1,1) + normal(2,2))/3 + 1e-12;
        for (unsigned int d = 0; d < 3; ++d) {normal(d, d) += lambda;}
        rhs += lambda*cm;
        return inverse3(normal)*rhs;
    }

    // Recover the group frame F from the measured copy rotations R[k] = F G[k] F^T (correspondence
    // R[k] <-> G[k] follows the copy ordering). An initial estimate is the null vector of the stacked
    // linear operator L (from R[k] F - F G[k] = 0). For groups whose 3D rotation representation is
    // irreducible (tetra/octa/icosa) that null space is one-dimensional and the estimate is exact;
    // for the reducible dihedral groups it is degenerate, so the raw estimate need not be a valid
    // frame. In all cases the estimate is refined by a conjugation-averaging fixed point,
    //     F <- nearest_rotation( sum_k R[k] F G[k]^T ),
    // whose fixed points are exactly the valid frames and which, going through optimal_rotation, only
    // ever visits proper rotations. Any valid frame reproduces the assembly identically (distinct
    // frames differ only by a symmetry of the group), so a discrete frame ambiguity is harmless.
    Matrix<double> solve_frame(const std::vector<Matrix<double>>& R, const std::vector<Matrix<double>>& G) {
        assert(R.size() == G.size() && "solve_frame: rotation/element count mismatch.");
        Matrix<double> Q(9, 9);
        for (std::size_t k = 1; k < R.size(); ++k) {         // element 0 is the identity: L is zero
            Matrix<double> L(9, 9);
            for (unsigned int a = 0; a < 3; ++a) {
                for (unsigned int b = 0; b < 3; ++b) {
                    unsigned int row = 3*a + b;
                    for (unsigned int i = 0; i < 3; ++i) {L(row, 3*i + b) += R[k](a, i);}   //  sum_i R(a,i) F(i,b)
                    for (unsigned int j = 0; j < 3; ++j) {L(row, 3*a + j) -= G[k](j, b);}   // -sum_j F(a,j) G(j,b)
                }
            }
            Q += L.transpose()*L;
        }

        auto eig = rbsd::symmetric_eigen(Q);
        const auto& vecF = eig.vectors.front();          // smallest eigenvalue -> initial null-vector estimate
        Matrix<double> F(3, 3);
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {F(i, j) = vecF[3*i + j];}
        }
        if (det3(F) < 0) {F *= -1;}                      // resolve the sign ambiguity of the null vector
        F = rbsd::optimal_rotation(F);                   // snap the seed onto SO(3)

        // refine to a valid frame by conjugation averaging
        for (int iteration = 0; iteration < 200; ++iteration) {
            Matrix<double> B(3, 3);
            for (std::size_t k = 1; k < R.size(); ++k) {B += R[k]*F*G[k].transpose();}
            Matrix<double> F_next = rbsd::optimal_rotation(B);
            double delta = 0;
            for (unsigned int i = 0; i < 3; ++i) {
                for (unsigned int j = 0; j < 3; ++j) {delta += std::abs(F_next(i, j) - F(i, j));}
            }
            F = std::move(F_next);
            if (delta < 1e-12) {break;}
        }
        return F;
    }

    double reconstruction_rmsd(
        const symmetry::ISymmetry& sym, const Vector3<double>& cm, const std::vector<std::vector<Vector3<double>>>& copies
    ) {
        double sum_sq = 0;
        std::size_t count = 0;
        for (unsigned int k = 1; k < copies.size(); ++k) {
            auto transform = sym.get_transform(cm, k);
            for (std::size_t i = 0; i < copies[0].size(); ++i) {
                Vector3<double> d = transform(copies[0][i]) - copies[k][i];
                sum_sq += d.dot(d);
                ++count;
            }
        }
        return count == 0 ? 0.0 : std::sqrt(sum_sq/static_cast<double>(count));
    }
}

namespace ausaxs::rigidbody::sequencer::detail {
SymmetryFitResult fit_symmetry(
    const symmetry::ISymmetry& template_symmetry, const Vector3<double>& cm, const std::vector<std::vector<Vector3<double>>>& copies
) {
    assert(copies.size() == template_symmetry.repetitions() + 1 && "fit_symmetry: body count must equal repetitions()+1.");
    assert(!copies.empty() && !copies[0].empty() && "fit_symmetry: empty input.");

    // per-copy alignment: exact Kabsch superposition of the reference onto each copy
    std::vector<Matrix<double>> R;      // R[k] maps copies[0] onto copies[k]; R[0] is the identity
    std::vector<Vector3<double>> t;
    R.reserve(copies.size());
    t.reserve(copies.size());
    R.emplace_back(Matrix<double>::identity(3));
    t.emplace_back(Vector3<double>{0, 0, 0});
    for (std::size_t k = 1; k < copies.size(); ++k) {
        assert(copies[k].size() == copies[0].size() && "fit_symmetry: all copies must have the same number of atoms.");
        auto s = superpose(copies[0], copies[k]);
        R.push_back(std::move(s.rotation));
        t.push_back(s.translation);
    }

    // shared rotation centre (skip R[0], the identity, which carries no information)
    Vector3<double> c = solve_common_centre(
        {R.begin() + 1, R.end()}, {t.begin() + 1, t.end()}, cm
    );

    SymmetryFitResult result;
    result.symmetry = template_symmetry.clone();
    symmetry::ISymmetry* sym = result.symmetry.get();

    if (auto* ps = dynamic_cast<symmetry::PointSymmetry*>(sym)) {
        // single copy: v -> R(v - cm) + cm + d, so d = t_1 - (I - R_1) cm
        ps->rotation = matrix::euler_angles(R[1]);
        ps->translation = t[1] - (Matrix<double>::identity(3) - R[1])*cm;
    }
    else if (auto* cs = dynamic_cast<symmetry::CyclicSymmetry*>(sym)) {
        // generator R_1 defines the axis/angle; the body offset is cm - c (point group, no screw)
        auto [axis, angle] = axis_angle(R[1]);
        cs->_repeat_relation.axis = axis;
        cs->_repeat_relation.angle = angle;
        cs->_repeat_relation.translation = {0, 0, 0};
        cs->_initial_relation.translation = cm - c;
    }
    else if (auto* poly = dynamic_cast<symmetry::IPolyhedralSymmetry*>(sym)) {
        Matrix<double> F = solve_frame(R, poly->group_elements());
        poly->rotation = matrix::euler_angles(F);
        if (poly->span_translation().size() == 2) {
            // planar dihedral: the offset lives in the group frame with the axial (z) component pinned to 0
            Vector3<double> v = F.transpose()*(c - cm);
            poly->translation = {v.x(), v.y(), 0};
        } else {
            poly->translation = c - cm;
        }
    }
    else {
        assert(false && "fit_symmetry: unsupported target symmetry type.");
    }

    result.rmsd = reconstruction_rmsd(*sym, cm, copies);
    return result;
}
}
