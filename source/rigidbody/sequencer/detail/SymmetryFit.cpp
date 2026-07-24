// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SymmetryFit.h>
#include <math/Superposition.h>
#include <math/SymmetricEigen.h>
#include <math/Matrix.h>
#include <data/symmetry/PointSymmetry.h>
#include <data/symmetry/CyclicSymmetry.h>
#include <data/symmetry/IPolyhedralSymmetry.h>
#include <data/symmetry/CompositeSymmetry.h>
#include <math/MatrixUtils.h>

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace ausaxs;

namespace {
    // Find the common rotation centre c by solving
    //     min_c Σ_k ||(I - R_k)c - t_k||² + λ||c - cm||².
    // The small regularisation towards cm fixes the otherwise-undetermined component along the symmetry axis.
    Vector3<double> solve_common_centre(const std::vector<Matrix<double>>& R, const std::vector<Vector3<double>>& t, const Vector3<double>& cm) {
        // Solve by forming the normal equations:
        //   (A^T A + lambda I)c = A^T t + lambda cm
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
        return matrix::inverse(normal)*rhs;
    }

    // Recover the group frame F relating the ideal symmetry operators G[k] to the observed copy rotations R[k] = F G[k] F^T.
    // An initial estimate is obtained from the linear constraints R[k]F - FG[k] = 0. This seed is projected onto SO(3) and then refined by iteratively 
    // averaging the conjugation relation above. The result is a rotation matrix F that best maps the ideal group frame onto the observed one.
    Matrix<double> solve_frame(const std::vector<Matrix<double>>& R, const std::vector<Matrix<double>>& G) {
        assert(R.size() == G.size() && "solve_frame: rotation/element count mismatch.");
        Matrix<double> Q(9, 9);
        for (std::size_t k = 1; k < R.size(); ++k) { // element 0 is the identity: L is zero
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

        auto eig = matrix::symmetric_eigen(Q);
        const auto& vecF = eig.vectors.front(); // smallest eigenvalue -> initial null-vector estimate
        Matrix<double> F(3, 3);
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {F(i, j) = vecF[3*i + j];}
        }
        if (matrix::det(F) < 0) {F *= -1;} // resolve the sign ambiguity of the null vector
        F = matrix::optimal_rotation(F); // snap the seed onto SO(3)

        // refine to a valid frame by conjugation averaging
        for (int iteration = 0; iteration < 200; ++iteration) {
            Matrix<double> B(3, 3);
            for (std::size_t k = 1; k < R.size(); ++k) {B += R[k]*F*G[k].transpose();}
            Matrix<double> F_next = matrix::optimal_rotation(B);
            double delta = 0;
            for (unsigned int i = 0; i < 3; ++i) {
                for (unsigned int j = 0; j < 3; ++j) {delta += std::abs(F_next(i, j) - F(i, j));}
            }
            F = std::move(F_next);
            if (delta < 1e-12) {break;}
        }
        return F;
    }

    // Detect the trivial cyclic solution (identity repeat transform). 
    // Rather than treating it as a valid symmetry, reconstruct using the identity map and retain the resulting large fit residual.    
    bool is_trivial_cyclic(const symmetry::ISymmetry& sym) {
        auto* cs = dynamic_cast<const symmetry::CyclicSymmetry*>(&sym);
        return cs
            && std::abs(cs->_repeat_relation.angle) < 1e-6
            && cs->_repeat_relation.translation.magnitude() < 1e-6;
    }

    double reconstruction_rmsd(
        const symmetry::ISymmetry& sym, const Vector3<double>& cm, const std::vector<std::vector<Vector3<double>>>& copies
    ) {
        auto reconstructed = rigidbody::sequencer::detail::reconstruct_copies(sym, cm, copies[0]);
        double sum_sq = 0;
        std::size_t count = 0;
        for (std::size_t k = 1; k < copies.size(); ++k) {
            for (std::size_t i = 0; i < copies[0].size(); ++i) {
                Vector3<double> d = reconstructed[k][i] - copies[k][i];
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

    // A composite factorises: copy (k, j) = outer_k(inner_j(reference)), and since inner_0 and outer_0 are the identity, copy (0, j) is the inner symmetry's 
    // own copy j and copy (k, 0) the outer symmetry's own copy k. So each sub-symmetry can be fitted independently from the appropriate slice of the assembly, 
    // recursing to arbitrary nesting depth.
    if (auto* comp = dynamic_cast<const symmetry::CompositeSymmetry*>(&template_symmetry)) {
        int stride = 1 + static_cast<int>(comp->inner->repetitions());
        int outer_reps = static_cast<int>(comp->outer->repetitions());

        std::vector<std::vector<Vector3<double>>> inner_copies(copies.begin(), copies.begin() + stride);
        auto inner_fit = fit_symmetry(*comp->inner, cm, inner_copies);

        std::vector<std::vector<Vector3<double>>> outer_copies;
        outer_copies.reserve(outer_reps + 1);
        for (int k = 0; k <= outer_reps; ++k) {outer_copies.push_back(copies[k*stride]);}
        auto outer_fit = fit_symmetry(*comp->outer, cm, outer_copies);

        SymmetryFitResult result;
        result.symmetry = std::make_unique<symmetry::CompositeSymmetry>(std::move(inner_fit.symmetry), std::move(outer_fit.symmetry));
        result.rmsd = reconstruction_rmsd(*result.symmetry, cm, copies);
        return result;
    }

    // per-copy alignment: exact Kabsch superposition of the reference onto each copy
    std::vector<Matrix<double>> R; // R[k] maps copies[0] onto copies[k]; R[0] is the identity
    std::vector<Vector3<double>> t;
    R.reserve(copies.size());
    t.reserve(copies.size());
    R.emplace_back(Matrix<double>::identity(3));
    t.emplace_back(Vector3<double>{0, 0, 0});
    for (std::size_t k = 1; k < copies.size(); ++k) {
        assert(copies[k].size() == copies[0].size() && "fit_symmetry: all copies must have the same number of atoms.");
        auto s = matrix::superpose(copies[0], copies[k]);
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
        auto [axis, angle] = matrix::axis_angle(R[1]);
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

SymmetryFitResult fit_symmetry_best_order(
    const symmetry::ISymmetry& template_symmetry, const Vector3<double>& cm,
    const std::vector<std::vector<Vector3<double>>>& copies, double accept_rmsd
) {
    SymmetryFitResult best = fit_symmetry(template_symmetry, cm, copies);
    if (best.rmsd <= accept_rmsd) {return best;}

    // Enumerating orderings of the non-reference bodies costs (n-1)! fits; cap it so pathological sizes fall back to the given order rather than hanging. 
    // Highly-constrained large-order groups (which is where n grows) are precisely the cases that do not need this feature.
    int n = static_cast<int>(copies.size());
    if (9 < n) {return best;}

    std::vector<int> perm(n - 1);
    for (int i = 0; i < n - 1; ++i) {perm[i] = i + 1;} // reference (index 0) stays fixed

    std::vector<std::vector<Vector3<double>>> reordered(n);
    reordered[0] = copies[0];
    while (std::next_permutation(perm.begin(), perm.end())) { // the identity order was already tried
        for (int i = 0; i < n - 1; ++i) {reordered[i + 1] = copies[perm[i]];}
        auto candidate = fit_symmetry(template_symmetry, cm, reordered);
        if (candidate.rmsd < best.rmsd) {best = std::move(candidate);}
        if (best.rmsd <= accept_rmsd) {break;}
    }
    return best;
}

std::vector<std::vector<Vector3<double>>> reconstruct_copies(
    const symmetry::ISymmetry& symmetry, const Vector3<double>& cm, const std::vector<Vector3<double>>& reference
) {
    bool trivial = is_trivial_cyclic(symmetry);
    unsigned int reps = symmetry.repetitions();
    std::vector<std::vector<Vector3<double>>> result;
    result.reserve(reps + 1);
    result.push_back(reference); // copy 0 is the reference itself
    for (unsigned int k = 1; k <= reps; ++k) {
        std::function<Vector3<double>(Vector3<double>)> transform =
            trivial ? [](Vector3<double> v) {return v;} : symmetry.get_transform(cm, k);
        std::vector<Vector3<double>> copy;
        copy.reserve(reference.size());
        for (const auto& p : reference) {copy.push_back(transform(p));}
        result.push_back(std::move(copy));
    }
    return result;
}
}
