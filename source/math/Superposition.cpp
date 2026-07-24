// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <math/Superposition.h>
#include <math/SymmetricEigen.h>

#include <cassert>
#include <cmath>

using namespace ausaxs;

namespace ausaxs::matrix {
// Orthogonal Procrustes solution (the rotation maximizing the Frobenius product <R, H>), computed via
// Horn's quaternion method: the optimal rotation is the largest eigenvector of a symmetric 4x4 matrix
// built from H. Going through a quaternion guarantees a proper rotation, never a reflection.
Matrix<double> optimal_rotation(const Matrix<double>& H) {
    assert(H.N == 3 && H.M == 3 && "optimal_rotation: H must be 3x3.");
    double Sxx = H(0,0), Sxy = H(0,1), Sxz = H(0,2);
    double Syx = H(1,0), Syy = H(1,1), Syz = H(1,2);
    double Szx = H(2,0), Szy = H(2,1), Szz = H(2,2);

    // Horn's symmetric 4x4 key matrix; its largest eigenvector is the optimal rotation quaternion
    Matrix<double> K = {
        {Sxx + Syy + Szz, Syz - Szy,        Szx - Sxz,        Sxy - Syx       },
        {Syz - Szy,       Sxx - Syy - Szz,  Sxy + Syx,        Szx + Sxz       },
        {Szx - Sxz,       Sxy + Syx,       -Sxx + Syy - Szz,  Syz + Szy       },
        {Sxy - Syx,       Szx + Sxz,        Syz + Szy,       -Sxx - Syy + Szz }
    };

    // largest eigenvalue is last (symmetric_eigen returns ascending order)
    auto eig = symmetric_eigen(K);
    const auto& q = eig.vectors.back();
    double w = q[0], x = q[1], y = q[2], z = q[3];

    // normalize defensively before forming the rotation matrix
    double norm = std::sqrt(w*w + x*x + y*y + z*z);
    if (norm < 1e-12) {return Matrix<double>::identity(3);}
    w /= norm; x /= norm; y /= norm; z /= norm;

    return Matrix<double>{
        {1 - 2*(y*y + z*z), 2*(x*y + w*z),     2*(x*z - w*y)    },
        {2*(x*y - w*z),     1 - 2*(x*x + z*z), 2*(y*z + w*x)    },
        {2*(x*z + w*y),     2*(y*z - w*x),     1 - 2*(x*x + y*y)}
    };
}

// Kabsch algorithm: centre both point sets, then recover the optimal rotation from their
// cross-covariance and the translation that lines the centroids up.
SuperpositionResult superpose(const std::vector<Vector3<double>>& from, const std::vector<Vector3<double>>& to) {
    assert(from.size() == to.size() && "superpose: point sets must have equal size.");
    assert(!from.empty() && "superpose: point sets must be non-empty.");
    std::size_t n = from.size();

    // centroids
    Vector3<double> from_c{0, 0, 0}, to_c{0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) {from_c += from[i]; to_c += to[i];}
    from_c /= static_cast<double>(n);
    to_c   /= static_cast<double>(n);

    // cross-covariance H = sum (to - to_c)(from - from_c)^T, so R maps from -> to
    Matrix<double> H(3, 3);
    for (std::size_t i = 0; i < n; ++i) {
        Vector3<double> f = from[i] - from_c;
        Vector3<double> t = to[i]   - to_c;
        for (unsigned int r = 0; r < 3; ++r) {
            for (unsigned int c = 0; c < 3; ++c) {H(r, c) += t[r]*f[c];}
        }
    }

    SuperpositionResult result;
    result.rotation = optimal_rotation(H);
    result.translation = to_c - result.rotation*from_c;

    // residual RMSD
    double sum_sq = 0;
    for (std::size_t i = 0; i < n; ++i) {
        Vector3<double> d = result.rotation*from[i] + result.translation - to[i];
        sum_sq += d.dot(d);
    }
    result.rmsd = std::sqrt(sum_sq/static_cast<double>(n));
    return result;
}
}
