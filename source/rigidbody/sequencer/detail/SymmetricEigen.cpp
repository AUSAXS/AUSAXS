// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/SymmetricEigen.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

using namespace ausaxs;

namespace ausaxs::rigidbody::sequencer::detail {
EigenResult symmetric_eigen(const Matrix<double>& A) {
    assert(A.N == A.M && "symmetric_eigen: matrix must be square.");
    const unsigned int n = A.N;

    // working copy of the (symmetric) matrix; it is driven towards diagonal form
    std::vector<double> a(A.data.begin(), A.data.end());
    auto at = [&](unsigned int i, unsigned int j) -> double& {return a[i*n + j];};

    // accumulated eigenvectors (columns of V)
    std::vector<double> v(n*n, 0);
    for (unsigned int i = 0; i < n; ++i) {v[i*n + i] = 1;}
    auto vt = [&](unsigned int i, unsigned int j) -> double& {return v[i*n + j];};

    // cyclic Jacobi sweeps
    constexpr int max_sweeps = 100;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        // convergence check: sum of squared off-diagonal entries
        double off = 0;
        for (unsigned int p = 0; p < n; ++p) {
            for (unsigned int q = p+1; q < n; ++q) {off += at(p, q)*at(p, q);}
        }
        if (off < 1e-30) {break;}

        for (unsigned int p = 0; p < n; ++p) {
            for (unsigned int q = p+1; q < n; ++q) {
                double apq = at(p, q);
                if (std::abs(apq) < 1e-300) {continue;}

                // Jacobi rotation angle zeroing the (p,q) entry
                double app = at(p, p), aqq = at(q, q);
                double phi = 0.5*(aqq - app)/apq;
                double t = (phi >= 0 ? 1.0 : -1.0)/(std::abs(phi) + std::sqrt(phi*phi + 1));
                double c = 1/std::sqrt(t*t + 1);
                double s = t*c;

                // apply rotation to columns then rows p and q of the working matrix
                for (unsigned int k = 0; k < n; ++k) {
                    double akp = at(k, p), akq = at(k, q);
                    at(k, p) = c*akp - s*akq;
                    at(k, q) = s*akp + c*akq;
                }
                for (unsigned int k = 0; k < n; ++k) {
                    double apk = at(p, k), aqk = at(q, k);
                    at(p, k) = c*apk - s*aqk;
                    at(q, k) = s*apk + c*aqk;
                }

                // accumulate the rotation into the eigenvector matrix
                for (unsigned int k = 0; k < n; ++k) {
                    double vkp = vt(k, p), vkq = vt(k, q);
                    vt(k, p) = c*vkp - s*vkq;
                    vt(k, q) = s*vkp + c*vkq;
                }
            }
        }
    }

    // gather (eigenvalue, eigenvector) pairs and sort ascending by eigenvalue
    std::vector<unsigned int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](unsigned int i, unsigned int j) {return at(i, i) < at(j, j);});

    EigenResult result;
    result.values.reserve(n);
    result.vectors.reserve(n);
    for (unsigned int idx : order) {
        result.values.push_back(at(idx, idx));
        std::vector<double> vec(n);
        for (unsigned int k = 0; k < n; ++k) {vec[k] = vt(k, idx);}
        result.vectors.push_back(std::move(vec));
    }
    return result;
}
}
