// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <mini/LevenbergMarquardt.h>

#include <mini/detail/FittedParameter.h>
#include <mini/detail/Parameter.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::mini;

namespace {
    constexpr double ftol = 1e-10;      // stop when an accepted step lowers chi2 by less than this fraction
    constexpr double gtol = 1e-10;      // stop when the scaled gradient is this small relative to chi2
    constexpr double xtol = 1e-10;      // stop when a step moves every parameter by less than this fraction of its scale
    constexpr double lambda_max = 1e16; // damping this strong means no step can lower chi2 any further
    const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());

    /**
     * @brief Solve A x = b for a small symmetric positive definite matrix by Cholesky decomposition.
     *
     * @param A The row-major n x n matrix. Taken by value since it is overwritten by its factor.
     * @param b The right-hand side, overwritten by the solution.
     * @return false if A is not positive definite, in which case b is left in an unspecified state.
     */
    bool cholesky_solve(std::vector<double> A, std::vector<double>& b, int n) {
        for (int j = 0; j < n; ++j) {
            double d = A[j*n+j];
            for (int k = 0; k < j; ++k) {d -= A[j*n+k]*A[j*n+k];}
            if (!(0 < d)) {return false;}
            d = std::sqrt(d);
            A[j*n+j] = d;
            for (int i = j+1; i < n; ++i) {
                double s = A[i*n+j];
                for (int k = 0; k < j; ++k) {s -= A[i*n+k]*A[j*n+k];}
                A[i*n+j] = s/d;
            }
        }

        // L y = b, then L^T x = y
        for (int i = 0; i < n; ++i) {
            double s = b[i];
            for (int k = 0; k < i; ++k) {s -= A[i*n+k]*b[k];}
            b[i] = s/A[i*n+i];
        }
        for (int i = n-1; 0 <= i; --i) {
            double s = b[i];
            for (int k = i+1; k < n; ++k) {s -= A[k*n+i]*b[k];}
            b[i] = s/A[i*n+i];
        }
        return true;
    }
}

LevenbergMarquardt::LevenbergMarquardt(residual_function func, const std::vector<Parameter>& params) : Minimizer(std::move(func)) {
    for (const auto& p : params) {add_parameter(p);}
}

Result LevenbergMarquardt::minimize_override() {
    const int n = static_cast<int>(parameters.size());

    // starting point, bounds, and the scale below which a parameter's magnitude is not used for its step sizes
    constexpr double inf = std::numeric_limits<double>::infinity();
    std::vector<double> x(n), lo(n, -inf), hi(n, inf), typ(n, 1);
    for (int j = 0; j < n; ++j) {
        const auto& p = parameters[j];
        if (p.bounds.has_value()) {
            lo[j] = std::min(p.bounds->min, p.bounds->max);
            hi[j] = std::max(p.bounds->min, p.bounds->max);
            typ[j] = 1e-2*(hi[j] - lo[j]);
        }
        x[j] = std::clamp(p.guess.has_value() ? *p.guess : p.bounds->center(), lo[j], hi[j]);
    }

    auto evaluate = [this] (const std::vector<double>& p) {
        auto r = residuals(p);
        double f = chi2(r);
        return std::make_pair(std::move(r), f);
    };

    auto [r, F] = evaluate(x);
    const int m = static_cast<int>(r.size());
    std::vector<double> J(m*n);         // column-major, J[j*m + i] = dr_i/dx_j
    std::vector<double> A(n*n, 0), g(n, 0);
    double lambda = 1e-3, nu = 2;
    bool converged = false;
    while (!converged && fevals < max_evals) {
        // forward-difference Jacobian, stepping inwards from an upper bound
        for (int j = 0; j < n; ++j) {
            double h = sqrt_eps*std::max(std::abs(x[j]), typ[j]);
            if (hi[j] < x[j] + h) {h = -h;}
            auto xh = x;
            xh[j] += h;
            auto [rh, _] = evaluate(xh);
            for (int i = 0; i < m; ++i) {J[j*m+i] = (rh[i] - r[i])/h;}
        }

        // normal equations: A = J^T J, g = J^T r, so the gradient of chi2 is 2g
        for (int a = 0; a < n; ++a) {
            g[a] = std::inner_product(J.begin()+a*m, J.begin()+(a+1)*m, r.begin(), 0.0);
            for (int b = 0; b <= a; ++b) {
                A[a*n+b] = A[b*n+a] = std::inner_product(J.begin()+a*m, J.begin()+(a+1)*m, J.begin()+b*m, 0.0);
            }
        }

        // freeze parameters sitting on a bound with the descent direction -g pointing out of the box
        std::vector<int> free;
        double gmax = 0, dmax = 0;
        for (int j = 0; j < n; ++j) {
            if ((x[j] <= lo[j] && 0 < g[j]) || (hi[j] <= x[j] && g[j] < 0)) {continue;}
            free.push_back(j);
            gmax = std::max(gmax, std::abs(g[j])*std::max(std::abs(x[j]), typ[j]));
            dmax = std::max(dmax, A[j*n+j]);
        }
        if (free.empty() || gmax <= gtol*F) {converged = true; break;}

        // find a step which lowers chi2, adjusting the damping as we go
        const int k = static_cast<int>(free.size());
        const double dmin = 1e-12*dmax;
        bool accepted = false;
        while (!accepted && !converged && fevals < max_evals) {
            // Marquardt's damping: (A + lambda diag(A)) d = -g
            std::vector<double> M(k*k), d(k);
            for (int a = 0; a < k; ++a) {
                for (int b = 0; b < k; ++b) {M[a*k+b] = A[free[a]*n+free[b]];}
                M[a*k+a] += lambda*std::max(A[free[a]*n+free[a]], dmin);
                d[a] = -g[free[a]];
            }
            if (!cholesky_solve(std::move(M), d, k)) {
                lambda *= nu; nu *= 2;
                converged = lambda_max < lambda;
                continue;
            }

            // project the step onto the box
            auto xn = x;
            std::vector<double> s(n, 0);
            bool small = true;
            for (int a = 0; a < k; ++a) {
                int j = free[a];
                xn[j] = std::clamp(x[j] + d[a], lo[j], hi[j]);
                s[j] = xn[j] - x[j];
                small &= std::abs(s[j]) <= xtol*(std::abs(x[j]) + typ[j]);
            }
            if (small) {converged = true; break;}

            // the reduction predicted by the linearized model, F - |r + Js|^2
            double pred = 0;
            for (int a = 0; a < n; ++a) {
                pred -= 2*g[a]*s[a];
                for (int b = 0; b < n; ++b) {pred -= s[a]*A[a*n+b]*s[b];}
            }

            auto [rn, Fn] = evaluate(xn);
            double rho = 0 < pred ? (F - Fn)/pred : -1;
            if (0 < rho) {
                converged = F - Fn <= ftol*F;
                x = std::move(xn);
                r = std::move(rn);
                F = Fn;
                lambda *= std::max(1./3, 1 - std::pow(2*rho - 1, 3)); // Nielsen's update
                nu = 2;
                accepted = true;
            } else {
                lambda *= nu; nu *= 2;
                converged = lambda_max < lambda;
            }
        }
    }

    // parameter errors from the curvature, cov = (J^T J)^-1, using the most recent Jacobian
    Result res;
    for (int j = 0; j < n; ++j) {
        std::vector<double> e(n, 0);
        e[j] = 1;
        double err = cholesky_solve(A, e, n) && 0 < e[j] ? std::sqrt(e[j]) : 0;
        res.add_parameter(FittedParameter(parameters[j], x[j], err));
    }
    res.fval = F;
    res.fevals = fevals;
    res.status = converged ? 0 : 1;
    return res;
}
