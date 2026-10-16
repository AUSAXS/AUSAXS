// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/DebyeGradient.h>

#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <settings/GeneralSettings.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <algorithm>
#include <cmath>

using namespace ausaxs;

namespace {
    constexpr double table_width = 0.01;
    constexpr double inv_table_width = 1/table_width;
    constexpr double series_limit = 0.1;

    /**
     * @brief Tabulate H(r) = 2 \sum_q v(q) (cos(qr) - sinc(qr))/r^2 at r = k*table_width for k in [0, size).
     */
    std::vector<float> tabulate_pair_derivative(const std::vector<double>& q, const std::vector<double>& v, int size) {
        std::vector<float> H(size);
        auto* pool = utility::multi_threading::get_global_pool();
        pool->detach_loop(0, size, [&] (int k) {
            double r = k*table_width;
            double sum = 0;
            for (unsigned int i = 0; i < q.size(); ++i) {
                double x = q[i]*r;
                double x2 = x*x;
                sum += v[i]*(x < series_limit
                    ? q[i]*q[i]*(-1.0/3 + x2*(1.0/30 - x2/840))
                    : (std::cos(x) - std::sin(x)/x)/(r*r)
                );
            }
            H[k] = static_cast<float>(2*sum);
        });
        pool->wait();
        return H;
    }
}

std::vector<Vector3<double>> hist::debye_raw_vjp(const data::Molecule& molecule, const std::vector<double>& q, const std::vector<double>& v) {
    if (q.size() != v.size()) {throw except::size_error("debye_raw_vjp: q and v must have the same length.");}

    // the same representation the binned calculation reads, so the weights are guaranteed to agree with it
    auto data = hist::detail::factory::construct_from_atoms<false>(&molecule);
    const int n = data.size();
    std::vector<Vector3<double>> gradient(n, {0, 0, 0});
    if (n < 2 || q.empty()) {return gradient;}

    // the margin covers the interpolation neighbour and rounding
    auto H = tabulate_pair_derivative(q, v, static_cast<int>(hist::detail::bin_estimate::max_distance(data)*inv_table_width) + 3);

    auto* pool = utility::multi_threading::get_global_pool();
    container::ThreadLocalWrapper<std::vector<float>> buffers(3*n, 0.f);
    int job_size = settings::general::detail::get_job_size(n);
    for (int i = 0; i < n; i += job_size) {
        pool->detach_task([&, imin = i, imax = std::min(i+job_size, n)] () {
            // the table is copied into an allocation local to the task, since the compiler can only prove that the gathers from it never alias the gradient 
            // buffer stores when it sees the allocation; otherwise the loop below is not vectorized through a raw pointer, gcc extracts each index once and 
            // reaches k+1 by an address offset; through table[k] it extracts both k and k+1 per lane, which is measurably slower
            std::vector<float> table = H; // NOLINT
            const float* h = table.data();
            auto& g = buffers.get();
            for (int i = imin; i < imax; ++i) {
                float xi = data.x(i), yi = data.y(i), zi = data.z(i), wi = data.get_weight(i);
                float sx = 0, sy = 0, sz = 0;
                for (int j = i+1; j < n; ++j) {
                    float dx = xi - data.x(j), dy = yi - data.y(j), dz = zi - data.z(j);
                    float t = std::sqrt(dx*dx + dy*dy + dz*dz)*static_cast<float>(inv_table_width);
                    int k = static_cast<int>(t);
                    float s = wi*data.get_weight(j)*(h[k] + (t - static_cast<float>(k))*(h[k+1] - h[k]));
                    float fx = dx*s, fy = dy*s, fz = dz*s;
                    sx += fx; sy += fy; sz += fz;
                    g[j] -= fx; g[n+j] -= fy; g[2*n+j] -= fz;
                }
                g[i] += sx; g[n+i] += sy; g[2*n+i] += sz;
            }
        });
    }
    pool->wait();

    for (const std::vector<float>& g : buffers.get_all()) {
        for (int i = 0; i < n; ++i) {gradient[i] += Vector3<double>{g[i], g[i+n], g[i+2*n]};}
    }
    return gradient;
}
