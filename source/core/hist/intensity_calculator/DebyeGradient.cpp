// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/DebyeGradient.h>

#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <algorithm>
#include <cmath>

using namespace ausaxs;

namespace {
    // H(r) only contains frequencies up to the largest q, so linear interpolation is accurate to about (q_max*table_width)^2/8
    constexpr double table_width = 0.01;
    constexpr double inv_table_width = 1/table_width;

    // each block seeds its own trigonometric recurrence, which bounds the accumulated rounding
    constexpr int table_block = 512;

    // below this qr, cos(x) - sinc(x) cancels too badly, and its Taylor series is used instead
    constexpr double series_limit = 0.1;

    /**
     * @brief Tabulate H(r) = 2 \sum_q v(q) (cos(qr) - sinc(qr))/r^2 at r = k*table_width for k in [0, size).
     */
    std::vector<double> tabulate_pair_derivative(const std::vector<double>& q, const std::vector<double>& v, int size) {
        int n_q = static_cast<int>(q.size());
        std::vector<double> q2(n_q), inv_q(n_q), v2(n_q), cd(n_q), sd(n_q);
        for (int i = 0; i < n_q; ++i) {
            q2[i] = q[i]*q[i];
            inv_q[i] = q[i] == 0 ? 0 : 1/q[i];
            v2[i] = 2*v[i];
            cd[i] = std::cos(q[i]*table_width);
            sd[i] = std::sin(q[i]*table_width);
        }

        std::vector<double> H(size);
        auto* pool = utility::multi_threading::get_global_pool();
        int n_blocks = (size + table_block - 1)/table_block;
        pool->detach_blocks(0, n_blocks, [&] (int block_start, int block_end) {
            std::vector<double> c(n_q), s(n_q);
            for (int block = block_start; block < block_end; ++block) {
                int k0 = block*table_block;
                int k1 = std::min(size, k0 + table_block);
                for (int i = 0; i < n_q; ++i) {
                    c[i] = std::cos(q[i]*k0*table_width);
                    s[i] = std::sin(q[i]*k0*table_width);
                }

                for (int k = k0; k < k1; ++k) {
                    double r = k*table_width;
                    double inv_r = k == 0 ? 0 : 1/r;
                    double inv_r2 = inv_r*inv_r;
                    double sum = 0;
                    for (int i = 0; i < n_q; ++i) {
                        double x = q[i]*r;
                        double x2 = x*x;
                        double series = q2[i]*(-1.0/3 + x2*(1.0/30 - x2/840));
                        double direct = (c[i] - s[i]*inv_q[i]*inv_r)*inv_r2;
                        sum += v2[i]*(x < series_limit ? series : direct);
                    }
                    H[k] = sum;

                    // advance to the next r by rotating (cos, sin) by q*table_width
                    for (int i = 0; i < n_q; ++i) {
                        double cn = c[i]*cd[i] - s[i]*sd[i];
                        s[i] = s[i]*cd[i] + c[i]*sd[i];
                        c[i] = cn;
                    }
                }
            }
        });
        pool->wait();
        return H;
    }

    /**
     * @brief Accumulate the contributions of all pairs (i, j > i) of row i into the gradient buffers.
     */
    void accumulate_row(
        int i, int n, const double* x, const double* y, const double* z, const double* w, const double* H,
        double* gx, double* gy, double* gz
    ) {
        const double xi = x[i], yi = y[i], zi = z[i], wi = w[i];
        double ax = 0, ay = 0, az = 0;
        for (int j = i+1; j < n; ++j) {
            double dx = xi - x[j];
            double dy = yi - y[j];
            double dz = zi - z[j];
            double t = std::sqrt(dx*dx + dy*dy + dz*dz)*inv_table_width;
            int k = static_cast<int>(t);
            double h = H[k] + (t - k)*(H[k+1] - H[k]);
            double f = wi*w[j]*h;
            ax += f*dx; ay += f*dy; az += f*dz;
            gx[j] -= f*dx; gy[j] -= f*dy; gz[j] -= f*dz;
        }
        gx[i] += ax; gy[i] += ay; gz[i] += az;
    }
}

std::vector<Vector3<double>> hist::debye_raw_vjp(const data::Molecule& molecule, const std::vector<double>& q, const std::vector<double>& v) {
    if (q.size() != v.size()) {throw except::size_error("debye_raw_vjp: q and v must have the same length.");}

    // the same representation the binned calculation reads, so the weights are guaranteed to agree with it
    auto data = hist::detail::factory::construct_from_atoms<false>(&molecule);
    const int n = data.size();
    std::vector<Vector3<double>> gradient(n, {0, 0, 0});
    if (n < 2 || q.empty()) {return gradient;}

    std::vector<double> x(n), y(n), z(n), w(n);
    Vector3<double> min{data.x(0), data.y(0), data.z(0)}, max = min;
    for (int i = 0; i < n; ++i) {
        x[i] = data.x(i); y[i] = data.y(i); z[i] = data.z(i); w[i] = data.get_weight(i);
        min = {std::min(min.x(), x[i]), std::min(min.y(), y[i]), std::min(min.z(), z[i])};
        max = {std::max(max.x(), x[i]), std::max(max.y(), y[i]), std::max(max.z(), z[i])};
    }

    // no pair is farther apart than the bounding box diagonal; the margin covers the interpolation neighbour and rounding
    int table_size = static_cast<int>(max.distance(min)*inv_table_width) + 3;
    auto H = tabulate_pair_derivative(q, v, table_size);

    // pairing row m with row n-1-m gives every task the same n-1 pairs
    auto* pool = utility::multi_threading::get_global_pool();
    container::ThreadLocalWrapper<std::vector<double>> buffers(3*n, 0.0);
    int folded = (n+1)/2;
    pool->detach_blocks(0, folded, [&] (int start, int end) {
        auto& buffer = buffers.get();
        double* gx = buffer.data();
        double* gy = gx + n;
        double* gz = gy + n;
        for (int m = start; m < end; ++m) {
            accumulate_row(m, n, x.data(), y.data(), z.data(), w.data(), H.data(), gx, gy, gz);
            if (int mirror = n-1-m; mirror != m) {
                accumulate_row(mirror, n, x.data(), y.data(), z.data(), w.data(), H.data(), gx, gy, gz);
            }
        }
    }, 4*pool->get_thread_count());
    pool->wait();

    for (const std::vector<double>& buffer : buffers.get_all()) {
        for (int i = 0; i < n; ++i) {
            gradient[i] += Vector3<double>{buffer[i], buffer[n+i], buffer[2*n+i]};
        }
    }
    return gradient;
}
