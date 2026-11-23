// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_shape.h>

#include <dataset/Dataset.h>
#include <form_factor/ExvFormFactor.h>
#include <grid/detail/GridExcludedVolume.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <hist/scattering_model/ScatteringModelPoint.h>
#include <math/Vector3.h>
#include <settings/InternalState.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace ausaxs;

void shape_debye_userq(
    const double* x, const double* y, const double* z, int n_cells, double spacing,
    const double* q, double* I, int n_q,
    int* status
) {
    shape_debye_rectangular_userq(x, y, z, n_cells, spacing, spacing, spacing, q, I, n_q, status);
}

void shape_debye_rectangular_userq(
    const double* x, const double* y, const double* z, int n_cells,
    double spacing_x, double spacing_y, double spacing_z,
    const double* q, double* I, int n_q,
    int* status
) {execute_with_catch([&]() {
#if !defined(POCKETFFT_AVAILABLE)
    throw except::disabled("shape_debye_rectangular_userq: this build has no lattice transform. Configure with -DPOCKETFFT=ON.");
#else
    if (n_cells <= 0) {throw except::invalid_argument("shape_debye_rectangular_userq: the shape has no cells.");}
    if (n_q <= 0) {throw except::invalid_argument("shape_debye_rectangular_userq: no q values were supplied.");}
    std::array<double, 3> spacing = {spacing_x, spacing_y, spacing_z};
    for (double step : spacing) {
        if (!std::isfinite(step) || step <= 0) {
            throw except::invalid_argument("shape_debye_rectangular_userq: each lattice spacing must be finite and positive.");
        }
    }
    double V = spacing[0]*spacing[1]*spacing[2];
    if (!std::isfinite(V) || V == 0) {
        throw except::invalid_argument("shape_debye_rectangular_userq: cell volume is not representable.");
    }

    // snap the centres to integer lattice sites relative to the first cell, so the lattice may have any origin; anything off it is a caller error
    std::array<const double*, 3> coords = {x, y, z};
    std::vector<Vector3<int>> sites(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        for (int k = 0; k < 3; ++k) {
            double s = (coords[k][i] - coords[k][0])/spacing[k];
            // Leave room for shifting negative sites and padding the resulting extent.
            constexpr double max_site = std::numeric_limits<int>::max()/4 - 1;
            if (!std::isfinite(s) || std::abs(s) > max_site) {
                throw except::invalid_argument("shape_debye_rectangular_userq: lattice coordinates are not finite or are too large.");
            }
            double site = std::round(s);
            if (1e-3 < std::abs(s - site)) {
                throw except::invalid_argument(
                    "shape_debye_rectangular_userq: cell " + std::to_string(i) + " is not on the lattice of cell 0 along axis " + std::to_string(k) + "."
                );
            }
            sites[i][k] = static_cast<int>(site);
        }
    }

    // the transform needs non-negative sites
    Vector3<int> min = sites[0], extent{0, 0, 0};
    for (const auto& s : sites) {for (int k = 0; k < 3; ++k) {min[k] = std::min(min[k], s[k]);}}
    for (auto& s : sites) {for (int k = 0; k < 3; ++k) {s[k] -= min[k]; extent[k] = std::max(extent[k], s[k]);}}
    {   // sort the sites and check for duplicates, which would be a caller error
        auto sorted = sites;
        auto order = [] (const Vector3<int>& a, const Vector3<int>& b) {return std::ranges::lexicographical_compare(a, b);};
        std::ranges::sort(sorted, order);
        if (std::ranges::adjacent_find(sorted) != sorted.end()) {throw except::invalid_argument("shape_debye_rectangular_userq: two cells occupy the same lattice site.");}
    }

    double inv_width = settings::internal_state::inv_bin_width;
    double max_distance = std::hypot(extent.x()*spacing[0], extent.y()*spacing[1], extent.z()*spacing[2]);
    if (!std::isfinite(max_distance*inv_width) || max_distance*inv_width > std::numeric_limits<int>::max() - 2.) {
        throw except::invalid_argument("shape_debye_rectangular_userq: distance histogram is too large.");
    }
    int bin_count = static_cast<int>(std::ceil(max_distance*inv_width)) + 2;
    auto p = hist::detail::lattice::self_correlation(sites, spacing, inv_width, bin_count);
    p.add_index(0, hist::detail::WeightedEntry(n_cells, n_cells, 0)); // self-pairs

    // every cell scatters with the Gaussian form factor of its volume, as in the Grid excluded volume model
    std::vector<double> q_vals(q, q + n_q);
    auto debye_I = hist::ScatteringModelPoint(std::move(p)).debye_transform(q_vals);
    form_factor::detail::ExvFormFactor ff(V, 1); // unit contrast: only the shape is modelled
    for (int j = 0; j < n_q; ++j) {
        double f = ff.evaluate(q[j]);
        I[j] = debye_I.y(j)*f*f;
    }
#endif
}, status);}
