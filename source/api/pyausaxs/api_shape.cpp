// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_shape.h>

#include <form_factor/ExvFormFactor.h>
#include <grid/detail/GridExcludedVolume.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <math/Vector3.h>
#include <settings/InternalState.h>
#include <utility/Exceptions.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace ausaxs;

void shape_debye_userq(
    const double* x, const double* y, const double* z, int n_cells, double spacing,
    const double* q, double* I, int n_q,
    int* status
) {execute_with_catch([&]() {
#if !defined(POCKETFFT_AVAILABLE)
    throw except::disabled("shape_debye_userq: this build has no lattice transform. Configure with -DPOCKETFFT=ON.");
#else
    if (n_cells <= 0) {throw except::invalid_argument("shape_debye_userq: the shape has no cells.");}
    if (!(0 < spacing)) {throw except::invalid_argument("shape_debye_userq: the lattice spacing must be positive.");}

    // snap the centres to integer lattice sites; anything off the lattice is a caller error, not something to round away
    const double* coords[3] = {x, y, z};
    std::vector<Vector3<int>> sites(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        for (int k = 0; k < 3; ++k) {
            double s = coords[k][i]/spacing;
            double site = std::round(s);
            if (1e-3 < std::abs(s - site)) {
                throw except::invalid_argument(
                    "shape_debye_userq: cell " + std::to_string(i) + " is not on a lattice of spacing " + std::to_string(spacing) + "."
                );
            }
            sites[i][k] = static_cast<int>(site);
        }
    }

    // the transform needs non-negative sites, and would count a doubly occupied site as a pair at distance zero
    Vector3<int> min = sites[0], extent{0, 0, 0};
    for (const auto& s : sites) {for (int k = 0; k < 3; ++k) {min[k] = std::min(min[k], s[k]);}}
    for (auto& s : sites) {for (int k = 0; k < 3; ++k) {s[k] -= min[k]; extent[k] = std::max(extent[k], s[k]);}}
    {
        auto sorted = sites;
        auto order = [] (const Vector3<int>& a, const Vector3<int>& b) {return std::ranges::lexicographical_compare(a, b);};
        std::ranges::sort(sorted, order);
        if (std::ranges::adjacent_find(sorted) != sorted.end()) {throw except::invalid_argument("shape_debye_userq: two cells occupy the same lattice site.");}
    }

    grid::exv::GridExcludedVolume exv;
    exv.interior.resize(n_cells); // only the sites and the spacing enter the transform
    exv.interior_sites = std::move(sites);
    exv.spacing = spacing;

    double inv_width = settings::internal_state::inv_bin_width;
    double max_distance = spacing*std::sqrt(static_cast<double>(extent.x())*extent.x() + static_cast<double>(extent.y())*extent.y() + static_cast<double>(extent.z())*extent.z());
    int bin_count = static_cast<int>(std::ceil(max_distance*inv_width)) + 2;
    auto p = hist::detail::lattice::self_correlation(exv, inv_width, bin_count);
    p.add_index(0, hist::detail::WeightedEntry(n_cells, n_cells, 0)); // self-pairs
    auto counts = p.get_content();
    auto d = p.get_weighted_axis();

    double V = spacing*spacing*spacing;
    form_factor::ExvFormFactor ff(V);
    for (int j = 0; j < n_q; ++j) {
        double sum = 0;
        for (int k = 0; k < static_cast<int>(counts.size()); ++k) {
            if (counts[k] == 0) {continue;}
            double qd = q[j]*d[k];
            sum += counts[k]*(qd < 1e-8 ? 1 : std::sin(qd)/qd);
        }
        double f = V*ff.evaluate_normalized(q[j]);
        I[j] = sum*f*f;
    }
#endif
}, status);}
