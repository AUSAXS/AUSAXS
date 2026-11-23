// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringAmplitude.h>
#include <constants/ConstantsAxes.h>
#include <math/ConstexprMath.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <cassert>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    /**
     * @brief The weights of a group of scatterers, binned by their distance from a centre. Each bin keeps the mean distance of its scatterers.
     *        Only the occupied bins are kept.
     */
    struct BinnedDistances {
        std::vector<double> weights, distances;
    };

    BinnedDistances bin_distances(const hist::detail::CompactCoordinates& coords, const Vector3<double>& centre) {
        constexpr double width = constants::axes::d_axis.width();

        // the bins are weighted by the absolute weights, so a bin of mixed signs still has its mean distance inside it
        std::vector<double> w, abs_w, abs_wr;
        for (int i = 0; i < coords.size(); ++i) {
            double r = Vector3<double>(coords.x(i), coords.y(i), coords.z(i)).distance(centre);
            auto bin = static_cast<std::size_t>(r/width);
            if (w.size() <= bin) {
                w.resize(bin+1, 0);
                abs_w.resize(bin+1, 0);
                abs_wr.resize(bin+1, 0);
            }
            double weight = coords.get_weight(i);
            w[bin] += weight;
            abs_w[bin] += std::abs(weight);
            abs_wr[bin] += std::abs(weight)*r;
        }

        BinnedDistances res;
        for (std::size_t bin = 0; bin < w.size(); ++bin) {
            if (abs_w[bin] == 0) {continue;}
            res.weights.emplace_back(w[bin]);
            res.distances.emplace_back(abs_wr[bin]/abs_w[bin]);
        }
        return res;
    }
}

ScatteringAmplitude::ScatteringAmplitude(std::vector<hist::detail::CompactCoordinates>&& groups) : groups(std::move(groups)) {}

std::vector<double> ScatteringAmplitude::evaluate(
    const std::vector<double>& q, const std::optional<Vector3<double>>& centre, const ScatteringParameters& params
) const {
    auto origin = centre.has_value() ? *centre : default_centre(params);
    std::vector<BinnedDistances> bins;
    bins.reserve(groups.size());
    for (const auto& group : groups) {bins.emplace_back(bin_distances(group, origin));}

    // the form factors are evaluated up front, since they are not required to be thread-safe
    std::vector<std::vector<double>> f(q.size());
    for (std::size_t k = 0; k < q.size(); ++k) {
        f[k] = form_factors(q[k], params);
        assert(f[k].size() == groups.size() && "ScatteringAmplitude::evaluate: form_factors must return one value per group.");
    }

    std::vector<double> A(q.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), [&q, &f, &bins, &A] (int start, int end) {
        for (int k = start; k < end; ++k) {
            for (std::size_t g = 0; g < bins.size(); ++g) {
                if (f[k][g] == 0) {continue;}
                const auto& b = bins[g];
                double sum = 0;
                for (std::size_t i = 0; i < b.distances.size(); ++i) {sum += b.weights[i]*math::sinc(q[k]*b.distances[i]);}
                A[k] += f[k][g]*sum;
            }
        }
    });
    pool->wait();
    return A;
}

Vector3<double> ScatteringAmplitude::default_centre(const ScatteringParameters& params) const {
    auto f0 = form_factors(0, params);
    assert(f0.size() == groups.size() && "ScatteringAmplitude::default_centre: form_factors must return one value per group.");
    Vector3<double> sum{0, 0, 0}, abs_sum{0, 0, 0};
    double total = 0, abs_total = 0;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        const auto& coords = groups[g];
        for (int i = 0; i < coords.size(); ++i) {
            double w = f0[g]*coords.get_weight(i);
            Vector3<double> r(coords.x(i), coords.y(i), coords.z(i));
            sum += r*w;
            abs_sum += r*std::abs(w);
            total += w;
            abs_total += std::abs(w);
        }
    }
    if (abs_total == 0) {throw except::invalid_argument("ScatteringAmplitude::default_centre: The molecule has no scattering length.");}
    if (std::abs(total) < 1e-6*abs_total) {return abs_sum/abs_total;}
    return sum/total;
}
