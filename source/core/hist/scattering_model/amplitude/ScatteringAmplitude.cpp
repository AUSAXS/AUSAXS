// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/amplitude/ScatteringAmplitude.h>

#include <constants/ConstantsAxes.h>
#include <math/ConstexprMath.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    /**
     * @brief The terms w_i sinc(q r_i) of the amplitude sum of one group, with its scatterers merged into thin shells about the centre.
     */
    struct SincTerms {
        std::vector<double> weights, distances;
    };

    SincTerms make_sinc_terms(const hist::detail::CompactCoordinates& coords, const Vector3<double>& centre) {
        constexpr double width = constants::axes::d_axis.width();

        // the bins are weighted by the absolute weights, so a bin of mixed signs still has its mean distance inside it
        std::vector<double> w, abs_w, abs_wr;
        for (int i = 0; i < coords.size(); ++i) {
            double r = Vector3<double>{coords.x(i), coords.y(i), coords.z(i)}.distance(centre);
            auto bin = static_cast<int>(r/width);
            if (static_cast<int>(w.size()) <= bin) { // grow the vectors to accommodate the new bin
                w.resize(bin+1, 0);
                abs_w.resize(bin+1, 0);
                abs_wr.resize(bin+1, 0);
            }
            double weight = coords.get_weight(i);
            w[bin] += weight;
            abs_w[bin] += std::abs(weight);
            abs_wr[bin] += std::abs(weight)*r;
        }

        SincTerms res;
        for (int bin = 0; bin < static_cast<int>(w.size()); ++bin) {
            if (abs_w[bin] == 0) {continue;} // an empty bin adds nothing, and has no mean distance
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
    std::vector<SincTerms> terms;
    terms.reserve(groups.size());
    std::ranges::transform(groups, std::back_inserter(terms), [&origin] (const auto& group) {return make_sinc_terms(group, origin);});

    // the form factors are evaluated up front, since they are not required to be thread-safe
    const auto f = get_excess_form_factors(q, params);
    assert(
        f.size_x() == static_cast<int>(groups.size()) && f.size_y() == static_cast<int>(q.size()) 
        && "ScatteringAmplitude::evaluate: get_excess_form_factors must return one row per group and one column per q value."
    );

    std::vector<double> A(q.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), [&q, &f, &terms, &A] (int start, int end) {
        for (int k = start; k < end; ++k) {
            for (int g = 0; g < static_cast<int>(terms.size()); ++g) {
                const auto& t = terms[g];
                double sum = std::transform_reduce(
                    t.weights.begin(), t.weights.end(), t.distances.begin(), 0.0, std::plus<>(),
                    [qk = q[k]] (double w, double r) {return w*math::sinc(qk*r);}
                );
                A[k] += f(g, k)*sum;
            }
        }
    });
    pool->wait();
    return A;
}

Vector3<double> ScatteringAmplitude::default_centre(const ScatteringParameters& params) const {
    const auto f0 = get_excess_form_factors(std::vector<double>{0}, params);
    assert(f0.size_x() == static_cast<int>(groups.size()) && "ScatteringAmplitude::default_centre: get_excess_form_factors must return one row per group.");
    Vector3<double> sum{0, 0, 0};
    double total = 0;
    for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
        const auto& coords = groups[g];
        for (int i = 0; i < coords.size(); ++i) {
            double w = f0(g, 0)*coords.get_weight(i);
            sum += Vector3<double>(coords.x(i), coords.y(i), coords.z(i))*w;
            total += w;
        }
    }
    assert(total != 0 && "ScatteringAmplitude::default_centre: The molecule has no contrast.");
    return sum/total;
}
