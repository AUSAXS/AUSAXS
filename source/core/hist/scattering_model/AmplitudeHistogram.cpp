// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/AmplitudeHistogram.h>
#include <constants/ConstantsAxes.h>
#include <math/ConstexprMath.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <cassert>
#include <cmath>

using namespace ausaxs;
using namespace ausaxs::hist;

AmplitudeHistogram::AmplitudeHistogram(const AmplitudeScatterers& scatterers, const Vector3<double>& centre) {
    constexpr double width = constants::axes::d_axis.width();
    classes.resize(scatterers.classes.size());
    for (std::size_t c = 0; c < scatterers.classes.size(); ++c) {
        const auto& coords = scatterers.classes[c];

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

        for (std::size_t bin = 0; bin < w.size(); ++bin) {
            if (abs_w[bin] == 0) {continue;}
            classes[c].weights.emplace_back(w[bin]);
            classes[c].distances.emplace_back(abs_wr[bin]/abs_w[bin]);
        }
    }
}

Vector3<double> AmplitudeHistogram::scattering_centre(const AmplitudeScatterers& scatterers, const std::vector<double>& g0) {
    assert(g0.size() == scatterers.classes.size() && "AmplitudeHistogram::scattering_centre: g0 must hold one value per class.");
    Vector3<double> sum{0, 0, 0}, abs_sum{0, 0, 0};
    double total = 0, abs_total = 0;
    for (std::size_t c = 0; c < scatterers.classes.size(); ++c) {
        const auto& coords = scatterers.classes[c];
        for (int i = 0; i < coords.size(); ++i) {
            double w = g0[c]*coords.get_weight(i);
            Vector3<double> r(coords.x(i), coords.y(i), coords.z(i));
            sum += r*w;
            abs_sum += r*std::abs(w);
            total += w;
            abs_total += std::abs(w);
        }
    }
    if (abs_total == 0) {throw except::invalid_argument("AmplitudeHistogram::scattering_centre: The molecule has no scattering length.");}
    if (std::abs(total) < 1e-6*abs_total) {return abs_sum/abs_total;}
    return sum/total;
}

std::vector<double> AmplitudeHistogram::transform(const std::vector<double>& q, const std::function<std::vector<double>(double)>& form_factors) const {
    // the form factors are evaluated up front, since a model is not required to evaluate them concurrently
    std::vector<std::vector<double>> g(q.size());
    for (std::size_t k = 0; k < q.size(); ++k) {
        g[k] = form_factors(q[k]);
        assert(g[k].size() == classes.size() && "AmplitudeHistogram::transform: form_factors must return one value per class.");
    }

    std::vector<double> A(q.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), [this, &q, &g, &A] (int start, int end) {
        for (int k = start; k < end; ++k) {
            for (std::size_t c = 0; c < classes.size(); ++c) {
                if (g[k][c] == 0) {continue;}
                const auto& cls = classes[c];
                double sum = 0;
                for (std::size_t i = 0; i < cls.distances.size(); ++i) {sum += cls.weights[i]*math::sinc(q[k]*cls.distances[i]);}
                A[k] += g[k][c]*sum;
            }
        }
    });
    pool->wait();
    return A;
}
