// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/distribution/Distribution2D.h>
#include <hist/distribution/Distribution3D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>

#include <algorithm>
#include <vector>

using namespace ausaxs;
using namespace ausaxs::hist;

namespace {
    void apply_averaged_exv(Distribution3D& p_aa, Distribution2D& p_aw, double Z_exv_avg) {
        const int n_active = static_cast<int>(form_factor::get_active_count());
        const int ff_start = form_factor::start_index_for_explicit_exv();
        const int nbins = static_cast<int>(p_aa.size_z());

        std::vector<double> xx(nbins, 0);
        for (int a = ff_start; a < n_active; ++a) {
            std::vector<double> ax(nbins, 0);
            for (int b = ff_start; b < n_active; ++b) {
                for (int d = 0; d < nbins; ++d) {
                    ax[d] += p_aa.index(a, b, d);
                    ax[d] += p_aa.index(b, a, d);
                    xx[d] += p_aa.index(a, b, d);
                }
            }
            for (int d = 0; d < nbins; ++d) {
                p_aa.index(a, form_factor::exv_bin, d) = 0.5*ax[d];
            }
        }
        for (int d = 0; d < nbins; ++d) {
            p_aa.index(form_factor::exv_bin, form_factor::exv_bin, d) = xx[d];
        }

        for (int a = ff_start; a < n_active; ++a) {
            for (int d = 0; d < nbins; ++d) {
                p_aw.index(form_factor::exv_bin, d) += p_aw.index(a, d);
            }
        }

        // multiply the excluded volume charge onto the excluded volume bins
        for (int ff1 = ff_start; ff1 < n_active; ++ff1) {
            std::transform(p_aa.begin(ff1, form_factor::exv_bin), p_aa.end(ff1, form_factor::exv_bin), p_aa.begin(ff1, form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg;});
        }
        std::transform(p_aa.begin(form_factor::exv_bin, form_factor::exv_bin), p_aa.end(form_factor::exv_bin, form_factor::exv_bin), p_aa.begin(form_factor::exv_bin, form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg*Z_exv_avg;});
        std::transform(p_aw.begin(form_factor::exv_bin), p_aw.end(form_factor::exv_bin), p_aw.begin(form_factor::exv_bin), [Z_exv_avg] (auto val) {return val*Z_exv_avg;});
    }
}

std::unique_ptr<CompositeDistanceHistogramFFAvg> CompositeDistanceHistogramFFAvg::with_averaged_exv(
    hist::Distribution3D&& p_aa,
    hist::Distribution2D&& p_aw,
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot,
    double Z_exv_avg
) {
    apply_averaged_exv(p_aa, p_aw, Z_exv_avg);
    return std::make_unique<CompositeDistanceHistogramFFAvg>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot));
}

std::unique_ptr<CompositeDistanceHistogramFFAvg> CompositeDistanceHistogramFFAvg::with_averaged_exv(
    hist::Distribution3D&& p_aa,
    hist::Distribution2D&& p_aw,
    hist::Distribution1D&& p_ww,
    hist::WeightedDistribution1D&& p_tot,
    double Z_exv_avg
) {
    apply_averaged_exv(p_aa, p_aw, Z_exv_avg);
    return std::make_unique<CompositeDistanceHistogramFFAvg>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot));
}
