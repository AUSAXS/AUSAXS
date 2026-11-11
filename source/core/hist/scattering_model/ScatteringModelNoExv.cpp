// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/scattering_model/ScatteringModelNoExv.h>

#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <utility/MultiThreading.h>

#include <numeric>

using namespace ausaxs;
using namespace ausaxs::hist;

const form_factor::lookup::table_t& ScatteringModelNoExv::get_ff_table() const {
    return form_factor::manager::get_active_product_tables()->raw_atomic_table;
}

std::vector<double> ScatteringModelNoExv::debye_transform_direct(const std::vector<double>& q) const {
    auto sinqd_table = sinc_table_for(q);
    const auto& ff_indices = form_factor::manager::get_active_product_tables()->ff_indices;
    const int n_ff = form_factor::get_active_count();
    constexpr int first = form_factor::start_index_for_explicit_exv();
    constexpr int water = form_factor::water_bin;
    const double cw = free_params.cw;

    std::vector<double> I(q.size(), 0);
    auto* pool = utility::multi_threading::get_global_pool();
    pool->detach_blocks(0, static_cast<int>(q.size()), [&] (int start, int end) {
        std::vector<double> f(n_ff), s(n_ff);
        for (int k = start; k < end; ++k) {
            const double* sinqd = sinqd_table->begin(k);
            auto debye_sum = [sinqd] (auto begin, auto end) {return std::transform_reduce(begin, end, sinqd, 0.0);};
            for (int ff = first; ff < n_ff; ++ff) {
                auto type = static_cast<form_factor::form_factor_t>(ff_indices[ff]);
                f[ff] = form_factor::manager::evaluate_amplitude(type, q[k]);
                s[ff] = form_factor::manager::evaluate_self(type, q[k]);
            }

            // as in cache_refresh_intensity_profiles: the first bin of a diagonal partial also holds the self-correlations, which the 
            // sums evaluate with the cross-correlation form factor at the distance of that bin, so they are swapped for the self-correlation one
            double aa = 0, aw = 0, ww = 0;
            for (int ff1 = first; ff1 < n_ff; ++ff1) {
                for (int ff2 = ff1; ff2 < n_ff; ++ff2) {
                    aa += f[ff1]*f[ff2]*debye_sum(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2));
                }
                aa += self_weights.aa[ff1]*(s[ff1] - sinqd[0]*f[ff1]*f[ff1]);
                aw += 2*cw*f[ff1]*f[water]*debye_sum(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1));
            }
            ww += cw*cw*f[water]*f[water]*debye_sum(distance_profiles.ww.begin(), distance_profiles.ww.end());
            ww += cw*cw*self_weights.ww*(s[water] - sinqd[0]*f[water]*f[water]);

            double B = get_atomic_debye_waller_factor(q[k], free_params.DW_sigma_atomic);
            I[k] = aa*B*B + aw*B + ww;
        }
    });
    pool->wait();
    return I;
}
