#include <hist/intensity_calculator/waxsis/CompositeDistanceHistogramWAXSiS.h>
#include <utility/MultiThreading.h>
#include <settings/HistogramSettings.h>

using namespace ausaxs;
using namespace ausaxs::hist;

CompositeDistanceHistogramWAXSiS::~CompositeDistanceHistogramWAXSiS() = default;

void CompositeDistanceHistogramWAXSiS::cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& ff_table = get_ff_table(); 

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (sinqd_changed) {
        cache.intensity_profiles.aa = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed) {
        cache.intensity_profiles.aw = std::vector<double>(debye_axis.bins, 0);
        cache.intensity_profiles.ww = std::vector<double>(debye_axis.bins, 0);
    }
    if (cx_changed) {
        cache.intensity_profiles.ax = std::vector<double>(debye_axis.bins, 0);
        cache.intensity_profiles.xx = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed || cx_changed) {
        cache.intensity_profiles.wx = std::vector<double>(debye_axis.bins, 0);
    }

    // calculate exv factor
    std::vector<double> cx(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
        cx[q-q0] = CompositeDistanceHistogramFFGrid::exv_factor(constants::axes::q_vals[q], free_params.cx);
    }

    if (sinqd_changed) {
        // aa
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                        cache.intensity_profiles.aa[q-q0] += cache.sinqd.aa.index(ff1, ff2, q-q0)*ff_table.index(ff1, ff2).evaluate(q);
                    }
                }
            }
        });
    }

    if (cx_changed) {
        // ax
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    cache.intensity_profiles.ax[q-q0] += 2*free_params.crho*cx[q-q0]*cache.sinqd.ax.index(ff1, q-q0)*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
                }
            }
        });

        // xx
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.xx[q-q0] += std::pow(cx[q-q0]*free_params.crho, 2)*cache.sinqd.xx.index(q-q0)*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }

    if (cw_changed) {
        // aw
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    cache.intensity_profiles.aw[q-q0] += 2*free_params.cw*cache.sinqd.aw.index(ff1, q-q0)*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
                }
            }
        });

        // ww
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.ww[q-q0] += free_params.cw*free_params.cw*cache.sinqd.ww.index(q-q0)*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }

    if (cw_changed || cx_changed) {
        // wx
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.wx[q-q0] += 2*free_params.crho*cx[q-q0]*free_params.cw*cache.sinqd.wx.index(q-q0)*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }

    cache.intensity_profiles.cached_cx = free_params.cx;
    cache.intensity_profiles.cached_crho = free_params.crho;
    cache.intensity_profiles.cached_cw = free_params.cw;
    pool->wait();
}