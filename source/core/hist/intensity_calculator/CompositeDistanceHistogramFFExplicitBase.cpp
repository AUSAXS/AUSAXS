/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicitBase.h>
#include <hist/Histogram.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <form_factor/PrecalculatedExvFormFactorProduct.h>
#include <settings/HistogramSettings.h>
#include <constants/ConstantsMath.h>
#include <math/ConstexprMath.h>
#include <dataset/SimpleDataset.h>
#include <utility/MultiThreading.h>

using namespace hist;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase() = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(const CompositeDistanceHistogramFFExplicitBase&) = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::CompositeDistanceHistogramFFExplicitBase(CompositeDistanceHistogramFFExplicitBase&&) noexcept = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>& CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::operator=(CompositeDistanceHistogramFFExplicitBase&&) noexcept = default;

template<typename AA, typename AX, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>& CompositeDistanceHistogramFFExplicitBase<AA, AX, XX>::operator=(const CompositeDistanceHistogramFFExplicitBase&) = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::~CompositeDistanceHistogramFFExplicitBase() = default;

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AA>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), exv_distance_profiles{.xx=std::move(p_xx), .ax=std::move(p_ax), .wx=std::move(p_wx)} {}

template<typename AA, typename AXFormFactorTableType, typename XX>
CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::CompositeDistanceHistogramFFExplicitBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution3D&& p_ax, 
    hist::Distribution3D&& p_xx, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution2D&& p_wx, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : CompositeDistanceHistogramFFAvgBase<AA>(std::move(p_aa), std::move(p_aw), std::move(p_ww), std::move(p_tot)), exv_distance_profiles{.xx=std::move(p_xx), .ax=std::move(p_ax), .wx=std::move(p_wx)} {}

template<typename AA, typename AXFormFactorTableType, typename XX>
const AA CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::get_ffaa_table() const {
    return this->get_ff_table();
}

template<typename AA, typename AXFormFactorTableType, typename XX>
double CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::exv_factor(double q) const {
    constexpr double rm2 = constants::radius::average_atomic_radius*constants::radius::average_atomic_radius/4;
    return std::pow(this->free_params.cx, 3)*std::exp(-rm2*(std::pow(this->free_params.cx, 2) - 1)*q*q);
}

template<typename AA, typename AXFormFactorTableType, typename XX>
void CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::cache_refresh_sinqd() const {
    auto pool = utility::multi_threading::get_global_pool();
    auto sinqd_table = this->get_sinc_table();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (exv_cache.sinqd.aa.empty()) {
        exv_cache.sinqd.aa = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        exv_cache.sinqd.ax = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        exv_cache.sinqd.xx = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        exv_cache.sinqd.aw = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        exv_cache.sinqd.wx = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        exv_cache.sinqd.ww = container::Container1D<double>(debye_axis.bins);
    }

    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            pool->detach_task([this, q0, bins=debye_axis.bins, ff1, ff2, sinqd_table] () {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    exv_cache.sinqd.aa.index(ff1, ff2, q-q0) = std::inner_product(this->distance_profiles.aa.begin(ff1, ff2), this->distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    exv_cache.sinqd.ax.index(ff1, ff2, q-q0) = std::inner_product(exv_distance_profiles.ax.begin(ff1, ff2), exv_distance_profiles.ax.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                    exv_cache.sinqd.xx.index(ff1, ff2, q-q0) = std::inner_product(exv_distance_profiles.xx.begin(ff1, ff2), exv_distance_profiles.xx.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                }
            });
        }
        pool->detach_task([this, q0, bins=debye_axis.bins, ff1, sinqd_table] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                exv_cache.sinqd.aw.index(ff1, q-q0) = std::inner_product(this->distance_profiles.aw.begin(ff1), this->distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
                exv_cache.sinqd.wx.index(ff1, q-q0) = std::inner_product(exv_distance_profiles.wx.begin(ff1), exv_distance_profiles.wx.end(ff1), sinqd_table->begin(q), 0.0);
            }
        });
    }
    pool->detach_task([&] () {
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            exv_cache.sinqd.ww.index(q-q0) = std::inner_product(this->distance_profiles.ww.begin(), this->distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
        }
    });
    exv_cache.sinqd.valid = true;
    pool->wait();
}

template<typename AA, typename AXFormFactorTableType, typename XX>
void CompositeDistanceHistogramFFExplicitBase<AA, AXFormFactorTableType, XX>::cache_refresh_intensity_profiles(bool sinqd_changed, bool cw_changed, bool cx_changed) const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& ff_aa_table = get_ffaa_table();
    const auto& ff_ax_table = get_ffax_table();
    const auto& ff_xx_table = get_ffxx_table();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin); // account for a possibly different qmin

    if (sinqd_changed) {
        this->cache.intensity_profiles.aa = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed) {
        this->cache.intensity_profiles.aw = std::vector<double>(debye_axis.bins, 0);
        this->cache.intensity_profiles.ww = std::vector<double>(debye_axis.bins, 0);
    }
    if (cx_changed) {
        this->cache.intensity_profiles.ax = std::vector<double>(debye_axis.bins, 0);
        this->cache.intensity_profiles.xx = std::vector<double>(debye_axis.bins, 0);
    }
    if (cw_changed || cx_changed) {
        this->cache.intensity_profiles.wx = std::vector<double>(debye_axis.bins, 0);
    }

    std::vector<double> cx(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {cx[q-q0] = exv_factor(constants::axes::q_vals[q]);}

    if (sinqd_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                        this->cache.intensity_profiles.aa[q-q0] += 
                            exv_cache.sinqd.aa.index(ff1, ff2, q-q0)*ff_aa_table.index(ff1, ff2).evaluate(q);
                    }
                }
            }
        });
    }

    if (cx_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                        this->cache.intensity_profiles.ax[q-q0] += 
                            2*cx[q]*exv_cache.sinqd.ax.index(ff1, ff2, q-q0)*ff_ax_table.index(ff1, ff2).evaluate(q);
                    }
                }
            }
        });
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                        this->cache.intensity_profiles.xx[q-q0] += 
                            cx[q]*cx[q]*exv_cache.sinqd.xx.index(ff1, ff2, q-q0)*ff_xx_table.index(ff1, ff2).evaluate(q);
                    }
                }
            }
        });
    }

    if (cw_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    this->cache.intensity_profiles.aw[q-q0] += 
                        2*this->free_params.cw*exv_cache.sinqd.aw.index(ff1, q-q0)
                        *ff_aa_table.index(ff1, form_factor::water_bin).evaluate(q);
                }
            }
        });
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                this->cache.intensity_profiles.ww[q-q0] += 
                    this->free_params.cw*this->free_params.cw*exv_cache.sinqd.ww.index(q-q0)
                    *ff_aa_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }

    if (cw_changed || cx_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    this->cache.intensity_profiles.wx[q-q0] += 
                        2*cx[q]*this->free_params.cw*exv_cache.sinqd.wx.index(ff1, q-q0)
                        *ff_ax_table.index(form_factor::water_bin, ff1).evaluate(q);
                }
            }
        });
    }
    this->cache.intensity_profiles.cached_cx = this->free_params.cx;
    this->cache.intensity_profiles.cached_cw = this->free_params.cw;
    pool->wait();    
}

template class hist::CompositeDistanceHistogramFFExplicitBase<
    form_factor::storage::atomic::table_t, form_factor::storage::cross::table_t, form_factor::storage::exv::table_t
>;