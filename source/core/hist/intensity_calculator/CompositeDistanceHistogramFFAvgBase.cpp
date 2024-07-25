/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvgBase.h>
#include <hist/Histogram.h>
#include <dataset/SimpleDataset.h>
#include <table/ArrayDebyeTable.h>
#include <form_factor/FormFactor.h>
#include <form_factor/PrecalculatedFormFactorProduct.h>
#include <utility/MultiThreading.h>
#include <settings/HistogramSettings.h>

#include <cassert>

using namespace hist;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase() = default;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(CompositeDistanceHistogramFFAvgBase&&) noexcept = default;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::operator=(CompositeDistanceHistogramFFAvgBase&&) noexcept = default;

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww,
    hist::Distribution1D&& p_tot
) : ICompositeDistanceHistogramExv(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::CompositeDistanceHistogramFFAvgBase(
    hist::Distribution3D&& p_aa, 
    hist::Distribution2D&& p_aw, 
    hist::Distribution1D&& p_ww, 
    hist::WeightedDistribution1D&& p_tot
) : ICompositeDistanceHistogramExv(std::move(p_tot)), distance_profiles{.aa=std::move(p_aa), .aw=std::move(p_aw), .ww=std::move(p_ww)} {}

template<typename FormFactorTableType>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::~CompositeDistanceHistogramFFAvgBase() = default;

template<typename FormFactorTableType>
double CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::exv_factor(double) const {
    return free_params.cx;
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform() const {
    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    std::vector<double> Iq(debye_axis.bins, 0);
    auto[aa, ax, aw, xx, wx, ww] = cache_get_intensity_profiles();
    assert(aa.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): aa.size() != Iq.size()");
    assert(ax.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): ax.size() != Iq.size()");
    assert(aw.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): aw.size() != Iq.size()");
    assert(xx.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): xx.size() != Iq.size()");
    assert(wx.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): wx.size() != Iq.size()");
    assert(ww.size() == Iq.size() && "CompositeDistanceHistogramFFAvgBase::debye_transform(): ww.size() != Iq.size()");
    std::transform(Iq.begin(), Iq.end(), aa.begin(), Iq.begin(), std::plus<>());
    std::transform(Iq.begin(), Iq.end(), ax.begin(), Iq.begin(), std::minus<>());
    std::transform(Iq.begin(), Iq.end(), aw.begin(), Iq.begin(), std::plus<>());
    std::transform(Iq.begin(), Iq.end(), xx.begin(), Iq.begin(), std::plus<>());
    std::transform(Iq.begin(), Iq.end(), wx.begin(), Iq.begin(), std::minus<>());
    std::transform(Iq.begin(), Iq.end(), ww.begin(), Iq.begin(), std::plus<>());
    return ScatteringProfile(Iq, debye_axis);
}

template<typename FormFactorTableType>
SimpleDataset CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::debye_transform(const std::vector<double>&) const {
    throw except::not_implemented("CompositeDistanceHistogramFFGrid::debye_transform(const std::vector<double>& q) const");
}

template<typename FormFactorTableType>
const std::vector<double>& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_counts() const {
    p = std::vector<double>(DistanceHistogram::get_counts().size(), 0);
    auto& p_pp = get_aa_counts();
    auto& p_hp = get_aw_counts();
    auto& p_hh = get_ww_counts();
    for (unsigned int i = 0; i < p.size(); ++i) {
        p[i] = p_pp.index(i) + 2*free_params.cw*p_hp.index(i) + free_params.cw*free_params.cw*p_hh.index(i);
    }
    return p.data;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() {
    return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aa_counts());
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts() const {
    auto[aa, _, __] = cache_get_distance_profiles();
    return aa;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_aw_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts() const {
    auto[_, aw, __] = cache_get_distance_profiles();
    return aw;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() {return const_cast<Distribution1D&>(const_cast<const CompositeDistanceHistogramFFAvgBase*>(this)->get_ww_counts());}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts() const {
    auto[_, __, ww] = cache_get_distance_profiles();
    return ww;
}

template<typename FormFactorTableType>
const Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() const {
    return distance_profiles.aa;
}

template<typename FormFactorTableType>
Distribution3D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aa_counts_ff() {
    return distance_profiles.aa;
}

template<typename FormFactorTableType>
const Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() const {
    return distance_profiles.aw;
}

template<typename FormFactorTableType>
Distribution2D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_aw_counts_ff() {
    return distance_profiles.aw;
}

template<typename FormFactorTableType>
const Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() const {
    return distance_profiles.ww;
}

template<typename FormFactorTableType>
Distribution1D& CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_ww_counts_ff() {
    return distance_profiles.ww;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_water_scaling_factor(double k) {
    free_params.cw = k;
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::apply_excluded_volume_scaling_factor(double k) {
    free_params.cx = k;
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aa() const {
    std::vector<double> aa;
    std::tie(aa, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore) = cache_get_intensity_profiles();
    return ScatteringProfile(aa, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ax() const {
    std::vector<double> ax;
    std::tie(std::ignore, ax, std::ignore, std::ignore, std::ignore, std::ignore) = cache_get_intensity_profiles();
    return ScatteringProfile(ax, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_xx() const {
    std::vector<double> xx;
    std::tie(std::ignore, std::ignore, std::ignore, xx, std::ignore, std::ignore) = cache_get_intensity_profiles();
    return ScatteringProfile(xx, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_wx() const {
    std::vector<double> wx;
    std::tie(std::ignore, std::ignore, std::ignore, std::ignore, wx, std::ignore) = cache_get_intensity_profiles();
    return ScatteringProfile(wx, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_aw() const {
    std::vector<double> aw;
    std::tie(std::ignore, std::ignore, aw, std::ignore, std::ignore, std::ignore) = cache_get_intensity_profiles();
    return ScatteringProfile(aw, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
ScatteringProfile CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::get_profile_ww() const {
    std::vector<double> ww;
    std::tie(std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, ww) = cache_get_intensity_profiles();
    return ScatteringProfile(ww, constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax));
}

template<typename FormFactorTableType>
std::tuple<const Distribution1D&, const Distribution1D&,const Distribution1D&> 
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::cache_get_distance_profiles() const {
    if (!cache.distance_profiles.valid) {cache_refresh_distance_profiles();}
    return std::forward_as_tuple(cache.distance_profiles.p_aa, cache.distance_profiles.p_aw, cache.distance_profiles.p_ww);
}

template<typename FormFactorTableType>
std::tuple<const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&>
CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::cache_get_intensity_profiles() const {
    if (!cache.sinqd.valid) {
        cache_refresh_sinqd();
        cache_refresh_intensity_profiles<true, true, true>();
        cache.sinqd.valid = true;
    } else {
        if (cache.intensity_profiles.cached_cx != free_params.cx) {
            cache.intensity_profiles.cached_cw != free_params.cw ? cache_refresh_intensity_profiles<false, true, true>() : cache_refresh_intensity_profiles<false, false, true>();
        } else if (cache.intensity_profiles.cached_cw != free_params.cw) {
            cache_refresh_intensity_profiles<false, true, false>();
        }
    }

    return std::forward_as_tuple(
        cache.intensity_profiles.aa, cache.intensity_profiles.ax, cache.intensity_profiles.aw, 
        cache.intensity_profiles.xx, cache.intensity_profiles.wx, cache.intensity_profiles.ww
    );
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::cache_refresh_distance_profiles() const {
    auto pool = utility::multi_threading::get_global_pool();

    cache.distance_profiles.p_aa = Distribution1D(axis.bins, 0);
    cache.distance_profiles.p_aw = Distribution1D(axis.bins, 0);
    cache.distance_profiles.p_ww = Distribution1D(axis.bins, 0);

    pool->detach_task([this] () {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                std::transform(cache.distance_profiles.p_aa.begin(), cache.distance_profiles.p_aa.end(), distance_profiles.aa.begin(ff1, ff2), cache.distance_profiles.p_aa.begin(), std::plus<>());
            }
        }
    });
    pool->detach_task([this] () {
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            std::transform(cache.distance_profiles.p_aw.begin(), cache.distance_profiles.p_aw.end(), distance_profiles.aw.begin(ff1), cache.distance_profiles.p_aw.begin(), std::plus<>());
        }
    });
    pool->detach_task([this] () {
        std::transform(cache.distance_profiles.p_ww.begin(), cache.distance_profiles.p_ww.end(), distance_profiles.ww.begin(), cache.distance_profiles.p_ww.begin(), std::plus<>());
    });
    cache.distance_profiles.valid = true;
    pool->wait();
}

template<typename FormFactorTableType>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::cache_refresh_sinqd() const {
    auto pool = utility::multi_threading::get_global_pool();
    auto sinqd_table = get_sinc_table();

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if (cache.sinqd.aa.empty()) {
        cache.sinqd.aa = container::Container3D<double>(form_factor::get_count(), form_factor::get_count(), debye_axis.bins);
        cache.sinqd.ax = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        cache.sinqd.aw = container::Container2D<double>(form_factor::get_count(), debye_axis.bins);
        cache.sinqd.xx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.wx = container::Container1D<double>(debye_axis.bins);
        cache.sinqd.ww = container::Container1D<double>(debye_axis.bins);
    }

    std::vector<double> cx(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {cx[q-q0] = exv_factor(q);}
    for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
        for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
            pool->detach_task([this, q0, bins=debye_axis.bins, ff1, ff2, sinqd_table] () {
                for (unsigned int q = q0; q < q0+bins; ++q) {
                    cache.sinqd.aa.index(ff1, ff2, q-q0) = std::inner_product(distance_profiles.aa.begin(ff1, ff2), distance_profiles.aa.end(ff1, ff2), sinqd_table->begin(q), 0.0);
                }
            });
        }
        pool->detach_task([this, q0, bins=debye_axis.bins, ff1, sinqd_table] () {
            for (unsigned int q = q0; q < q0+bins; ++q) {
                cache.sinqd.ax.index(ff1, q-q0) = std::inner_product(distance_profiles.aa.begin(ff1, form_factor::exv_bin), distance_profiles.aa.end(ff1, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
                cache.sinqd.aw.index(ff1, q-q0) = std::inner_product(distance_profiles.aw.begin(ff1), distance_profiles.aw.end(ff1), sinqd_table->begin(q), 0.0);
            }
        });
    }
    pool->detach_task([&] () {
        for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
            cache.sinqd.xx.index(q-q0) = std::inner_product(distance_profiles.aa.begin(form_factor::exv_bin, form_factor::exv_bin), distance_profiles.aa.end(form_factor::exv_bin, form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            cache.sinqd.wx.index(q-q0) = std::inner_product(distance_profiles.aw.begin(form_factor::exv_bin), distance_profiles.aw.end(form_factor::exv_bin), sinqd_table->begin(q), 0.0);
            cache.sinqd.ww.index(q-q0) = std::inner_product(distance_profiles.ww.begin(), distance_profiles.ww.end(), sinqd_table->begin(q), 0.0);
        }
    });
    cache.sinqd.valid = true;
    pool->wait();
}

template<typename FormFactorTableType> template<bool sinqd_changed, bool cw_changed, bool cx_changed>
void CompositeDistanceHistogramFFAvgBase<FormFactorTableType>::cache_refresh_intensity_profiles() const {
    auto pool = utility::multi_threading::get_global_pool();
    const auto& ff_table = get_ff_table(); 

    Axis debye_axis = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax);
    unsigned int q0 = constants::axes::q_axis.get_bin(settings::axes::qmin);

    if constexpr (sinqd_changed) {
        cache.intensity_profiles.aa = std::vector<double>(debye_axis.bins, 0);
    }
    if constexpr (cw_changed) {
        cache.intensity_profiles.aw = std::vector<double>(debye_axis.bins, 0);
        cache.intensity_profiles.ww = std::vector<double>(debye_axis.bins, 0);
    }
    if constexpr (cx_changed) {
        cache.intensity_profiles.ax = std::vector<double>(debye_axis.bins, 0);
        cache.intensity_profiles.xx = std::vector<double>(debye_axis.bins, 0);
    }
    if constexpr (cw_changed || cx_changed) {
        cache.intensity_profiles.wx = std::vector<double>(debye_axis.bins, 0);
    }

    std::vector<double> cx(debye_axis.bins, 0);
    for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {cx[q-q0] = exv_factor(q);}

    if constexpr (sinqd_changed) {
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

    if constexpr (cx_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    cache.intensity_profiles.ax[q-q0] += 2*cx[q-q0]*cache.sinqd.ax.index(ff1, q-q0)*ff_table.index(ff1, form_factor::exv_bin).evaluate(q);
                }
            }
        });
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.xx[q-q0] += cx[q-q0]*cx[q-q0]*cache.sinqd.xx.index(q-q0)*ff_table.index(form_factor::exv_bin, form_factor::exv_bin).evaluate(q);
            }
        });
    }

    if constexpr (cw_changed) {
        pool->detach_task([&] () {
            for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
                for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                    cache.intensity_profiles.aw[q-q0] += 2*free_params.cw*cache.sinqd.aw.index(ff1, q-q0)*ff_table.index(ff1, form_factor::water_bin).evaluate(q);
                }
            }
        });
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.ww[q-q0] += free_params.cw*free_params.cw*cache.sinqd.ww.index(q-q0)*ff_table.index(form_factor::water_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }

    if constexpr (cw_changed || cx_changed) {
        pool->detach_task([&] () {
            for (unsigned int q = q0; q < q0+debye_axis.bins; ++q) {
                cache.intensity_profiles.wx[q-q0] += 2*cx[q-q0]*free_params.cw*cache.sinqd.wx.index(q-q0)*ff_table.index(form_factor::exv_bin, form_factor::water_bin).evaluate(q);
            }
        });
    }
    cache.intensity_profiles.cached_cx = free_params.cx;
    cache.intensity_profiles.cached_cw = free_params.cw;
    pool->wait();
}
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<true,  true,  true>()  const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<true,  true,  false>() const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<true,  false, true>()  const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<true,  false, false>() const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<false, true,  true>()  const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<false, true,  false>() const;
template void CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>::cache_refresh_intensity_profiles<false, false, true>()  const;
template class hist::CompositeDistanceHistogramFFAvgBase<form_factor::storage::atomic::table_t>;