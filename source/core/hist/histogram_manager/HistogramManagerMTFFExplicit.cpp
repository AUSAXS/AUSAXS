// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/detail/CompactCoordinates.h>
#include <form_factor/FormFactorType.h>
#include <data/Molecule.h>
#include <settings/ExvSettings.h>
#include <settings/GeneralSettings.h>
#include <container/ThreadLocalWrapper.h>
#include <utility/MultiThreading.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::container;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

namespace {
    // Local evaluation helpers for FFExplicit - atom-atom (3 x 3D histograms: aa, ax, xx)
    template<bool vbw, int factor>
    void evaluate_aa8(
        WeightedDistribution3D& p_aa, WeightedDistribution3D& p_ax, WeightedDistribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            int ff_j = data_a.get_ff_type(j+k);
            p_aa.increment<factor>(ff_i, ff_j, res.distances[k]);
            p_ax.increment<factor>(ff_i, ff_j, res.distances[k]);
            p_xx.increment<factor>(ff_i, ff_j, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_aa8(
        Distribution3D& p_aa, Distribution3D& p_ax, Distribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            int ff_j = data_a.get_ff_type(j+k);
            p_aa.increment_index<factor>(ff_i, ff_j, res.distances[k]);
            p_ax.increment_index<factor>(ff_i, ff_j, res.distances[k]);
            p_xx.increment_index<factor>(ff_i, ff_j, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_aa4(
        WeightedDistribution3D& p_aa, WeightedDistribution3D& p_ax, WeightedDistribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            int ff_j = data_a.get_ff_type(j+k);
            p_aa.increment<factor>(ff_i, ff_j, res.distances[k]);
            p_ax.increment<factor>(ff_i, ff_j, res.distances[k]);
            p_xx.increment<factor>(ff_i, ff_j, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_aa4(
        Distribution3D& p_aa, Distribution3D& p_ax, Distribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            int ff_j = data_a.get_ff_type(j+k);
            p_aa.increment_index<factor>(ff_i, ff_j, res.distances[k]);
            p_ax.increment_index<factor>(ff_i, ff_j, res.distances[k]);
            p_xx.increment_index<factor>(ff_i, ff_j, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_aa1(
        WeightedDistribution3D& p_aa, WeightedDistribution3D& p_ax, WeightedDistribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        int ff_j = data_a.get_ff_type(j);
        p_aa.increment<factor>(ff_i, ff_j, res.distance);
        p_ax.increment<factor>(ff_i, ff_j, res.distance);
        p_xx.increment<factor>(ff_i, ff_j, res.distance);
    }

    template<bool vbw, int factor>
    void evaluate_aa1(
        Distribution3D& p_aa, Distribution3D& p_ax, Distribution3D& p_xx,
        const CompactCoordinatesFF<vbw>& data_a, int i, int j
    ) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_a, data_a, i, j);
        int ff_i = data_a.get_ff_type(i);
        int ff_j = data_a.get_ff_type(j);
        p_aa.increment_index<factor>(ff_i, ff_j, res.distance);
        p_ax.increment_index<factor>(ff_i, ff_j, res.distance);
        p_xx.increment_index<factor>(ff_i, ff_j, res.distance);
    }

    // Local evaluation helpers for FFExplicit - atom-water (2 x 2D histograms: wa, wx)
    template<bool vbw, int factor>
    void evaluate_wa8(
        WeightedDistribution2D& p_wa, WeightedDistribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::OctoEvaluatedResult res = add8::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            p_wa.increment<factor>(ff_i, res.distances[k]);
            p_wx.increment<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_wa8(
        Distribution2D& p_wa, Distribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::OctoEvaluatedResultRounded res = add8::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 8; ++k) {
            p_wa.increment_index<factor>(ff_i, res.distances[k]);
            p_wx.increment_index<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_wa4(
        WeightedDistribution2D& p_wa, WeightedDistribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::QuadEvaluatedResult res = add4::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            p_wa.increment<factor>(ff_i, res.distances[k]);
            p_wx.increment<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_wa4(
        Distribution2D& p_wa, Distribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::QuadEvaluatedResultRounded res = add4::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        for (int k = 0; k < 4; ++k) {
            p_wa.increment_index<factor>(ff_i, res.distances[k]);
            p_wx.increment_index<factor>(ff_i, res.distances[k]);
        }
    }

    template<bool vbw, int factor>
    void evaluate_wa1(
        WeightedDistribution2D& p_wa, WeightedDistribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::EvaluatedResult res = add1::evaluate_weighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        p_wa.increment<factor>(ff_i, res.distance);
        p_wx.increment<factor>(ff_i, res.distance);
    }

    template<bool vbw, int factor>
    void evaluate_wa1(
        Distribution2D& p_wa, Distribution2D& p_wx,
        const CompactCoordinatesFF<vbw>& data_a, const CompactCoordinatesFF<vbw>& data_w, int i, int j
    ) {
        xyzff::EvaluatedResultRounded res = add1::evaluate_unweighted(data_a, data_w, i, j);
        int ff_i = data_a.get_ff_type(i);
        p_wa.increment_index<factor>(ff_i, res.distance);
        p_wx.increment_index<factor>(ff_i, res.distance);
    }

    // Note: evaluate_ww functions removed - use TemplateHelperAvg's evaluate8/4/1 directly
}

template<bool wb, bool vbw>
HistogramManagerMTFFExplicit<wb, vbw>::~HistogramManagerMTFFExplicit() = default;

template<bool wb, bool vbw>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit<wb, vbw>::calculate() {return calculate_all();}

template<bool wb, bool vbw>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFExplicit<wb, vbw>::calculate_all() {
    logging::log("HistogramManagerMTFFExplicit::calculate: starting calculation");
    using GenericDistribution1D_t = typename GenericDistribution1D<wb>::type;
    using GenericDistribution2D_t = typename GenericDistribution2D<wb>::type;
    using GenericDistribution3D_t = typename GenericDistribution3D<wb>::type;
    data_a_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_bodies());
    data_w_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_waters());
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    auto pool = utility::multi_threading::get_global_pool();

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_aa_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), settings::axes::bin_count
    ); // ff_type1, ff_type2, distance

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_ax_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), settings::axes::bin_count
    ); // ff_type1, ff_type2, distance

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_xx_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), settings::axes::bin_count
    ); // ff_type1, ff_type2, distance

    auto calc_aa = [&data_a, &p_aa_all, &p_ax_all, &p_xx_all, data_a_size] (int imin, int imax) {
        auto& p_aa = p_aa_all.get();
        auto& p_ax = p_ax_all.get();
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = i+1;                    // atom
            for (; j+7 < data_a_size; j+=8) {
                evaluate_aa8<vbw, 2>(p_aa, p_ax, p_xx, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate_aa4<vbw, 2>(p_aa, p_ax, p_xx, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate_aa1<vbw, 2>(p_aa, p_ax, p_xx, data_a, i, j);
             }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution2D_t> p_wa_all(form_factor::get_count_without_excluded_volume(), settings::axes::bin_count); // ff_type, distance
    container::ThreadLocalWrapper<GenericDistribution2D_t> p_wx_all(form_factor::get_count_without_excluded_volume(), settings::axes::bin_count); // ff_type, distance
    auto calc_wa = [&data_w, &data_a, &p_wa_all, &p_wx_all, data_w_size] (int imin, int imax) {
        auto& p_wa = p_wa_all.get();
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = 0;                      // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate_wa8<vbw, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate_wa4<vbw, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate_wa1<vbw, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_ww_all(settings::axes::bin_count); // distance
    auto calc_ww = [&data_w, &p_ww_all, data_w_size] (int imin, int imax) {
        auto& p_ww = p_ww_all.get();
        for (int i = imin; i < imax; ++i) { // water
            int j = i+1;                    // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<vbw, 2>(p_ww, data_w, data_w, i, j);
            }
        }
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    int job_size = settings::general::detail::job_size;
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pool->detach_task(
            [&calc_aa, i, job_size, data_a_size] () {calc_aa(i, std::min(i+job_size, data_a_size));}
        );
    }
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pool->detach_task(
            [&calc_wa, i, job_size, data_a_size] () {calc_wa(i, std::min(i+job_size, data_a_size));}
        ); 
    }
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        pool->detach_task(
            [&calc_ww, i, job_size, data_w_size] () {calc_ww(i, std::min(i+job_size, data_w_size));}
        );
    }

    pool->wait();
    auto p_aa = p_aa_all.merge();
    auto p_ax = p_ax_all.merge();
    auto p_xx = p_xx_all.merge();
    auto p_wa = p_wa_all.merge();
    auto p_wx = p_wx_all.merge();
    auto p_ww = p_ww_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    for (int i = 0; i < data_a_size; ++i) {
        p_aa.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0);
        p_xx.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0);
    }
    if constexpr (wb) {
        p_ww.add_index(0, WeightedEntry(data_w_size, data_w_size, 0));
    } else {
        p_ww.add_index(0, data_w_size);
    }

    // this is counter-intuitive, but splitting the loop into separate parts is likely faster since it allows both SIMD optimizations and better cache usage
    GenericDistribution1D_t p_tot(settings::axes::bin_count);
    {   // sum all elements to the total
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            for (unsigned int ff2 = 0; ff2 < form_factor::get_count_without_excluded_volume(); ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<>());
            }
        }
        for (unsigned int ff1 = 0; ff1 < form_factor::get_count_without_excluded_volume(); ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_wa.begin(ff1), p_tot.begin(), std::plus<>());
        }
        std::transform(p_tot.begin(), p_tot.end(), p_ww.begin(), p_tot.begin(), std::plus<>());
    }

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (unsigned int i = p_tot.size()-1; i >= 10; --i) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    pool->detach_task([&p_aa, max_bin] () { p_aa.resize(max_bin); });
    pool->detach_task([&p_ax, max_bin] () { p_ax.resize(max_bin); });
    pool->detach_task([&p_xx, max_bin] () { p_xx.resize(max_bin); });
    pool->detach_task([&p_wa, max_bin] () { p_wa.resize(max_bin); });
    pool->detach_task([&p_wx, max_bin] () { p_wx.resize(max_bin); });
    pool->detach_task([&p_ww, max_bin] () { p_ww.resize(max_bin); });
    pool->detach_task([&p_tot, max_bin] () { p_tot.resize(max_bin); });
    pool->wait();

    auto displaced_avg = [&] () {
        auto V = std::accumulate(
            data_a.get_data().begin(), 
            data_a.get_data().end(), 
            0.0, 
            [] (double sum, const auto& data_point) {
                auto ff = static_cast<form_factor::form_factor_t>(data_point.value.ff);
                switch (ff) {
                    case form_factor::form_factor_t::H: return constants::exv::standard.H + sum;
                    case form_factor::form_factor_t::C: return constants::exv::standard.C + sum;
                    case form_factor::form_factor_t::CH: return constants::exv::standard.C + sum;
                    case form_factor::form_factor_t::CH2: return constants::exv::standard.C + sum;
                    case form_factor::form_factor_t::CH3: return constants::exv::standard.C + sum;
                    case form_factor::form_factor_t::N: return constants::exv::standard.N + sum;
                    case form_factor::form_factor_t::NH: return constants::exv::standard.N + sum;
                    case form_factor::form_factor_t::NH2: return constants::exv::standard.N + sum;
                    case form_factor::form_factor_t::NH3: return constants::exv::standard.N + sum;
                    case form_factor::form_factor_t::O: return constants::exv::standard.O + sum;
                    case form_factor::form_factor_t::OH: return constants::exv::standard.O + sum;
                    case form_factor::form_factor_t::S: return constants::exv::standard.S + sum;
                    case form_factor::form_factor_t::SH: return constants::exv::standard.S + sum;
                    default: return sum + constants::exv::standard.OH;
                }
            }
        );
        return V / data_a_size;
    };

    switch (settings::exv::exv_method) {
        case settings::exv::ExvMethod::FoXS:
            return std::make_unique<CompositeDistanceHistogramFoXS>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot)
            );
        case settings::exv::ExvMethod::Pepsi:
            return std::make_unique<CompositeDistanceHistogramPepsi>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot),
                displaced_avg()
            );
        case settings::exv::ExvMethod::CRYSOL:
            return std::make_unique<CompositeDistanceHistogramCrysol>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot),
                displaced_avg()
            );
        default:
            return std::make_unique<CompositeDistanceHistogramFFExplicit>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot)
            );
    }
}

template class hist::HistogramManagerMTFFExplicit<false, false>;
template class hist::HistogramManagerMTFFExplicit<false, true>;
template class hist::HistogramManagerMTFFExplicit<true, false>;
template class hist::HistogramManagerMTFFExplicit<true, true>;