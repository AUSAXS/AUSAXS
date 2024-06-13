/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/intensity_calculator/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/intensity_calculator/pepsi/CompositeDistanceHistogramPepsi.h>
#include <hist/intensity_calculator/crysol/CompositeDistanceHistogramCrysol.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <container/ThreadLocalWrapper.h>
#include <constants/Axes.h>
#include <utility/MultiThreading.h>

using namespace container;
using namespace hist;

template<bool use_weighted_distribution>
HistogramManagerMTFFExplicit<use_weighted_distribution>::~HistogramManagerMTFFExplicit() = default;

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFExplicit<use_weighted_distribution>::calculate() {return calculate_all();}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFExplicit<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    using GenericDistribution2D_t = typename hist::GenericDistribution2D<use_weighted_distribution>::type;
    using GenericDistribution3D_t = typename hist::GenericDistribution3D<use_weighted_distribution>::type;
    data_a_ptr = std::make_unique<hist::detail::CompactCoordinatesFF>(this->protein->get_bodies());
    data_w_ptr = std::make_unique<hist::detail::CompactCoordinatesFF>(this->protein->get_waters());
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    auto pool = utility::multi_threading::get_global_pool();

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_aa_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins
    ); // ff_type1, ff_type2, distance

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_ax_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins
    ); // ff_type1, ff_type2, distance

    container::ThreadLocalWrapper<GenericDistribution3D_t> p_xx_all(
        form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins
    ); // ff_type1, ff_type2, distance

    auto calc_aa = [&data_a, &p_aa_all, &p_ax_all, &p_xx_all, data_a_size] (int imin, int imax) {
        auto& p_aa = p_aa_all.get();
        auto& p_ax = p_ax_all.get();
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = i+1;                    // atom
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_aa, p_ax, p_xx, data_a, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_aa, p_ax, p_xx, data_a, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_aa, p_ax, p_xx, data_a, data_a, i, j);
             }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution2D_t> p_wa_all(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins); // ff_type, distance
    container::ThreadLocalWrapper<GenericDistribution2D_t> p_wx_all(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins); // ff_type, distance
    auto calc_wa = [&data_w, &data_a, &p_wa_all, &p_wx_all, data_w_size] (int imin, int imax) {
        auto& p_wa = p_wa_all.get();
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = 0;                      // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_wa, p_wx, data_a, data_w, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_ww_all(constants::axes::d_axis.bins); // distance
    auto calc_ww = [&data_w, &p_ww_all, data_w_size] (int imin, int imax) {
        auto& p_ww = p_ww_all.get();
        for (int i = imin; i < imax; ++i) { // water
            int j = i+1;                    // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
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
        p_aa.add(data_a.get_ff_type(i), data_a.get_ff_type(i), 0, std::pow(data_a[i].value.w, 2));
        p_xx.add(data_a.get_ff_type(i), data_a.get_ff_type(i), 0, 1);
    }
    p_ww.add(0, std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& data) {return sum + std::pow(data.value.w, 2);}));

    // this is counter-intuitive, but splitting the loop into separate parts is likely faster since it allows both SIMD optimizations and better cache usage
    GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
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

    switch (settings::hist::histogram_manager) {
        case settings::hist::HistogramManagerChoice::FoXSManager:
            return std::make_unique<CompositeDistanceHistogramFoXS>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot)
            );
        case settings::hist::HistogramManagerChoice::PepsiManager:
            return std::make_unique<CompositeDistanceHistogramPepsi>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot)
            );
        case settings::hist::HistogramManagerChoice::CrysolManager:
            return std::make_unique<CompositeDistanceHistogramCrysol>(
                std::move(Distribution3D(std::move(p_aa))), 
                std::move(Distribution3D(std::move(p_ax))), 
                std::move(Distribution3D(std::move(p_xx))),
                std::move(Distribution2D(std::move(p_wa))), 
                std::move(Distribution2D(std::move(p_wx))), 
                std::move(Distribution1D(std::move(p_ww))),
                std::move(p_tot)
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

template class hist::HistogramManagerMTFFExplicit<false>;
template class hist::HistogramManagerMTFFExplicit<true>;