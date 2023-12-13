#include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFExplicit.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFExplicit.h>
#include <hist/foxs/CompositeDistanceHistogramFoXS.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <form_factor/FormFactorType.h>
#include <data/Molecule.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <container/Container3D.h>
#include <container/Container2D.h>
#include <container/Container1D.h>
#include <constants/Constants.h>
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
    auto calc_pp = [&data_a, data_a_size] (int imin, int imax) {
        GenericDistribution3D_t p_aa(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
        GenericDistribution3D_t p_ax(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
        GenericDistribution3D_t p_xx(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
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
        return std::tuple(std::move(p_aa), std::move(p_ax), std::move(p_xx));
    };

    auto calc_hp = [&data_w, &data_a, data_w_size] (int imin, int imax) {
        GenericDistribution2D_t p_ww(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
        GenericDistribution2D_t p_wx(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
        for (int i = imin; i < imax; ++i) { // atom
            int j = 0;                      // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_ww, p_wx, data_a, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_ww, p_wx, data_a, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_ww, p_wx, data_a, data_w, i, j);
            }
        }
        return std::tuple(std::move(p_ww), std::move(p_wx));
    };

    auto calc_hh = [&data_w, data_w_size] (int imin, int imax) {
        GenericDistribution1D_t p_hh(constants::axes::d_axis.bins, 0);
        for (int i = imin; i < imax; ++i) { // water
            int j = i+1;                    // water
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_hh, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_hh, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_hh, data_w, data_w, i, j);
            }
        }
        return p_hh;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    int job_size = settings::general::detail::job_size;
    BS::multi_future<std::tuple<GenericDistribution3D_t, GenericDistribution3D_t, GenericDistribution3D_t>> pp;
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pp.push_back(pool->submit(calc_pp, i, std::min<int>(i+job_size, (int) data_a_size)));
    }
    BS::multi_future<std::tuple<GenericDistribution2D_t, GenericDistribution2D_t>> hp;
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        hp.push_back(pool->submit(calc_hp, i, std::min<int>(i+job_size, (int) data_a_size)));
    }
    BS::multi_future<GenericDistribution1D_t> hh;
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        hh.push_back(pool->submit(calc_hh, i, std::min<int>(i+job_size, (int) data_w_size)));
    }

    //#################//
    // COLLECT RESULTS //
    //#################//
    auto p_pp_future = pool->submit(
        [&]() {
            GenericDistribution3D_t p_aa(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            GenericDistribution3D_t p_ax(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            GenericDistribution3D_t p_xx(form_factor::get_count_without_excluded_volume(), form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type1, ff_type2, distance
            for (const auto& [p_aa_tmp, p_ax_tmp, p_xx_tmp] : pp.get()) {
                std::transform(p_aa.begin(), p_aa.end(), p_aa_tmp.begin(), p_aa.begin(), std::plus<>());
                std::transform(p_ax.begin(), p_ax.end(), p_ax_tmp.begin(), p_ax.begin(), std::plus<>());
                std::transform(p_xx.begin(), p_xx.end(), p_xx_tmp.begin(), p_xx.begin(), std::plus<>());
            }
            return std::tuple(std::move(p_aa), std::move(p_ax), std::move(p_xx));
        }
    );

    auto p_hp_future = pool->submit(
        [&]() {
            GenericDistribution2D_t p_wa(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
            GenericDistribution2D_t p_wx(form_factor::get_count_without_excluded_volume(), constants::axes::d_axis.bins, 0); // ff_type, distance
            for (const auto& [p_wa_tmp, p_wx_tmp] : hp.get()) {
                std::transform(p_wa.begin(), p_wa.end(), p_wa_tmp.begin(), p_wa.begin(), std::plus<>());
                std::transform(p_wx.begin(), p_wx.end(), p_wx_tmp.begin(), p_wx.begin(), std::plus<>());
            }
            return std::tuple(std::move(p_wa), std::move(p_wx));
        }
    );

    auto p_hh_future = pool->submit(
        [&]() {
            GenericDistribution1D_t p_ww(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hh.get()) {
                std::transform(p_ww.begin(), p_ww.end(), tmp.begin(), p_ww.begin(), std::plus<>());
            }
            return p_ww;
        }
    );
    pool->wait_for_tasks();

    // we cannot use structured bindings due to a Clang bug
    // auto [p_aa, p_ax, p_xx] = p_pp_future.get();
    // auto [p_wa, p_wx] = p_hp_future.get();
    // auto p_ww = p_hh_future.get();
    auto p3 = p_pp_future.get();
    auto p2 = p_hp_future.get();
    auto p_ww = p_hh_future.get();
    auto& p_aa = std::get<0>(p3);
    auto& p_ax = std::get<1>(p3);
    auto& p_xx = std::get<2>(p3);
    auto& p_wa = std::get<0>(p2);
    auto& p_wx = std::get<1>(p2);

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    for (int i = 0; i < data_a_size; ++i) {
        p_aa.index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0) += std::pow(data_a[i].value.w, 2);
        p_xx.index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0) += 1;
    }
    p_ww.index(0) = std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& data) {return sum + std::pow(data.value.w, 2);});

    // this is counter-intuitive, but splitting the loop into separate parts is likely faster since it allows both SIMD optimizations and better cache usage
    hist::Distribution1D p_tot(constants::axes::d_axis.bins, 0);
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

    pool->push_task([&] () { p_aa.resize(max_bin); });
    pool->push_task([&] () { p_ax.resize(max_bin); });
    pool->push_task([&] () { p_xx.resize(max_bin); });
    pool->push_task([&] () { p_wa.resize(max_bin); });
    pool->push_task([&] () { p_wx.resize(max_bin); });
    pool->push_task([&] () { p_ww.resize(max_bin); });
    pool->push_task([&] () { p_tot.resize(max_bin); });
    pool->wait_for_tasks();

    if (settings::hist::use_foxs_method) {
        return std::make_unique<CompositeDistanceHistogramFoXS>(
            std::move(p_aa), 
            std::move(p_ax), 
            std::move(p_xx),
            std::move(p_wa), 
            std::move(p_wx), 
            std::move(p_ww),
            Axis(0, max_bin*constants::axes::d_axis.width(), max_bin)
        );
    } else {
        return std::make_unique<CompositeDistanceHistogramFFExplicit>(
            std::move(p_aa), 
            std::move(p_ax), 
            std::move(p_xx),
            std::move(p_wa), 
            std::move(p_wx), 
            std::move(p_ww),
            Axis(0, max_bin*constants::axes::d_axis.width(), max_bin)
        );
    }
}

template class hist::HistogramManagerMTFFExplicit<false>;
template class hist::HistogramManagerMTFFExplicit<true>;