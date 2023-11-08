#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <data/Molecule.h>
#include <settings/GeneralSettings.h>
#include <constants/Constants.h>
#include <utility/MultiThreading.h>

using namespace hist;

template<bool use_weighted_distribution>
HistogramManagerMT<use_weighted_distribution>::~HistogramManagerMT() = default;

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> HistogramManagerMT<use_weighted_distribution>::calculate() {return calculate_all();}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMT<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time (class is not meant for serial use anyway?)
    data_a_ptr = std::make_unique<hist::detail::CompactCoordinates>(this->protein->get_bodies());
    data_w_ptr = std::make_unique<hist::detail::CompactCoordinates>(this->protein->get_waters());
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    auto calc_pp = [&data_a, data_a_size] (int imin, int imax) {
        GenericDistribution1D_t p_pp(constants::axes::d_axis.bins, 0);
        for (int i = imin; i < imax; ++i) {
            int j = i+1;
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_pp, data_a, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_pp, data_a, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_pp, data_a, data_a, i, j);
            }
        }
        return p_pp;
    };

    auto calc_hh = [&data_w, data_w_size] (int imin, int imax) {
        GenericDistribution1D_t p_hh(constants::axes::d_axis.bins, 0);
        for (int i = imin; i < imax; ++i) {
            int j = i+1;
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

    auto calc_hp = [&data_w, &data_a, data_a_size] (int imin, int imax) {
        GenericDistribution1D_t p_hp(constants::axes::d_axis.bins, 0);
        for (int i = imin; i < imax; ++i) {
            int j = 0;
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_hp, data_w, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_hp, data_w, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_hp, data_w, data_a, i, j);
            }
        }
        return p_hp;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    int job_size = settings::general::detail::job_size;
    BS::multi_future<GenericDistribution1D_t> pp;
    for (int i = 0; i < (int) data_a_size; i+=job_size) {
        pp.push_back(pool->submit(calc_pp, i, std::min<int>(i+job_size, (int) data_a_size)));
    }
    BS::multi_future<GenericDistribution1D_t> hh;
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        hh.push_back(pool->submit(calc_hh, i, std::min<int>(i+job_size, (int) data_w_size)));
    }
    BS::multi_future<GenericDistribution1D_t> hp;
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        hp.push_back(pool->submit(calc_hp, i, std::min<int>(i+job_size, (int) data_w_size)));
    }

    //#################//
    // COLLECT RESULTS //
    //#################//
    auto p_pp_future = pool->submit(
        [&](){
            GenericDistribution1D_t p_pp(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : pp.get()) {
                std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<>());
            }
            return p_pp;
        }
    );

    auto p_hp_future = pool->submit(
        [&](){
            GenericDistribution1D_t p_hp(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hp.get()) {
                std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<>());
            }
            return p_hp;
        }
    );

    auto p_hh_future = pool->submit(
        [&](){
            GenericDistribution1D_t p_hh(constants::axes::d_axis.bins, 0);
            for (const auto& tmp : hh.get()) {
                std::transform(p_hh.begin(), p_hh.end(), tmp.begin(), p_hh.begin(), std::plus<>());
            }
            return p_hh;            
        }
    );

    pool->wait_for_tasks();
    auto p_pp = p_pp_future.get();
    auto p_hp = p_hp_future.get();
    auto p_hh = p_hh_future.get();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_pp.index(0) = std::accumulate(data_a.get_data().begin(), data_a.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );
    p_hh.index(0) = std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );

    // calculate p_tot    
    Distribution1D p_tot(constants::axes::d_axis.bins, 0);
    for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_pp.index(i) + p_hh.index(i) + 2*p_hp.index(i);}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_pp.resize(max_bin);
    p_hh.resize(max_bin);
    p_hp.resize(max_bin);
    p_tot.resize(max_bin);
    return std::make_unique<CompositeDistanceHistogram>(
        std::move(p_pp), 
        std::move(p_hp), 
        std::move(p_hh), 
        std::move(p_tot), 
        Axis(0, max_bin*constants::axes::d_axis.width(), max_bin)
    );
}

template class hist::HistogramManagerMT<false>;
template class hist::HistogramManagerMT<true>;