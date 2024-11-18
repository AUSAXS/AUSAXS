/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::hist;

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
    data_a.implicit_excluded_volume(this->protein->get_volume_grid()/data_a.size());

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<GenericDistribution1D_t> p_aa_all(constants::axes::d_axis.bins);
    auto calc_aa = [&data_a, &p_aa_all, data_a_size] (int imin, int imax) {
        auto& p_aa = p_aa_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = i+1;                    // atom
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_aw_all(constants::axes::d_axis.bins);
    auto calc_aw = [&data_w, &data_a, &p_aw_all, data_a_size] (int imin, int imax) {
        auto& p_aw = p_aw_all.get();
        for (int i = imin; i < imax; ++i) { // water
            int j = 0;                      // atom
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }
        }
    };
    
    container::ThreadLocalWrapper<GenericDistribution1D_t> p_ww_all(constants::axes::d_axis.bins);
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
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        pool->detach_task(
            [&calc_aw, i, job_size, data_w_size] () {calc_aw(i, std::min(i+job_size, data_w_size));}
        );
    }
    for (int i = 0; i < (int) data_w_size; i+=job_size) {
        pool->detach_task(
            [&calc_ww, i, job_size, data_w_size] () {calc_ww(i, std::min(i+job_size, data_w_size));}
        );
    }

    pool->wait();
    GenericDistribution1D_t p_aa = p_aa_all.merge();
    GenericDistribution1D_t p_aw = p_aw_all.merge();
    GenericDistribution1D_t p_ww = p_ww_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_aa.add(0, std::accumulate(data_a.get_data().begin(), data_a.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} ));
    p_ww.add(0, std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} ));

    // calculate p_tot
    GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
    for (unsigned int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + 2*p_aw.index(i);}

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (use_weighted_distribution) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(std::move(p_aa))), 
            std::move(Distribution1D(std::move(p_aw))), 
            std::move(Distribution1D(std::move(p_ww))), 
            std::move(p_tot)
        );
    } else {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(p_aa), 
            std::move(p_aw), 
            std::move(p_ww), 
            std::move(p_tot)
        );
    }
}

template class hist::HistogramManagerMT<false>;
template class hist::HistogramManagerMT<true>;