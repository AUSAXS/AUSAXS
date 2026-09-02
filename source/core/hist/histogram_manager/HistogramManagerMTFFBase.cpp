// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFBase.h>
#include <hist/histogram_manager/detail/HistogramManagerMTFFHelpers.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/detail/BinEstimate.h>
#include <container/ThreadLocalWrapper.h>
#include <form_factor/FormFactorType.h>
#include <form_factor/lookup/FormFactorManager.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::container;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool wb, bool vbw>
HistogramManagerMTFFBase<wb, vbw>::~HistogramManagerMTFFBase() = default;

template<bool wb, bool vbw>
typename HistogramManagerMTFFBase<wb, vbw>::RawDistributions HistogramManagerMTFFBase<wb, vbw>::compute_raw_distributions() {
    assert(this->protein != nullptr && "HistogramManagerMTFFBase::compute_raw_distributions: Molecule is not set.");

    using GenericDistribution1D_t = typename GenericDistribution1D<wb>::type;
    using GenericDistribution2D_t = typename GenericDistribution2D<wb>::type;
    using GenericDistribution3D_t = typename GenericDistribution3D<wb>::type;
    auto pool = utility::multi_threading::get_global_pool();

    data_a_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_bodies());
    data_w_ptr = std::make_unique<CompactCoordinatesFF<vbw>>(this->protein->get_waters());
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();
    unsigned int bin_count = hist::detail::required_bin_count<vbw>(data_a, data_w);

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<GenericDistribution3D_t> p_aa_all(form_factor::get_active_count(), form_factor::get_active_count(), bin_count); // ff_type1, ff_type2, distance
    auto calc_aa = [&data_a, &p_aa_all, data_a_size] (int imin, int imax) {
        auto& p_aa = p_aa_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = i+1;                    // atom
            for (; j+15 < data_a_size; j+=16) {
                evaluate_aa16<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j+7 < data_a_size; j+=8) {
                evaluate_aa8<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate_aa4<vbw, 2>(p_aa, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate_aa1<vbw, 2>(p_aa, data_a, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution2D_t> p_aw_all(form_factor::get_active_count(), bin_count); // ff_type, distance
    auto calc_aw = [&data_w, &data_a, &p_aw_all, data_w_size] (int imin, int imax) {
        auto& p_aw = p_aw_all.get();
        for (int i = imin; i < imax; ++i) { // atom
            int j = 0;                      // water
            for (; j+15 < data_w_size; j+=16) {
                evaluate_aw16<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j+7 < data_w_size; j+=8) {
                evaluate_aw8<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate_aw4<vbw, 1>(p_aw, data_a, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate_aw1<vbw, 1>(p_aw, data_a, data_w, i, j);
            }
        }
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_ww_all(bin_count); // distance
    auto calc_ww = [&data_w, &p_ww_all, data_w_size] (int imin, int imax) {
        auto& p_ww = p_ww_all.get();
        for (int i = imin; i < imax; ++i) { // water
            int j = i+1;                    // water
            for (; j+15 < data_w_size; j+=16) {
                evaluate16<vbw, 2>(p_ww, data_w, data_w, i, j);
            }

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
    int job_size_a = settings::general::detail::get_job_size(data_a_size);
    int job_size_w = settings::general::detail::get_job_size(data_w_size);
    for (int i = 0; i < (int) data_a_size; i+=job_size_a) {
        pool->detach_task(
            [&calc_aa, i, job_size_a, data_a_size] () {calc_aa(i, std::min(i+job_size_a, data_a_size));}
        );
    }
    for (int i = 0; i < (int) data_a_size; i+=job_size_a) {
        pool->detach_task(
            [&calc_aw, i, job_size_a, data_a_size] () {calc_aw(i, std::min(i+job_size_a, data_a_size));}
        );
    }
    for (int i = 0; i < (int) data_w_size; i+=job_size_w) {
        pool->detach_task(
            [&calc_ww, i, job_size_w, data_w_size] () {calc_ww(i, std::min(i+job_size_w, data_w_size));}
        );
    }

    pool->wait();
    auto p_aa = p_aa_all.merge();
    auto p_aw = p_aw_all.merge();
    auto p_ww = p_ww_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    // save the self-correlations for later use in the intensity calculation
    for (int i = 0; i < data_a_size; ++i) {
        if constexpr (wb) {
            p_aa.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0, 0.0f);
        } else {
            p_aa.increment_index(data_a.get_ff_type(i), data_a.get_ff_type(i), 0);
        }
    }
    if constexpr (wb) {
        p_ww.add_index(0, WeightedEntry(data_w_size, data_w_size, 0));
    } else {
        p_ww.add_index(0, data_w_size);
    }

    GenericDistribution1D_t p_tot(bin_count);
    {   // sum all elements to the total
        unsigned int n_active = form_factor::get_active_count();
        for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            for (unsigned int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < n_active; ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<>());
            }
        }
        for (unsigned int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_aw.begin(ff1), p_tot.begin(), std::plus<>());
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
    pool->detach_task([&p_aw, max_bin] () { p_aw.resize(max_bin); });
    pool->detach_task([&p_ww, max_bin] () { p_ww.resize(max_bin); });
    pool->detach_task([&p_tot, max_bin] () { p_tot.resize(max_bin); });
    pool->wait();

    return RawDistributions{
        .p_aa = std::move(p_aa),
        .p_aw = std::move(p_aw),
        .p_ww = std::move(p_ww),
        .p_tot = std::move(p_tot),
        .max_bin = max_bin
    };
}

template class hist::HistogramManagerMTFFBase<false, false>;
template class hist::HistogramManagerMTFFBase<false, true>;
template class hist::HistogramManagerMTFFBase<true, false>;
template class hist::HistogramManagerMTFFBase<true, true>;
