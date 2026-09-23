// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFBase.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/CalculatorFF.h>
#include <utility/MultiThreading.h>

#include <cassert>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool wb, bool vbw>
HistogramManagerMTFFBase<wb, vbw>::~HistogramManagerMTFFBase() = default;

template<bool wb, bool vbw>
typename HistogramManagerMTFFBase<wb, vbw>::RawDistributions HistogramManagerMTFFBase<wb, vbw>::compute_raw_distributions() {
    assert(this->protein != nullptr && "HistogramManagerMTFFBase::compute_raw_distributions: Molecule is not set.");

    using GenericDistribution1D_t = typename GenericDistribution1D<wb>::type;
    auto* pool = utility::multi_threading::get_global_pool();

    data_a_ptr = std::make_unique<std::vector<CompactCoordinates<vbw>>>(hist::detail::factory::construct_unit_weight_by_ff_from_atoms<vbw>(this->protein));
    data_w_ptr = std::make_unique<CompactCoordinates<vbw>>(hist::detail::factory::construct_unit_weight_from_waters<vbw>(this->protein));
    auto& data_a = *data_a_ptr;
    auto& data_w = *data_w_ptr;
    int bin_count = hist::detail::required_bin_count<vbw>(data_a, data_w);

    // the self-correlations are part of what the kernel evaluates, so they do not have to be added separately here
    hist::distance_calculator::CalculatorFF<wb, vbw> calculator(bin_count);
    calculator.enqueue_self_by_ff(data_a);
    calculator.enqueue_cross_by_ff(data_a, data_w);
    calculator.enqueue_self_flat(data_w);
    auto res = calculator.run();

    auto p_aa = std::move(res.aa);
    auto p_aw = std::move(res.aw);
    auto p_ww = std::move(res.ww);

    GenericDistribution1D_t p_tot(bin_count);
    {   // sum all elements to the total
        int n_active = form_factor::get_active_count();
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            for (int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < n_active; ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<>());
            }
        }
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_active; ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_aw.begin(ff1), p_tot.begin(), std::plus<>());
        }
        std::transform(p_tot.begin(), p_tot.end(), p_ww.begin(), p_tot.begin(), std::plus<>());
    }

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
    for (int i = p_tot.size()-1; i >= 10; --i) {
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
