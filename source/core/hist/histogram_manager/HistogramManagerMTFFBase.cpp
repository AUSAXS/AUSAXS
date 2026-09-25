// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFBase.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/CalculatorFF.h>
#include <hist/distance_calculator/HistogramStore.h>
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

    // the rows of the store: aa (ff1, ff2) from 0, aw (ff) from aw, and ww last
    int n_ff = form_factor::get_active_count();
    int aw = n_ff*n_ff, ww = aw + n_ff;
    hist::distance_calculator::HistogramStore<wb> store(ww + 1, bin_count);
    hist::distance_calculator::Calculator<wb, vbw> calculator(store);
    hist::distance_calculator::CalculatorFF<wb, vbw> calculator_ff(calculator);

    // the self-correlations are part of what the kernel evaluates, so they do not have to be added separately here.
    // all of them are known up front, so they are held and dispatched as one unit
    calculator.hold();
    calculator_ff.enqueue_self_by_ff(data_a, 0);
    calculator_ff.enqueue_cross_by_ff(data_a, data_w, aw);
    calculator_ff.enqueue_self_flat(data_w, ww);
    calculator.release_hold();
    calculator.run();

    auto p_aa = store.export_3d(0, n_ff, n_ff);
    auto p_aw = store.export_2d(aw, n_ff);
    auto p_ww = store.export_1d(ww);

    GenericDistribution1D_t p_tot(bin_count);
    {   // sum all elements to the total
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_ff; ++ff1) {
            for (int ff2 = form_factor::start_index_for_explicit_exv(); ff2 < n_ff; ++ff2) {
                std::transform(p_tot.begin(), p_tot.end(), p_aa.begin(ff1, ff2), p_tot.begin(), std::plus<>());
            }
        }
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_ff; ++ff1) {
            std::transform(p_tot.begin(), p_tot.end(), p_aw.begin(ff1), p_tot.begin(), std::plus<>());
        }
        std::transform(p_tot.begin(), p_tot.end(), p_ww.begin(), p_tot.begin(), std::plus<>());
    }

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(p_tot);

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
