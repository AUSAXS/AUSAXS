// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTBase.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <hist/detail/AtomOrdering.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/CompactCoordinatesFactoryFF.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>

#include <algorithm>
#include <cassert>
#include <functional>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool wb, bool ff>
HistogramManagerMTBase<wb, ff>::~HistogramManagerMTBase() = default;

template<bool wb, bool ff>
typename HistogramManagerMTBase<wb, ff>::Distributions HistogramManagerMTBase<wb, ff>::compute_distributions() {
    assert(this->protein != nullptr && "HistogramManagerMTBase::compute_distributions: Molecule is not set.");
    using GenericDistribution1D_t = typename GenericDistribution1D<wb>::type;

    // the waters are a single set either way; with form factors their weights are simply never read
    data_w_ptr = std::make_unique<CompactCoordinates>(factory::construct_from_waters(this->protein));
    auto& data_w = *data_w_ptr;
    if constexpr (ff) {
        data_a_ptr = std::make_unique<std::vector<CompactCoordinates>>(factory::construct_by_ff_from_atoms(this->protein));
    } else {
        data_a_ptr = std::make_unique<CompactCoordinates>(factory::construct_from_atoms(this->protein));
        SimpleExvModel::apply_simple_excluded_volume(*data_a_ptr, this->protein);
    }
    auto& data_a = *data_a_ptr;
    int bin_count = hist::detail::required_bin_count(data_a, data_w);
    if constexpr (!ff) {hist::detail::decorrelate_order<wb>(bin_count, data_a, data_w);}

    // the shape of each result follows from its sets: with form factors the atoms are partitioned by type, the waters never are
    int classes = 1;
    if constexpr (ff) {classes = static_cast<int>(data_a.size());}
    hist::distance_calculator::HistogramStore<wb> store(bin_count, classes);
    int aa = ff ? store.allocate_3d() : store.allocate_1d();
    int aw = ff ? store.allocate_2d() : store.allocate_1d();
    int ww = store.allocate_1d();

    // with form factors, the pairs are only counted, and the form factors are applied later by the intensity calculator.
    // the form factor results store each atom-water pair once, while the weighted one stores both orders
    hist::distance_calculator::Calculator<wb, ff> calculator(store);
    constexpr int aw_pair_factor = ff ? 1 : 2;

    // the self-correlations are part of what the kernel evaluates, so they do not have to be added separately here.
    // all of them are known up front, so they are held and dispatched as one unit
    calculator.hold();
    calculator.enqueue_calculate_self(data_a, aa);
    calculator.enqueue_calculate_cross(data_a, data_w, aw, aw_pair_factor);
    calculator.enqueue_calculate_self(data_w, ww);
    calculator.release_hold();
    calculator.run();

    Distributions res;
    res.p_ww = store.export_1d(ww);
    res.p_tot = GenericDistribution1D_t(bin_count);
    auto add = [&p_tot = res.p_tot] (auto begin) {std::transform(p_tot.begin(), p_tot.end(), begin, p_tot.begin(), std::plus<>());};
    if constexpr (ff) {
        res.p_aa = store.export_3d(aa);
        res.p_aw = store.export_2d(aw);

        // no atom has the excluded volume type, so its rows are empty; they are skipped to make that explicit
        int n_ff = form_factor::get_active_count();
        for (int ff1 = form_factor::start_index_for_explicit_exv(); ff1 < n_ff; ++ff1) {
            for (int ff2 = ff1; ff2 < n_ff; ++ff2) {add(res.p_aa.begin(ff1, ff2));}
            add(res.p_aw.begin(ff1));
        }
    } else {
        res.p_aa = store.export_1d(aa);
        res.p_aw = store.export_1d(aw);
        add(res.p_aa.begin());
        add(res.p_aw.begin());
    }
    add(res.p_ww.begin());

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(res.p_tot);
    res.p_aa.resize(max_bin);
    res.p_aw.resize(max_bin);
    res.p_ww.resize(max_bin);
    res.p_tot.resize(max_bin);
    return res;
}

template class hist::HistogramManagerMTBase<false, false>;
template class hist::HistogramManagerMTBase<false, true>;
template class hist::HistogramManagerMTBase<true, false>;
template class hist::HistogramManagerMTBase<true, true>;
