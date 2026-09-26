// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMT.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <hist/detail/AtomOrdering.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool wb>
HistogramManagerMT<wb>::~HistogramManagerMT() = default;

template<bool wb>
std::unique_ptr<DistanceHistogram> HistogramManagerMT<wb>::calculate() {return calculate_all();}

template<bool wb>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMT<wb>::calculate_all() {
    logging::log("HistogramManagerMT::calculate: starting calculation");
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<wb>::type;

    // create a more compact representation of the coordinates
    // extremely wasteful to calculate this from scratch every time (class is not meant for serial use anyway?)
    auto data_a = hist::detail::factory::construct_from_atoms(this->protein);
    auto data_w = hist::detail::factory::construct_from_waters(this->protein);
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, this->protein);
    int bin_count = hist::detail::required_bin_count(data_a, data_w);
    hist::detail::decorrelate_order<wb>(bin_count, data_a, data_w);

    hist::distance_calculator::HistogramStore<wb> store(bin_count);
    int aa = store.allocate_1d(), ww = store.allocate_1d(), aw = store.allocate_1d();
    hist::distance_calculator::Calculator<wb> calculator(store);
    // all three are known up front, so they are held and dispatched as one unit
    calculator.hold();
    calculator.enqueue_calculate_self(data_a, aa);
    calculator.enqueue_calculate_self(data_w, ww);
    calculator.enqueue_calculate_cross(data_a, data_w, aw, 2);
    calculator.release_hold();
    calculator.run();

    auto p_aa = store.export_1d(aa);
    auto p_ww = store.export_1d(ww);
    auto p_aw = store.export_1d(aw);

    // calculate p_tot
    GenericDistribution1D_t p_tot(bin_count);
    for (int i = 0; i < p_tot.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + p_aw.index(i);}

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(p_tot);
    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (wb) {
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