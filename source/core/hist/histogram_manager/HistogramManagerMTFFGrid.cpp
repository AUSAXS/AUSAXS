// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <grid/exv/RawGridExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/GridExvHelpers.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>  // IWYU pragma: keep
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

HistogramManagerMTFFGrid::~HistogramManagerMTFFGrid() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGrid::calculate() {
    return calculate_all();
}

ausaxs::grid::exv::GridExcludedVolume HistogramManagerMTFFGrid::get_exv() const {
    return grid::exv::RawGridExv::create(this->protein->get_grid());
}

std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGrid::calculate_all() {
    logging::log("HistogramManagerMTFFGrid::calculate: starting calculation");

    auto base_res = HistogramManagerMTFFAvg<true>::calculate_all(); // make sure everything is initialized
    auto exv = get_exv();
    auto data_x = hist::detail::factory::construct(exv.interior);
    const auto& data_a = *this->data_a_ptr;
    const auto& data_w = *this->data_w_ptr;
    int bin_count = hist::detail::required_bin_count(data_a, data_w, data_x);

    //##############//
    // SUBMIT TASKS //
    //##############//
    // the atoms are partitioned by form factor, the waters and excluded volume cells are not
    distance_calculator::HistogramStore<true> store(bin_count, static_cast<int>(data_a.size()));
    int ax = store.allocate_2d(), wx = store.allocate_1d();
#if !defined(POCKETFFT_AVAILABLE)
    int xx = store.allocate_1d();
#endif
    distance_calculator::Calculator<true, UNIT_WEIGHTS> calculator(store);
    calculator.hold();
    calculator.enqueue_calculate_cross(data_a, data_x, ax, 1);
    calculator.enqueue_calculate_cross(data_w, data_x, wx, 1);
    calculator.release_hold();

#if defined(POCKETFFT_AVAILABLE)
    // use the more efficient lattice transform for the self-correlation. it runs on the calling thread, overlapping with the jobs above.
    WeightedDistribution1D p_xx_generic = detail::lattice::self_correlation(
        exv, detail::inv_bin_width(), bin_count
    );
    p_xx_generic.add_index(0, detail::WeightedEntry(data_x.size(), data_x.size(), 0)); // self-correlations
    calculator.run();
#else
    calculator.enqueue_calculate_self(data_x, xx);
    calculator.run();
    WeightedDistribution1D p_xx_generic = store.export_1d(xx);
#endif
    WeightedDistribution2D p_ax_generic = store.export_2d(ax);
    WeightedDistribution1D p_wx_generic = store.export_1d(wx);

    // the excluded volume may reach further than the atoms, in which case the atomic distributions are grown to match
    auto atomic = grid_exv::AtomicDistributions::take(static_cast<CompositeDistanceHistogramFFAvg&>(*base_res));
    atomic.grow(hist::detail::trimmed_bin_count(p_xx_generic, p_wx_generic));
    return grid_exv::splice(std::move(atomic), p_ax_generic, p_wx_generic, std::move(p_xx_generic));
}