// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <grid/exv/RawGridExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/histogram_manager/detail/GridExvHelpers.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool variable_bin_width>
HistogramManagerMTFFGridScalableExv<variable_bin_width>::~HistogramManagerMTFFGridScalableExv() = default;

template<bool variable_bin_width>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGridScalableExv<variable_bin_width>::calculate() {
    return calculate_all();
}

template<bool variable_bin_width>
grid::exv::GridExcludedVolume HistogramManagerMTFFGridScalableExv<variable_bin_width>::get_exv() const {
    return grid::exv::RawGridExv::create(this->protein->get_grid());
}

template<bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGridScalableExv<variable_bin_width>::calculate_all() {
    logging::log("HistogramManagerMTFFGridScalableExv::calculate: starting calculation");
    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized

    // wrap all calculations into a lambda which we can later pass to the intensity calculator to allow it to rescale the excluded volume and easily reevaluate the histograms
    auto eval_scaled_exv = [
        atomic = hist::detail::grid_exv::AtomicDistributions::take(static_cast<CompositeDistanceHistogramFFAvg&>(*base_res)),
        data_a = *this->data_a_ptr,
        data_w = *this->data_w_ptr,
        data_x = hist::detail::factory::construct<variable_bin_width>(get_exv().interior)] 
        (double scale) 
    {
        // stretch the excluded volume cells by the given scale factor
        auto scaled_x = data_x;
        scaled_x.scale_coordinates(scale);
        int bin_count = hist::detail::required_bin_count<variable_bin_width>(data_a, data_w, scaled_x);

        //##############//
        // SUBMIT TASKS //
        //##############//
        distance_calculator::HistogramStore<true> store(bin_count, static_cast<int>(data_a.size()));
        int ax = store.allocate_2d(), wx = store.allocate_1d(), xx = store.allocate_1d();
        distance_calculator::Calculator<true, variable_bin_width, UNIT_WEIGHTS> calculator(store);
        calculator.enqueue_calculate_self(scaled_x, xx);
        calculator.enqueue_calculate_cross(data_a, scaled_x, ax, 1);
        calculator.enqueue_calculate_cross(data_w, scaled_x, wx, 1);

        calculator.run();
        WeightedDistribution1D p_xx_generic = store.export_1d(xx);
        WeightedDistribution2D p_ax_generic = store.export_2d(ax);
        WeightedDistribution1D p_wx_generic = store.export_1d(wx);

        // the excluded volume may reach further than the atoms, in which case the atomic distributions are grown to match
        auto scaled = atomic;
        scaled.grow(hist::detail::trimmed_bin_count(p_xx_generic, p_wx_generic));
        return hist::detail::grid_exv::splice(std::move(scaled), p_ax_generic, p_wx_generic, std::move(p_xx_generic));
    };

    return std::make_unique<CompositeDistanceHistogramFFGridScalableExv>(std::move(*eval_scaled_exv(1)), std::move(eval_scaled_exv));
}

template class ausaxs::hist::HistogramManagerMTFFGridScalableExv<true>;
template class ausaxs::hist::HistogramManagerMTFFGridScalableExv<false>;