// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <grid/exv/RawGridExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/CalculatorFF.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool variable_bin_width>
HistogramManagerMTFFGrid<variable_bin_width>::~HistogramManagerMTFFGrid() = default;

template<bool variable_bin_width>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGrid<variable_bin_width>::calculate() {
    return calculate_all();
}

template<bool variable_bin_width>
ausaxs::grid::exv::GridExcludedVolume HistogramManagerMTFFGrid<variable_bin_width>::get_exv() const {
    return grid::exv::RawGridExv::create(this->protein->get_grid());
}

template<bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGrid<variable_bin_width>::calculate_all() {
    logging::log("HistogramManagerMTFFGrid::calculate: starting calculation");

    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized
    auto exv = get_exv();
    auto data_x = hist::detail::factory::construct_unit_weight<variable_bin_width>(exv.interior);
    const auto& data_a = *this->data_a_ptr;
    const auto& data_w = *this->data_w_ptr;
    int bin_count = hist::detail::required_bin_count<variable_bin_width>(data_a, data_w, data_x);

    //##############//
    // SUBMIT TASKS //
    //##############//
    // the rows of the store: ax (ff) from 0, then wx and xx. the atoms are resolved by form factor on their own side only
    int n_ff = form_factor::get_active_count();
    int wx = n_ff, xx = n_ff + 1;
    distance_calculator::HistogramStore<true> store(xx + 1, bin_count);
    distance_calculator::Calculator<true, variable_bin_width> calculator(store);
    calculator.hold();
    distance_calculator::CalculatorFF<true, variable_bin_width>(calculator).enqueue_cross_by_ff(data_a, data_x, 0);
    calculator.enqueue_calculate_cross(data_w, data_x, wx, 1);
    calculator.release_hold();

#if defined(POCKETFFT_AVAILABLE)
    // use the more efficient lattice transform for the self-correlation. it runs on the calling thread, overlapping with the jobs above.
    WeightedDistribution1D p_xx_generic = detail::lattice::self_correlation(
        exv, detail::WidthController<variable_bin_width>::get_inv_width(), bin_count
    );
    p_xx_generic.add_index(0, detail::WeightedEntry(data_x.size(), data_x.size(), 0)); // self-correlations
    calculator.run();
#else
    calculator.enqueue_calculate_self(data_x, xx);
    calculator.run();
    WeightedDistribution1D p_xx_generic = store.export_1d(xx);
#endif
    WeightedDistribution2D p_ax_generic = store.export_2d(0, n_ff);
    WeightedDistribution1D p_wx_generic = store.export_1d(wx);

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(p_xx_generic, p_wx_generic);

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto* cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    WeightedDistribution1D p_tot = cast_res->get_weighted_counts();
    p_tot.set_bin_centers(cast_res->get_d_axis());

    Distribution3D p_aa = std::move(cast_res->get_raw_aa_counts_by_ff());
    Distribution2D p_aw = std::move(cast_res->get_raw_aw_counts_by_ff());
    Distribution1D p_ww = std::move(cast_res->get_raw_ww_counts_by_ff());

    // either xx or ww are largest of all components
    max_bin = std::max<int>(max_bin, p_tot.size());

    // downsize the axes to only the relevant area
    if (static_cast<int>(base_res->get_d_axis().size()) < max_bin) {
        p_aa.resize(max_bin);
        p_aw.resize(max_bin);
        p_ww.resize(max_bin);
    } else {
        max_bin = base_res->get_d_axis().size(); // make sure we overwrite anything which may already be stored
    }

    // calculate weighted distance bins
    p_tot.resize(max_bin);
    WeightedDistribution1D p_tot_ax = std::max<int>(max_bin, p_wx_generic.size());
    for (int i = 0; i < max_bin; ++i) {
        p_tot_ax.add_index(i, p_wx_generic.index(i));
    }

    for (int i = 0; i < p_ax_generic.size_x(); ++i) {
        for (int j = 0; j < max_bin; ++j) {
            p_tot_ax.add_index(j, p_ax_generic.index(i, j));
        }
    }

    // overwrite the excluded volume calculations from the HistogramManagerMTFFAvg calculations with our new grid-based ones
    // first cast the weighted distributions to make iteration simpler
    Distribution2D p_ax(p_ax_generic);
    Distribution1D p_wx(p_wx_generic);
    Distribution1D p_xx(p_xx_generic);

    // replace the calculations
    for (int i = 0; i < p_aa.size_x(); ++i) {
        std::move(p_ax.begin(i), p_ax.begin(i)+max_bin, p_aa.begin(i, form_factor::exv_bin));
    }
    std::move(p_wx.begin(), p_wx.begin()+max_bin, p_aw.begin(form_factor::exv_bin));
    std::move(p_xx.begin(), p_xx.begin()+max_bin, p_aa.begin(form_factor::exv_bin, form_factor::exv_bin));

    return std::make_unique<CompositeDistanceHistogramFFGrid>(
        std::move(p_aa), 
        std::move(p_aw), 
        std::move(p_ww), 
        std::move(p_tot),
        std::move(p_tot_ax),
        std::move(p_xx_generic)
    );
}

template class ausaxs::hist::HistogramManagerMTFFGrid<false>;
template class ausaxs::hist::HistogramManagerMTFFGrid<true>;