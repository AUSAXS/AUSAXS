// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>

#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/GridExvFFT.h>
#include <hist/distance_calculator/Calculator.h>
#include <hist/distance_calculator/CalculatorFF.h>
#include <hist/distance_calculator/HistogramStore.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool variable_bin_width>
HistogramManagerMTFFGridSurface<variable_bin_width>::~HistogramManagerMTFFGridSurface() = default;

template<bool variable_bin_width>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGridSurface<variable_bin_width>::calculate() {
    return calculate_all();
}

template<bool variable_bin_width>
grid::exv::GridExcludedVolume HistogramManagerMTFFGridSurface<variable_bin_width>::get_exv() const {
    return grid::exv::RawGridWithSurfaceExv::create(this->protein->get_grid());
}

template<bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGridSurface<variable_bin_width>::calculate_all() {
    logging::log("HistogramManagerMTFFGridSurface::calculate: starting calculation");
    using XXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::XXContainer;
    using AXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::AXContainer;
    using WXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::WXContainer;

    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized
    auto exv = get_exv();
    auto data_x_i = hist::detail::factory::construct_unit_weight<variable_bin_width>(exv.interior);
    auto data_x_s = hist::detail::factory::construct_unit_weight<variable_bin_width>(exv.surface);
    const auto& data_a = *this->data_a_ptr;
    const auto& data_w = *this->data_w_ptr;
    int bin_count = hist::detail::required_bin_count<variable_bin_width>(data_a, data_w, data_x_i, data_x_s);

    //##############//
    // SUBMIT TASKS //
    //##############//
    // the rows of the store: ax (ff) of the interior from 0 and of the surface from n_ff, then wx of the interior and
    // surface, then xx of the interior, the surface, and their cross term. the atoms are resolved by form factor on their
    // own side only
    int n_ff = form_factor::get_active_count();
    int ax_i = 0, ax_s = n_ff, wx_i = 2*n_ff, wx_s = wx_i + 1, xx = wx_i + 2;
    distance_calculator::HistogramStore<true> store(xx + 3, bin_count);
    distance_calculator::Calculator<true, variable_bin_width> calculator(store);
    distance_calculator::CalculatorFF<true, variable_bin_width> calculator_ff(calculator);
    calculator.hold();
    calculator_ff.enqueue_cross_by_ff(data_a, data_x_i, ax_i);
    calculator_ff.enqueue_cross_by_ff(data_a, data_x_s, ax_s);
    calculator.enqueue_calculate_cross(data_w, data_x_i, wx_i, 1);
    calculator.enqueue_calculate_cross(data_w, data_x_s, wx_s, 1);
    calculator.release_hold();

    XXContainer p_xx(0);
#if defined(POCKETFFT_AVAILABLE)
    // use the more efficient lattice transform for the self-correlation. it runs on the calling thread, overlapping with the jobs above.
    auto p_xx_lattice = detail::lattice::correlations(
        exv, hist::detail::WidthController<variable_bin_width>::get_inv_width(), bin_count
    );
    p_xx.interior = std::move(p_xx_lattice.first);
    p_xx.surface  = std::move(p_xx_lattice.second);
    p_xx.cross    = std::move(p_xx_lattice.cross);
    p_xx.interior.add_index(0, detail::WeightedEntry(data_x_i.size(), data_x_i.size(), 0)); // self-correlations
    p_xx.surface.add_index(0, detail::WeightedEntry(data_x_s.size(), data_x_s.size(), 0));  // self-correlations
    calculator.run();
#else
    calculator.enqueue_calculate_self(data_x_i, xx);
    calculator.enqueue_calculate_self(data_x_s, xx + 1);
    calculator.enqueue_calculate_cross(data_x_i, data_x_s, xx + 2, 2);
    calculator.run();
    p_xx.interior = store.export_1d(xx);
    p_xx.surface  = store.export_1d(xx + 1);
    p_xx.cross    = store.export_1d(xx + 2);
#endif
    AXContainer p_ax(0, 0);
    p_ax.interior = store.export_2d(ax_i, n_ff);
    p_ax.surface  = store.export_2d(ax_s, n_ff);
    WXContainer p_wx(0);
    p_wx.interior = store.export_1d(wx_i);
    p_wx.surface  = store.export_1d(wx_s);

    // downsize our axes to only the relevant area
    int max_bin = hist::detail::trimmed_bin_count(p_xx.surface, p_xx.interior, p_wx.surface, p_wx.interior);

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
    WeightedDistribution1D p_tot_ax = std::max<int>(max_bin, p_wx.surface.size());
    for (int i = 0; i < max_bin; ++i) {
        p_tot_ax.add_index(i, p_wx.interior.index(i));
        p_tot_ax.add_index(i, p_wx.surface.index(i));
    }

    for (int i = 0; i < p_ax.surface.size_x(); ++i) {
        for (int j = 0; j < max_bin; ++j) {
            p_tot_ax.add_index(j, p_ax.interior.index(i, j));
            p_tot_ax.add_index(j, p_ax.surface.index(i, j));
        }
    }

    WeightedDistribution1D p_tot_xx = std::max<int>(max_bin, p_xx.surface.size());
    for (int i = 0; i < max_bin; ++i) {
        p_tot_xx.add_index(i, p_xx.interior.index(i));
        p_tot_xx.add_index(i, p_xx.surface.index(i));
        p_tot_xx.add_index(i, p_xx.cross.index(i));
    }

    {   // delete the exv information from the HistogramManagerMTFFAvg data
        // we delegate this work to the DistanceHistogram class, since it must be able to do this anyway to vary the surface contribution
        for (int i = 0; i < p_aa.size_x(); ++i) {
            std::for_each(p_aa.begin(i, form_factor::exv_bin), p_aa.end(i, form_factor::exv_bin) , [](auto& x) {x = 0;});
        }
        std::for_each(p_aw.begin(form_factor::exv_bin), p_aw.end(form_factor::exv_bin), [](auto& x) {x = 0;});
        std::for_each(p_aa.begin(form_factor::exv_bin, form_factor::exv_bin), p_aa.end(form_factor::exv_bin, form_factor::exv_bin), [](auto& x) {x = 0;});
    }

    return std::make_unique<CompositeDistanceHistogramFFGridSurface>(
        std::move(p_aa), 
        std::move(p_aw), 
        std::move(p_ww), 
        std::move(p_xx),
        std::move(p_ax),
        std::move(p_wx),
        std::move(p_tot),
        std::move(p_tot_ax),
        std::move(p_tot_xx)
    );
}

template class ausaxs::hist::HistogramManagerMTFFGridSurface<true>;
template class ausaxs::hist::HistogramManagerMTFFGridSurface<false>;