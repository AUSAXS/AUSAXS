// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>

#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>  // IWYU pragma: keep
#include <form_factor/FormFactorType.h>
#include <grid/exv/RawGridExv.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/distance_calculator/detail/CPUKernel.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::hist;
namespace kernel = ausaxs::hist::distance_calculator::detail;

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
    auto* pool = utility::multi_threading::get_global_pool();
    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto* cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    WeightedDistribution1D p_tot = cast_res->get_weighted_counts();
    p_tot.set_bin_centers(cast_res->get_d_axis());

    // wrap all calculations into a lambda which we can later pass to the intensity calculator to allow it to rescale the excluded volume and easily reevaluate the histograms
    // the atoms are resolved by form factor on their own side only
    int n_ff = form_factor::get_active_count();
    auto eval_scaled_exv = [
        p_tot = std::move(p_tot),
        p_aa = std::move(cast_res->get_raw_aa_counts_by_ff()),
        p_aw = std::move(cast_res->get_raw_aw_counts_by_ff()),
        p_ww = std::move(cast_res->get_raw_ww_counts_by_ff()),
        data_a = *this->data_a_ptr,
        data_w = *this->data_w_ptr,
        data_x = hist::detail::factory::construct_unit_weight<variable_bin_width>(get_exv().interior),
        n_ff,
        pool] 
        (double scale) 
    {
        // stretch the excluded volume cells by the given scale factor
        auto scaled_x = data_x;
        scaled_x.scale_coordinates(scale);
        int bin_count = hist::detail::required_bin_count<variable_bin_width>(data_a, data_w, scaled_x);

        //##############//
        // SUBMIT TASKS //
        //##############//
        container::ThreadLocalWrapper<WeightedDistribution1D> p_xx_all(bin_count);
        container::ThreadLocalWrapper<WeightedDistribution2D> p_ax_all(n_ff, bin_count);
        container::ThreadLocalWrapper<WeightedDistribution1D> p_wx_all(bin_count);
        kernel::enqueue_self<true, variable_bin_width, 2, 1>(scaled_x, kernel::row_target(p_xx_all, bin_count));
        for (int ff = 0; ff < n_ff; ++ff) {
            kernel::enqueue_balanced_cross<variable_bin_width, 1>(data_a[ff], scaled_x, kernel::row_target(p_ax_all, bin_count, ff));
        }
        kernel::enqueue_balanced_cross<variable_bin_width, 1>(data_w, scaled_x, kernel::row_target(p_wx_all, bin_count));

        pool->wait();
        WeightedDistribution1D p_xx_generic = p_xx_all.merge();
        WeightedDistribution2D p_ax_generic = p_ax_all.merge();
        WeightedDistribution1D p_wx_generic = p_wx_all.merge();

        // downsize our axes to only the relevant area
        int max_bin = 10; // minimum size is 10
        for (int i = p_xx_generic.size()-1; i >= 10; --i) {
            if (p_xx_generic.index(i) != 0 || p_wx_generic.index(i) != 0) {
                max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
                break;
            }
        }

        // downsize the axes to only the relevant area
        auto new_p_aa = p_aa;
        auto new_p_aw = p_aw;
        auto new_p_ww = p_ww;
        if (p_tot.size() < max_bin) {
            new_p_aa.resize(max_bin);
            new_p_aw.resize(max_bin);
            new_p_ww.resize(max_bin);
        } else {
            max_bin = p_tot.size(); // make sure we overwrite anything which may already be stored
        }

        // calculate weighted distance bins
        auto new_p_tot = p_tot;
        new_p_tot.resize(max_bin);
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
        Distribution2D p_ax = Distribution2D(p_ax_generic);
        Distribution1D p_wx = Distribution1D(p_wx_generic);
        Distribution1D p_xx = Distribution1D(p_xx_generic);

        // replace the calculations
        for (int i = 0; i < p_aa.size_x(); ++i) {
            std::move(p_ax.begin(i), p_ax.begin(i)+max_bin, new_p_aa.begin(i, form_factor::exv_bin));
        }
        std::move(p_wx.begin(), p_wx.begin()+max_bin, new_p_aw.begin(form_factor::exv_bin));
        std::move(p_xx.begin(), p_xx.begin()+max_bin, new_p_aa.begin(form_factor::exv_bin, form_factor::exv_bin));

        return std::make_unique<CompositeDistanceHistogramFFGrid>(
            std::move(new_p_aa), 
            std::move(new_p_aw), 
            std::move(new_p_ww), 
            std::move(new_p_tot),
            std::move(p_tot_ax),
            std::move(p_xx_generic)
        );
    };

    return std::make_unique<CompositeDistanceHistogramFFGridScalableExv>(std::move(*eval_scaled_exv(1)), std::move(eval_scaled_exv));
}

template class ausaxs::hist::HistogramManagerMTFFGridScalableExv<true>;
template class ausaxs::hist::HistogramManagerMTFFGridScalableExv<false>;