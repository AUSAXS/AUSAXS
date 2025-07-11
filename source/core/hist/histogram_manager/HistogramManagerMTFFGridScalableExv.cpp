// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridScalableExv.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <grid/exv/RawGridExv.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFAvg.h>
#include <form_factor/FormFactorType.h>
#include <utility/MultiThreading.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;

// custom evaluates for the grid since we don't want to account for the excluded volume
namespace ausaxs::grid {
    template<bool use_weighted_distribution, int factor>
    inline void evaluate8(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ausaxs::detail::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool use_weighted_distribution, int factor>
    inline void evaluate4(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ausaxs::detail::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool use_weighted_distribution, int factor>
    inline void evaluate1(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ausaxs::detail::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        p.add(data_i.get_ff_type(i), res.distance, factor*res.weight);
    }
}

HistogramManagerMTFFGridScalableExv::~HistogramManagerMTFFGridScalableExv() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGridScalableExv::calculate() {
    return calculate_all();
}

grid::exv::GridExcludedVolume HistogramManagerMTFFGridScalableExv::get_exv() const {
    return grid::exv::RawGridExv::create(this->protein->get_grid());
}

std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGridScalableExv::calculate_all() {
    logging::log("HistogramManagerMTFFGridScalableExv::calculate: starting calculation");
    auto pool = utility::multi_threading::get_global_pool();
    auto base_res = HistogramManagerMTFFAvg<true>::calculate_all(); // make sure everything is initialized

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    WeightedDistribution1D p_tot = std::move(cast_res->get_counts());
    p_tot.set_bin_centers(cast_res->get_d_axis());

    // wrap all calculations into a lambda which we can later pass to the intensity calculator to allow it to rescale the excluded volume and easily reevaluate the histograms
    auto eval_scaled_exv = [
        p_tot = std::move(p_tot),
        p_aa = std::move(cast_res->get_aa_counts_ff()),
        p_aw = std::move(cast_res->get_aw_counts_ff()),
        p_ww = std::move(cast_res->get_ww_counts_ff()),
        data_a = *this->data_a_ptr, 
        data_w = *this->data_w_ptr, 
        data_x = hist::detail::CompactCoordinates(std::move(get_exv().interior), 1),
        pool] 
        (double scale) 
    {
        int data_a_size = (int) data_a.size();
        int data_w_size = (int) data_w.size();
        int data_x_size = (int) data_x.size();

        // stretch the excluded volume cells by the given scale factor
        auto scaled_data_x = data_x;
        for (auto& coord : scaled_data_x.get_data()) {
            coord.value.pos.x() *= scale;
            coord.value.pos.y() *= scale;
            coord.value.pos.z() *= scale;
        }

        //########################//
        // PREPARE MULTITHREADING //
        //########################//
        container::ThreadLocalWrapper<WeightedDistribution1D> p_xx_all(constants::axes::d_axis.bins);
        auto calc_xx = [&scaled_data_x, &p_xx_all, data_x_size] (int imin, int imax) {
            auto& p_xx = p_xx_all.get();
            for (int i = imin; i < imax; ++i) { // exv
                int j = i+1;                    // exv
                for (; j+7 < data_x_size; j+=8) {
                    evaluate8<true, 2>(p_xx, scaled_data_x, scaled_data_x, i, j);
                }

                for (; j+3 < data_x_size; j+=4) {
                    evaluate4<true, 2>(p_xx, scaled_data_x, scaled_data_x, i, j);
                }

                for (; j < data_x_size; ++j) {
                    evaluate1<true, 2>(p_xx, scaled_data_x, scaled_data_x, i, j);
                }
            }
            return p_xx;
        };

        container::ThreadLocalWrapper<WeightedDistribution2D> p_ax_all(form_factor::get_count(), constants::axes::d_axis.bins);
        auto calc_ax = [&data_a, &scaled_data_x, &p_ax_all, data_x_size] (int imin, int imax) {
            auto& p_ax = p_ax_all.get();
            for (int i = imin; i < imax; ++i) { // atoms
                int j = 0;                      // exv
                for (; j+7 < data_x_size; j+=8) {
                    grid::evaluate8<true, 1>(p_ax, data_a, scaled_data_x, i, j);
                }

                for (; j+3 < data_x_size; j+=4) {
                    grid::evaluate4<true, 1>(p_ax, data_a, scaled_data_x, i, j);
                }

                for (; j < data_x_size; ++j) {
                    grid::evaluate1<true, 1>(p_ax, data_a, scaled_data_x, i, j);
                }
            }
            return p_ax;
        };

        container::ThreadLocalWrapper<WeightedDistribution1D> p_wx_all(constants::axes::d_axis.bins);
        auto calc_wx = [&data_w, &scaled_data_x, &p_wx_all, data_x_size] (int imin, int imax) {
            auto& p_wx = p_wx_all.get();
            for (int i = imin; i < imax; ++i) { // waters
                int j = 0;                      // exv
                for (; j+7 < data_x_size; j+=8) {
                    evaluate8<true, 1>(p_wx, data_w, scaled_data_x, i, j);
                }

                for (; j+3 < data_x_size; j+=4) {
                    evaluate4<true, 1>(p_wx, data_w, scaled_data_x, i, j);
                }

                for (; j < data_x_size; ++j) {
                    evaluate1<true, 1>(p_wx, data_w, scaled_data_x, i, j);
                }
            }
            return p_wx;
        };

        //##############//
        // SUBMIT TASKS //
        //##############//
        int job_size = settings::general::detail::job_size;
        for (int i = 0; i < (int) data_x_size; i+=job_size) {
            pool->detach_task(
                [&calc_xx, i, job_size, data_x_size] () {return calc_xx(i, std::min(i+job_size, data_x_size));}
            );
        }
        for (int i = 0; i < (int) data_a_size; i+=job_size) {
            pool->detach_task(
                [&calc_ax, i, job_size, data_a_size] () {return calc_ax(i, std::min(i+job_size, data_a_size));}
            );
        }
        for (int i = 0; i < (int) data_w_size; i+=job_size) {
            pool->detach_task(
                [&calc_wx, i, job_size, data_w_size] () {return calc_wx(i, std::min(i+job_size, data_w_size));}
            );
        }

        pool->wait();
        WeightedDistribution1D p_xx_generic = p_xx_all.merge();
        WeightedDistribution2D p_ax_generic = p_ax_all.merge();
        WeightedDistribution1D p_wx_generic = p_wx_all.merge();

        p_xx_generic.add(0, data_x_size); // self-correlations

        // downsize our axes to only the relevant area
        unsigned int max_bin = 10; // minimum size is 10
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
        for (unsigned int i = 0; i < max_bin; ++i) {
            p_tot_ax.add_index(i, p_wx_generic.index(i));
        }

        for (unsigned int i = 0; i < p_ax_generic.size_x(); ++i) {
            for (unsigned int j = 0; j < max_bin; ++j) {
                p_tot_ax.add_index(j, p_ax_generic.index(i, j));
            }
        }

        // overwrite the excluded volume calculations from the HistogramManagerMTFFAvg calculations with our new grid-based ones
        // first cast the weighted distributions to make iteration simpler
        Distribution2D p_ax = Distribution2D(std::move(p_ax_generic));
        Distribution1D p_wx = Distribution1D(std::move(p_wx_generic));
        Distribution1D p_xx = Distribution1D(p_xx_generic);

        // replace the calculations
        for (unsigned int i = 0; i < p_aa.size_x(); ++i) {
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