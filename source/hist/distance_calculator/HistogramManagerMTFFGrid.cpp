/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "hist/distribution/Distribution1D.h"
#include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <constants/Axes.h>
#include <hist/distance_calculator/detail/TemplateHelpersFFAvg.h>
#include <form_factor/FormFactorType.h>
#include <utility/MultiThreading.h>

using namespace hist;

// custom evaluates for the grid since we don't want to account for the excluded volume
namespace grid {
    template<bool use_weighted_distribution, int factor>
    inline void evaluate8(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ::detail::add8::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 8; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool use_weighted_distribution, int factor>
    inline void evaluate4(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ::detail::add4::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        for (unsigned int k = 0; k < 4; ++k) {p.add(data_i.get_ff_type(i), res.distances[k], factor*res.weights[k]);}
    }

    template<bool use_weighted_distribution, int factor>
    inline void evaluate1(typename hist::GenericDistribution2D<use_weighted_distribution>::type& p, const hist::detail::CompactCoordinatesFF& data_i, const hist::detail::CompactCoordinates& data_j, int i, int j) {
        auto res = ::detail::add1::evaluate<use_weighted_distribution>(data_i, data_j, i, j);
        p.add(data_i.get_ff_type(i), res.distance, factor*res.weight);
    }
}

template<bool use_weighted_distribution> 
HistogramManagerMTFFGrid<use_weighted_distribution>::~HistogramManagerMTFFGrid() = default;

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGrid<use_weighted_distribution>::calculate() {
    return calculate_all();
}

template<bool use_weighted_distribution> 
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGrid<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    using GenericDistribution2D_t = typename hist::GenericDistribution2D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();

    auto base_res = HistogramManagerMTFFAvg<use_weighted_distribution>::calculate_all(); // make sure everything is initialized
    hist::detail::CompactCoordinates data_x = hist::detail::CompactCoordinates(this->protein->get_grid()->generate_excluded_volume(false).interior, 1);
    auto& data_a = *this->data_a_ptr;
    auto& data_w = *this->data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();
    int data_x_size = (int) data_x.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<GenericDistribution1D_t> p_xx_all(constants::axes::d_axis.bins);
    auto calc_xx = [&data_x, &p_xx_all, data_x_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv
            int j = i+1;                    // exv
            for (; j+7 < data_x_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_xx, data_x, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_xx, data_x, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_xx, data_x, data_x, i, j);
            }
        }
        return p_xx;
    };

    container::ThreadLocalWrapper<GenericDistribution2D_t> p_ax_all(form_factor::get_count(), constants::axes::d_axis.bins);
    auto calc_ax = [&data_a, &data_x, &p_ax_all, data_x_size] (int imin, int imax) {
        auto& p_ax = p_ax_all.get();
        for (int i = imin; i < imax; ++i) { // atoms
            int j = 0;                      // exv
            for (; j+7 < data_x_size; j+=8) {
                grid::evaluate8<use_weighted_distribution, 1>(p_ax, data_a, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                grid::evaluate4<use_weighted_distribution, 1>(p_ax, data_a, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                grid::evaluate1<use_weighted_distribution, 1>(p_ax, data_a, data_x, i, j);
            }
        }
        return p_ax;
    };

    container::ThreadLocalWrapper<GenericDistribution1D_t> p_wx_all(constants::axes::d_axis.bins);
    auto calc_wx = [&data_w, &data_x, &p_wx_all, data_x_size] (int imin, int imax) {
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // waters
            int j = 0;                      // exv
            for (; j+7 < data_x_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_wx, data_w, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_wx, data_w, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_wx, data_w, data_x, i, j);
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
    GenericDistribution1D_t p_xx_generic = p_xx_all.merge();
    GenericDistribution2D_t p_ax_generic = p_ax_all.merge();
    GenericDistribution1D_t p_wx_generic = p_wx_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_xx_generic.add(0, data_x_size);

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_xx_generic.size()-1; i >= 10; i--) {
        if (p_xx_generic.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    GenericDistribution1D_t p_tot = hist::Distribution1D(std::move(cast_res->get_counts()));
    if constexpr (use_weighted_distribution) {
        p_tot.set_bin_centers(cast_res->get_d_axis());
    }

    Distribution3D p_aa = std::move(cast_res->get_aa_counts_ff());
    Distribution2D p_aw = std::move(cast_res->get_aw_counts_ff());
    Distribution1D p_ww = std::move(cast_res->get_ww_counts_ff());

    // either xx or ww are largest of all components
    max_bin = std::max<unsigned int>(max_bin, p_tot.size());

    // update p_tot
    p_tot.resize(max_bin);
    GenericDistribution1D_t p_tot_ax = hist::Distribution1D(std::max<int>(max_bin, p_wx_generic.size()));
    for (int i = 0; i < (int) max_bin; ++i) {
        p_tot_ax.add_index(i, p_wx_generic.index(i));
    }

    for (int i = 0; i < (int) p_ax_generic.size_x(); ++i) {
        for (int j = 0; j < (int) max_bin; ++j) {
            p_tot_ax.add_index(j, p_ax_generic.index(i, j));
        }
    }

    // downsize the axes to only the relevant area
    if (base_res->get_d_axis().size() < max_bin) {
        p_aa.resize(max_bin);
        p_aw.resize(max_bin);
        p_ww.resize(max_bin);
    } else {
        max_bin = base_res->get_d_axis().size(); // make sure we overwrite anything which may already be stored
    }

    // redefine the distributions without weights
    Distribution2D p_ax = Distribution2D(p_ax_generic);
    Distribution1D p_wx = Distribution1D(p_wx_generic);
    Distribution1D p_xx = Distribution1D(p_xx_generic);

    for (unsigned int i = 0; i < p_aa.size_x(); ++i) {
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

template class hist::HistogramManagerMTFFGrid<false>;
template class hist::HistogramManagerMTFFGrid<true>;