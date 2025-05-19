/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
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

HistogramManagerMTFFGridSurface::~HistogramManagerMTFFGridSurface() = default;

std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGridSurface::calculate() {
    return calculate_all();
}

std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGridSurface::calculate_all() {
    logging::log("HistogramManagerMTFFGridSurface::calculate: starting calculation");
    using XXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::XXContainer;
    using AXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::AXContainer;
    using WXContainer = typename hist::CompositeDistanceHistogramFFGridSurface::WXContainer;
    auto pool = utility::multi_threading::get_global_pool();

    auto base_res = HistogramManagerMTFFAvg<true>::calculate_all(); // make sure everything is initialized
    hist::detail::CompactCoordinates data_x_i, data_x_s;

    {   // generate the excluded volume representation
        auto exv = this->protein->get_grid()->generate_excluded_volume(true);
        data_x_i = hist::detail::CompactCoordinates(std::move(exv.interior), 1);
        data_x_s = hist::detail::CompactCoordinates(std::move(exv.surface), 1);
    }

    auto& data_a = *this->data_a_ptr;
    auto& data_w = *this->data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();
    int data_x_i_size = (int) data_x_i.size();
    int data_x_s_size = (int) data_x_s.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<XXContainer> p_xx_all(constants::axes::d_axis.bins);
    auto calc_xx_ii = [&data_x_i, &p_xx_all, data_x_i_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv interior
            int j = i+1;                    // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                evaluate8<true, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                evaluate4<true, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                evaluate1<true, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }
        }
        return p_xx;
    };

    auto calc_xx_ss = [&data_x_s, &p_xx_all, data_x_s_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv surface
            int j = i+1;                    // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<true, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<true, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<true, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }
        }
        return p_xx;
    };

    auto calc_xx_si = [&data_x_i, &data_x_s, &p_xx_all, data_x_s_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv interior
            int j = 0;                      // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<true, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<true, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<true, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }
        }
        return p_xx;
    };

    container::ThreadLocalWrapper<AXContainer> p_ax_all(form_factor::get_count(), constants::axes::d_axis.bins);
    auto calc_ax = [&data_a, &data_x_i, &data_x_s, &p_ax_all, data_x_i_size, data_x_s_size] (int imin, int imax) {
        auto& p_ax = p_ax_all.get();
        for (int i = imin; i < imax; ++i) { // atoms
            int j = 0;                      // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                grid::evaluate8<true, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                grid::evaluate4<true, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                grid::evaluate1<true, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            j = 0;                          // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                grid::evaluate8<true, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                grid::evaluate4<true, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                grid::evaluate1<true, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }
        }
        return p_ax;
    };

    container::ThreadLocalWrapper<WXContainer> p_wx_all(constants::axes::d_axis.bins);
    auto calc_wx = [&data_w, &data_x_i, &data_x_s, &p_wx_all, data_x_i_size, data_x_s_size] (int imin, int imax) {
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // waters
            int j = 0;                      // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                evaluate8<true, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                evaluate4<true, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                evaluate1<true, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            j = 0;                          // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<true, 1>(p_wx.surface, data_w, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<true, 1>(p_wx.surface, data_w, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<true, 1>(p_wx.surface, data_w, data_x_s, i, j);
            }
        }
        return p_wx;
    };

    //##############//
    // SUBMIT TASKS //
    //##############//
    int job_size = settings::general::detail::job_size;
    for (int i = 0; i < (int) data_x_i_size; i+=job_size) {
        pool->detach_task(
            [&calc_xx_ii, i, job_size, data_x_i_size] () {return calc_xx_ii(i, std::min(i+job_size, data_x_i_size));}
        );
    }

    for (int i = 0; i < (int) data_x_s_size; i+=job_size) {
        pool->detach_task(
            [&calc_xx_ss, i, job_size, data_x_s_size] () {return calc_xx_ss(i, std::min(i+job_size, data_x_s_size));}
        );
    }

    for (int i = 0; i < (int) data_x_i_size; i+=job_size) {
        pool->detach_task(
            [&calc_xx_si, i, job_size, data_x_i_size] () {return calc_xx_si(i, std::min(i+job_size, data_x_i_size));}
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
    XXContainer p_xx = p_xx_all.merge();
    AXContainer p_ax = p_ax_all.merge();
    WXContainer p_wx = p_wx_all.merge();

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_xx.interior.add(0, data_x_i_size);
    p_xx.surface.add(0, data_x_s_size);

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (unsigned int i = p_xx.surface.size()-1; i >= 10; i--) {
        if (p_xx.surface.index(i) != 0 || p_xx.interior.index(i) != 0 || p_wx.surface.index(i) != 0 || p_wx.interior.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    WeightedDistribution1D p_tot = std::move(cast_res->get_counts());
    p_tot.set_bin_centers(cast_res->get_d_axis());

    Distribution3D p_aa = std::move(cast_res->get_aa_counts_ff());
    Distribution2D p_aw = std::move(cast_res->get_aw_counts_ff());
    Distribution1D p_ww = std::move(cast_res->get_ww_counts_ff());

    // either xx or ww are largest of all components
    max_bin = std::max<unsigned int>(max_bin, p_tot.size());

    // downsize the axes to only the relevant area
    if (base_res->get_d_axis().size() < max_bin) {
        p_aa.resize(max_bin);
        p_aw.resize(max_bin);
        p_ww.resize(max_bin);
    } else {
        max_bin = base_res->get_d_axis().size(); // make sure we overwrite anything which may already be stored
    }

    // calculate weighted distance bins
    p_tot.resize(max_bin);
    WeightedDistribution1D p_tot_ax = std::max<int>(max_bin, p_wx.surface.size());
    for (unsigned int i = 0; i < max_bin; ++i) {
        p_tot_ax.add_index(i, p_wx.interior.index(i));
        p_tot_ax.add_index(i, p_wx.surface.index(i));
    }

    for (unsigned int i = 0; i < p_ax.surface.size_x(); ++i) {
        for (unsigned int j = 0; j < max_bin; ++j) {
            p_tot_ax.add_index(j, p_ax.interior.index(i, j));
            p_tot_ax.add_index(j, p_ax.surface.index(i, j));
        }
    }

    WeightedDistribution1D p_tot_xx = std::max<int>(max_bin, p_xx.surface.size());
    for (unsigned int i = 0; i < max_bin; ++i) {
        p_tot_xx.add_index(i, p_xx.interior.index(i));
        p_tot_xx.add_index(i, p_xx.surface.index(i));
        p_tot_xx.add_index(i, p_xx.cross.index(i));
    }

    {   // delete the exv information from the HistogramManagerMTFFAvg data
        // we delegate this work to the DistanceHistogram class, since it must be able to do this anyway to vary the surface contribution
        for (unsigned int i = 0; i < p_aa.size_x(); ++i) {
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