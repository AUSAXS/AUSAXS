// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGridSurface.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <grid/exv/RawGridWithSurfaceExv.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <hist/distance_calculator/detail/TemplateHelperAvg.h>
#include <hist/distance_calculator/detail/TemplateHelperGrid.h>
#include <form_factor/FormFactorType.h>
#include <utility/MultiThreading.h>
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
    auto pool = utility::multi_threading::get_global_pool();

    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized
    hist::detail::CompactCoordinatesFF<variable_bin_width> data_x_i, data_x_s;

    {   // generate the excluded volume representation
        auto exv = get_exv();
        std::vector<data::AtomFF> interior(exv.interior.size()), surface(exv.surface.size());
        std::transform(
            exv.interior.begin(), exv.interior.end(), interior.begin(),
            [] (const Vector3<double>& atom) {return data::AtomFF{atom, form_factor::form_factor_t::EXCLUDED_VOLUME};}
        );
        std::transform(
            exv.surface.begin(), exv.surface.end(), surface.begin(),
            [] (const Vector3<double>& atom) {return data::AtomFF{atom, form_factor::form_factor_t::EXCLUDED_VOLUME};}
        );
        data_x_i = hist::detail::CompactCoordinatesFF<variable_bin_width>(std::move(interior));
        data_x_s = hist::detail::CompactCoordinatesFF<variable_bin_width>(std::move(surface));
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
    container::ThreadLocalWrapper<XXContainer> p_xx_all(settings::axes::bin_count);
    auto calc_xx_ii = [&data_x_i, &p_xx_all, data_x_i_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv interior
            int j = i+1;                    // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                evaluate8<variable_bin_width, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                evaluate4<variable_bin_width, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                evaluate1<variable_bin_width, 2>(p_xx.interior, data_x_i, data_x_i, i, j);
            }
        }
        return p_xx;
    };

    auto calc_xx_ss = [&data_x_s, &p_xx_all, data_x_s_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv surface
            int j = i+1;                    // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<variable_bin_width, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<variable_bin_width, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<variable_bin_width, 2>(p_xx.surface, data_x_s, data_x_s, i, j);
            }
        }
        return p_xx;
    };

    auto calc_xx_si = [&data_x_i, &data_x_s, &p_xx_all, data_x_s_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv interior
            int j = 0;                      // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<variable_bin_width, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<variable_bin_width, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<variable_bin_width, 2>(p_xx.cross, data_x_i, data_x_s, i, j);
            }
        }
        return p_xx;
    };

    container::ThreadLocalWrapper<AXContainer> p_ax_all(form_factor::get_count(), settings::axes::bin_count);
    auto calc_ax = [&data_a, &data_x_i, &data_x_s, &p_ax_all, data_x_i_size, data_x_s_size] (int imin, int imax) {
        auto& p_ax = p_ax_all.get();
        for (int i = imin; i < imax; ++i) { // atoms
            int j = 0;                      // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                detail::grid::evaluate8<variable_bin_width, false, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                detail::grid::evaluate4<variable_bin_width, false, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                detail::grid::evaluate1<variable_bin_width, false, 1>(p_ax.interior, data_a, data_x_i, i, j);
            }

            j = 0;                          // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                detail::grid::evaluate8<variable_bin_width, false, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                detail::grid::evaluate4<variable_bin_width, false, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                detail::grid::evaluate1<variable_bin_width, false, 1>(p_ax.surface, data_a, data_x_s, i, j);
            }
        }
        return p_ax;
    };

    container::ThreadLocalWrapper<WXContainer> p_wx_all(settings::axes::bin_count);
    auto calc_wx = [&data_w, &data_x_i, &data_x_s, &p_wx_all, data_x_i_size, data_x_s_size] (int imin, int imax) {
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // waters
            int j = 0;                      // exv interior
            for (; j+7 < data_x_i_size; j+=8) {
                evaluate8<variable_bin_width, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            for (; j+3 < data_x_i_size; j+=4) {
                evaluate4<variable_bin_width, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            for (; j < data_x_i_size; ++j) {
                evaluate1<variable_bin_width, 1>(p_wx.interior, data_w, data_x_i, i, j);
            }

            j = 0;                          // exv surface
            for (; j+7 < data_x_s_size; j+=8) {
                evaluate8<variable_bin_width, 1>(p_wx.surface, data_w, data_x_s, i, j);
            }

            for (; j+3 < data_x_s_size; j+=4) {
                evaluate4<variable_bin_width, 1>(p_wx.surface, data_w, data_x_s, i, j);
            }

            for (; j < data_x_s_size; ++j) {
                evaluate1<variable_bin_width, 1>(p_wx.surface, data_w, data_x_s, i, j);
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
    p_xx.interior.add_index(0, 0, data_x_i_size);
    p_xx.surface.add_index(0, 0, data_x_s_size);

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

    Distribution3D p_aa = std::move(cast_res->get_raw_aa_counts_by_ff());
    Distribution2D p_aw = std::move(cast_res->get_raw_aw_counts_by_ff());
    Distribution1D p_ww = std::move(cast_res->get_raw_ww_counts_by_ff());

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

template class ausaxs::hist::HistogramManagerMTFFGridSurface<true>;
template class ausaxs::hist::HistogramManagerMTFFGridSurface<false>;