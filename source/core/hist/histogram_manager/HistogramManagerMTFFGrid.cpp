// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/distance_calculator/detail/TemplateHelperBase.h> // For ausaxs::detail::add8/4/1::evaluate
#include <hist/distance_calculator/detail/TemplateHelperAvg.h> // For ausaxs::evaluate8/4/1 with factor
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFAvg.h>
#include <hist/intensity_calculator/CompositeDistanceHistogramFFGrid.h>
#include <container/ThreadLocalWrapper.h>
#include <data/Molecule.h>
#include <grid/Grid.h>
#include <grid/exv/RawGridExv.h>
#include <settings/GeneralSettings.h>
#include <settings/GridSettings.h>
#include <settings/HistogramSettings.h>
#include <form_factor/FormFactorType.h>
#include <utility/MultiThreading.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::container;

// Atom-excluded volume (2D: ff, distance)
// Extract FF before calling evaluate, then use distances only (FF info comes from atom)
namespace ff2d {
template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate8(
    typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_a, 
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    int ff_i = data_a.get_ff_type(i);
    // Create temporary CompactCoordinates views to use existing evaluate infrastructure
    // CompactCoordinatesTemplate::get_data() returns the underlying vector
    // We evaluate using CompactCoordinatesXYZW (from data_a template) with CompactCoordinatesXYZW (from data_x)
    auto& data_a_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_a);
    auto res = ausaxs::detail::add8::evaluate<weighted_bins>(data_a_base, data_x, i, j);
    for (unsigned int k = 0; k < 8; ++k) {p.add(ff_i, res.distances[k], factor*res.weights[k]);}
}

template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate4(
    typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_a, 
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    int ff_i = data_a.get_ff_type(i);
    auto& data_a_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_a);
    auto res = ausaxs::detail::add4::evaluate<weighted_bins>(data_a_base, data_x, i, j);
    for (unsigned int k = 0; k < 4; ++k) {p.add(ff_i, res.distances[k], factor*res.weights[k]);}
}

template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate1(
    typename hist::GenericDistribution2D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_a, 
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    int ff_i = data_a.get_ff_type(i);
    auto& data_a_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_a);
    auto res = ausaxs::detail::add1::evaluate<weighted_bins>(data_a_base, data_x, i, j);
    p.add(ff_i, res.distance, factor*res.weight);
}
} // namespace ff2d

// Water-excluded volume (1D: distance) - water has constant FF, cast to base and ignore FF info
namespace ff1d {
template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate8(
    typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_w,
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    auto& data_w_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_w);
    auto res = ausaxs::detail::add8::evaluate<weighted_bins>(data_w_base, data_x, i, j);
    for (unsigned int k = 0; k < 8; ++k) {p.template add<factor>(res.distances[k], res.weights[k]);}
}

template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate4(
    typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_w,
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    auto& data_w_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_w);
    auto res = ausaxs::detail::add4::evaluate<weighted_bins>(data_w_base, data_x, i, j);
    for (unsigned int k = 0; k < 4; ++k) {p.template add<factor>(res.distances[k], res.weights[k]);}
}

template<bool weighted_bins, bool variable_bin_width, int factor>
void evaluate1(
    typename hist::GenericDistribution1D<weighted_bins>::type& p, const hist::detail::CompactCoordinatesFF<variable_bin_width>& data_w,
    const hist::detail::CompactCoordinates<variable_bin_width>& data_x, int i, int j
) {
    auto& data_w_base = reinterpret_cast<const hist::detail::CompactCoordinates<variable_bin_width>&>(data_w);
    auto res = ausaxs::detail::add1::evaluate<weighted_bins>(data_w_base, data_x, i, j);
    p.template add<factor>(res.distance, res.weight);
}
} // namespace ff1d

template<bool variable_bin_width>
HistogramManagerMTFFGrid<variable_bin_width>::~HistogramManagerMTFFGrid() = default;

template<bool variable_bin_width>
std::unique_ptr<DistanceHistogram> HistogramManagerMTFFGrid<variable_bin_width>::calculate() {
    return calculate_all();
}

template<bool variable_bin_width>
grid::exv::GridExcludedVolume HistogramManagerMTFFGrid<variable_bin_width>::get_exv() const {
    return grid::exv::RawGridExv::create(this->protein->get_grid());
}

template<bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManagerMTFFGrid<variable_bin_width>::calculate_all() {
    logging::log("HistogramManagerMTFFGrid::calculate: starting calculation");
    auto pool = utility::multi_threading::get_global_pool();

    auto base_res = HistogramManagerMTFFAvg<true, variable_bin_width>::calculate_all(); // make sure everything is initialized
    auto data_x = hist::detail::CompactCoordinates<variable_bin_width>(std::move(get_exv().interior), 1);
    auto& data_a = *this->data_a_ptr;
    auto& data_w = *this->data_w_ptr;
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();
    int data_x_size = (int) data_x.size();

    //########################//
    // PREPARE MULTITHREADING //
    //########################//
    container::ThreadLocalWrapper<WeightedDistribution1D> p_xx_all(settings::axes::bin_count);
    auto calc_xx = [&data_x, &p_xx_all, data_x_size] (int imin, int imax) {
        auto& p_xx = p_xx_all.get();
        for (int i = imin; i < imax; ++i) { // exv
            int j = i+1;                    // exv
            for (; j+7 < data_x_size; j+=8) {
                evaluate8<true, variable_bin_width, 2>(p_xx, data_x, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                evaluate4<true, variable_bin_width, 2>(p_xx, data_x, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                evaluate1<true, variable_bin_width, 2>(p_xx, data_x, data_x, i, j);
            }
        }
        return p_xx;
    };

    container::ThreadLocalWrapper<WeightedDistribution2D> p_ax_all(form_factor::get_count(), settings::axes::bin_count);
    auto calc_ax = [&data_a, &data_x, &p_ax_all, data_x_size] (int imin, int imax) {
        auto& p_ax = p_ax_all.get();
        for (int i = imin; i < imax; ++i) { // atoms
            int j = 0;                      // exv
            for (; j+7 < data_x_size; j+=8) {
                ff2d::evaluate8<true, variable_bin_width, 1>(p_ax, data_a, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                ff2d::evaluate4<true, variable_bin_width, 1>(p_ax, data_a, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                ff2d::evaluate1<true, variable_bin_width, 1>(p_ax, data_a, data_x, i, j);
            }
        }
        return p_ax;
    };

    container::ThreadLocalWrapper<WeightedDistribution1D> p_wx_all(settings::axes::bin_count);
    auto calc_wx = [&data_w, &data_x, &p_wx_all, data_x_size] (int imin, int imax) {
        auto& p_wx = p_wx_all.get();
        for (int i = imin; i < imax; ++i) { // waters
            int j = 0;                      // exv
            for (; j+7 < data_x_size; j+=8) {
                ff1d::evaluate8<true, variable_bin_width, 1>(p_wx, data_w, data_x, i, j);
            }

            for (; j+3 < data_x_size; j+=4) {
                ff1d::evaluate4<true, variable_bin_width, 1>(p_wx, data_w, data_x, i, j);
            }

            for (; j < data_x_size; ++j) {
                ff1d::evaluate1<true, variable_bin_width, 1>(p_wx, data_w, data_x, i, j);
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

    //###################//
    // SELF-CORRELATIONS //
    //###################//
    p_xx_generic.add(0, data_x_size);

    // downsize our axes to only the relevant area
    unsigned int max_bin = 10; // minimum size is 10
    for (int i = p_xx_generic.size()-1; i >= 10; i--) {
        if (p_xx_generic.index(i) != 0 || p_wx_generic.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    // ensure that our new vectors are compatible with those from the base class
    // also note that the order matters here, since we move data away from the cast_res object. Thus p_tot *must* be moved first. 
    auto cast_res = static_cast<CompositeDistanceHistogramFFAvg*>(base_res.get());
    WeightedDistribution1D p_tot = std::move(cast_res->get_counts());
    p_tot.set_bin_centers(cast_res->get_d_axis());

    Distribution3D p_aa = std::move(cast_res->get_aa_counts_by_ff());
    Distribution2D p_aw = std::move(cast_res->get_aw_counts_by_ff());
    Distribution1D p_ww = std::move(cast_res->get_ww_counts_by_ff());

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