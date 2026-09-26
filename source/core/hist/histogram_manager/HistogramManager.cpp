// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/histogram_manager/HistogramManager.h>

#include <data/Molecule.h>
#include <hist/detail/AtomOrdering.h>
#include <hist/detail/BinEstimate.h>
#include <hist/detail/CompactCoordinatesFactory.h>
#include <hist/detail/SimpleExvModel.h>
#include <hist/detail/data/BinWidth.h>
#include <hist/distance_calculator/detail/Evaluators.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool weighted_bins>
HistogramManager<weighted_bins>::HistogramManager(observer_ptr<const data::Molecule> protein) : protein(protein) {
    logging::log("initializing HistogramManager");
}

template<bool weighted_bins>
HistogramManager<weighted_bins>::~HistogramManager() = default;

template<bool weighted_bins>
std::unique_ptr<DistanceHistogram> HistogramManager<weighted_bins>::calculate() {return calculate_all();}

template<bool weighted_bins>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManager<weighted_bins>::calculate_all() {
    logging::log("HistogramManager::calculate: starting calculation");

    using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;

    auto data_a = hist::detail::factory::construct_from_atoms(protein);
    auto data_w = hist::detail::factory::construct_from_waters(protein);
    int data_a_size = data_a.size();
    int data_w_size = data_w.size();
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, protein);
    int bin_count = hist::detail::required_bin_count(data_a, data_w);
    hist::detail::decorrelate_order<weighted_bins>(bin_count, data_a, data_w);

    GenericDistribution1D_t p_aa(bin_count);
    GenericDistribution1D_t p_ww(bin_count);
    GenericDistribution1D_t p_aw(bin_count);
    auto b_aa = hist::detail::bins(p_aa);
    auto b_ww = hist::detail::bins(p_ww);
    auto b_aw = hist::detail::bins(p_aw);
    float inv_width = hist::detail::inv_bin_width();

    // calculate aa distances
    for (int i = 0; i < data_a_size; ++i) {
        int j = i+1;
        for (; j+15 < data_a_size; j+=16) {
            evaluate16<2>(b_aa, data_a, data_a, i, j, inv_width);
        }

        for (; j+7 < data_a_size; j+=8) {
            evaluate8<2>(b_aa, data_a, data_a, i, j, inv_width);
        }

        for (; j+3 < data_a_size; j+=4) {
            evaluate4<2>(b_aa, data_a, data_a, i, j, inv_width);
        }

        for (; j < data_a_size; ++j) {
            evaluate1<2>(b_aa, data_a, data_a, i, j, inv_width);
        }
    }

    for (int i = 0; i < data_w_size; ++i) {
        {   // calculate ww distances
            int j = i+1;
            for (; j+15 < data_w_size; j+=16) {
                evaluate16<2>(b_ww, data_w, data_w, i, j, inv_width);
            }

            for (; j+7 < data_w_size; j+=8) {
                evaluate8<2>(b_ww, data_w, data_w, i, j, inv_width);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<2>(b_ww, data_w, data_w, i, j, inv_width);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<2>(b_ww, data_w, data_w, i, j, inv_width);
            }
        }

        {   // calculate aw distances
            int j = 0;
            for (; j+15 < data_a_size; j+=16) {
                evaluate16<2>(b_aw, data_w, data_a, i, j, inv_width);
            }

            for (; j+7 < data_a_size; j+=8) {
                evaluate8<2>(b_aw, data_w, data_a, i, j, inv_width);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<2>(b_aw, data_w, data_a, i, j, inv_width);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<2>(b_aw, data_w, data_a, i, j, inv_width);
            }
        }
    }

    // add self-correlation
    auto sum_squared_weights = [] (const auto& set) {
        double sum = 0;
        for (int i = 0; i < set.size(); ++i) {sum += std::pow(set.get_weight(i), 2);}
        return sum;
    };
    double total_weight_aa = sum_squared_weights(data_a);
    double total_weight_ww = sum_squared_weights(data_w);
    if constexpr (weighted_bins) {
        p_aa.add_index(0, WeightedEntry(total_weight_aa, static_cast<int>(total_weight_aa), 0));
        p_ww.add_index(0, WeightedEntry(total_weight_ww, static_cast<int>(total_weight_ww), 0));
    } else {
        p_aa.add_index(0, total_weight_aa);
        p_ww.add_index(0, total_weight_ww);
    }

    // calculate p_tot
    GenericDistribution1D_t p_tot(bin_count);
    for (int i = 0; i < (int) p_aa.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + p_aw.index(i);}

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
    for (int i = (int) p_tot.size()-1; i >= 10; --i) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (weighted_bins) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(std::move(p_aa))), 
            std::move(Distribution1D(std::move(p_aw))), 
            std::move(Distribution1D(std::move(p_ww))), 
            std::move(p_tot)
        );
    } else {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(p_aa), 
            std::move(p_aw), 
            std::move(p_ww), 
            std::move(p_tot)
        );
    }
}

template class hist::HistogramManager<false>;
template class hist::HistogramManager<true>;