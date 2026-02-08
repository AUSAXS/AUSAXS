// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <hist/distance_calculator/detail/TemplateHelperSimple.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/SimpleExvModel.h>
#include <settings/HistogramSettings.h>
#include <constants/ConstantsAxes.h>
#include <utility/Logging.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool weighted_bins, bool variable_bin_width>
HistogramManager<weighted_bins, variable_bin_width>::HistogramManager(observer_ptr<const data::Molecule> protein) : protein(protein) {
    logging::log("initializing HistogramManager");
}

template<bool weighted_bins, bool variable_bin_width>
HistogramManager<weighted_bins, variable_bin_width>::~HistogramManager() = default;

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<DistanceHistogram> HistogramManager<weighted_bins, variable_bin_width>::calculate() {return calculate_all();}

template<bool weighted_bins, bool variable_bin_width>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManager<weighted_bins, variable_bin_width>::calculate_all() {
    logging::log("HistogramManager::calculate: starting calculation");

    using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
    GenericDistribution1D_t p_aa(settings::axes::bin_count);
    GenericDistribution1D_t p_ww(settings::axes::bin_count);
    GenericDistribution1D_t p_aw(settings::axes::bin_count);

    hist::detail::CompactCoordinates<variable_bin_width> data_a(protein->get_bodies());
    hist::detail::CompactCoordinates<variable_bin_width> data_w(protein->get_waters());
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(data_a, protein);

    // calculate aa distances
    for (int i = 0; i < data_a_size; ++i) {
        int j = i+1;
        for (; j+7 < data_a_size; j+=8) {
            evaluate8<variable_bin_width, 2>(p_aa, data_a, data_a, i, j);
        }

        for (; j+3 < data_a_size; j+=4) {
            evaluate4<variable_bin_width, 2>(p_aa, data_a, data_a, i, j);
        }

        for (; j < data_a_size; ++j) {
            evaluate1<variable_bin_width, 2>(p_aa, data_a, data_a, i, j);
        }
    }

    for (int i = 0; i < data_w_size; ++i) {
        {   // calculate ww distances
            int j = i+1;
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<variable_bin_width, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<variable_bin_width, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<variable_bin_width, 2>(p_ww, data_w, data_w, i, j);
            }
        }

        {   // calculate aw distances
            int j = 0;
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<variable_bin_width, 2>(p_aw, data_w, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<variable_bin_width, 2>(p_aw, data_w, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<variable_bin_width, 2>(p_aw, data_w, data_a, i, j);
            }
        }
    }

    // add self-correlation
    double total_weight_aa = std::accumulate(data_a.get_data().begin(), data_a.get_data().end(), 0.0, [](double sum, const auto& val) {return sum + std::pow(val.value.w, 2);});
    double total_weight_ww = std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const auto& val) {return sum + std::pow(val.value.w, 2);});
    if constexpr (weighted_bins) {
        p_aa.add_index(0, WeightedEntry(total_weight_aa, total_weight_aa, 0));
        p_ww.add_index(0, WeightedEntry(total_weight_ww, total_weight_ww, 0));
    } else {
        p_aa.add_index(0, total_weight_aa);
        p_ww.add_index(0, total_weight_ww);
    }

    // calculate p_tot
    GenericDistribution1D_t p_tot(settings::axes::bin_count);
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

template class hist::HistogramManager<false, false>;
template class hist::HistogramManager<false, true>;
template class hist::HistogramManager<true, false>;
template class hist::HistogramManager<true, true>;