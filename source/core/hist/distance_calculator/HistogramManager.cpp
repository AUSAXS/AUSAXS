/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/state/StateManager.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <settings/HistogramSettings.h>
#include <constants/Axes.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>

using namespace hist;

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::HistogramManager(observer_ptr<const data::Molecule> protein) : IHistogramManager(protein), protein(protein) {}

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::HistogramManager(const data::Molecule& protein) : HistogramManager(&protein) {}

template<bool use_weighted_distribution>
HistogramManager<use_weighted_distribution>::~HistogramManager() = default;

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> HistogramManager<use_weighted_distribution>::calculate() {return calculate_all();}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> HistogramManager<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    GenericDistribution1D_t p_aa(constants::axes::d_axis.bins);
    GenericDistribution1D_t p_ww(constants::axes::d_axis.bins);
    GenericDistribution1D_t p_aw(constants::axes::d_axis.bins);

    hist::detail::CompactCoordinates data_a(protein->get_bodies());
    hist::detail::CompactCoordinates data_w = hist::detail::CompactCoordinates(protein->get_waters());
    int data_a_size = (int) data_a.size();
    int data_w_size = (int) data_w.size();

    // calculate aa distances
    for (int i = 0; i < data_a_size; ++i) {
        int j = i+1;
        for (; j+7 < data_a_size; j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
        }

        for (; j+3 < data_a_size; j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
        }

        for (; j < data_a_size; ++j) {
            evaluate1<use_weighted_distribution, 2>(p_aa, data_a, data_a, i, j);
        }
    }

    for (int i = 0; i < data_w_size; ++i) {
        // calculate ww distances
        {
            int j = i+1;
            for (; j+7 < data_w_size; j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j+3 < data_w_size; j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
            }

            for (; j < data_w_size; ++j) {
                evaluate1<use_weighted_distribution, 2>(p_ww, data_w, data_w, i, j);
            }
        }
        
        // calculate aw distances
        {
            int j = 0;
            for (; j+7 < data_a_size; j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }

            for (; j+3 < data_a_size; j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }

            for (; j < data_a_size; ++j) {
                evaluate1<use_weighted_distribution, 1>(p_aw, data_w, data_a, i, j);
            }
        }
    }

    // add self-correlation
    p_aa.add(0, std::accumulate(data_a.get_data().begin(), data_a.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);}));
    p_ww.add(0, std::accumulate(data_w.get_data().begin(), data_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + std::pow(val.value.w, 2);}));

    // calculate p_tot
    GenericDistribution1D_t p_tot(constants::axes::d_axis.bins);
    for (int i = 0; i < (int) p_aa.size(); ++i) {p_tot.index(i) = p_aa.index(i) + p_ww.index(i) + 2*p_aw.index(i);}

    // downsize our axes to only the relevant area
    int max_bin = 10; // minimum size is 10
    for (int i = (int) p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }

    p_aa.resize(max_bin);
    p_ww.resize(max_bin);
    p_aw.resize(max_bin);
    p_tot.resize(max_bin);

    if constexpr (use_weighted_distribution) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(p_aa)), 
            std::move(Distribution1D(p_aw)), 
            std::move(Distribution1D(p_ww)), 
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