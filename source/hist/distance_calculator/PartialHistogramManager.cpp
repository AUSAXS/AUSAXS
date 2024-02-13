#include "hist/distribution/detail/WeightedEntry.h"
#include <algorithm>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <data/state/StateManager.h>
#include <data/Molecule.h>
#include <settings/HistogramSettings.h>
#include <constants/Axes.h>

using namespace hist;

template<bool use_weighted_distribution> 
PartialHistogramManager<use_weighted_distribution>::PartialHistogramManager(observer_ptr<const data::Molecule> protein) 
    : HistogramManager<use_weighted_distribution>(protein), coords_a(this->body_size), partials_aa(this->body_size, this->body_size), partials_aw(this->body_size) {}

template<bool use_weighted_distribution> 
PartialHistogramManager<use_weighted_distribution>::~PartialHistogramManager() = default;

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> PartialHistogramManager<use_weighted_distribution>::calculate() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    const std::vector<bool> externally_modified = this->statemanager->get_externally_modified_bodies();
    const std::vector<bool> internally_modified = this->statemanager->get_internally_modified_bodies();

    // check if the object has already been initialized
    if (this->master.size() == 0) [[unlikely]] {
        initialize(); 
    } 
    
    // if not, we must first check if the coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < this->body_size; i++) {
            if (internally_modified[i]) {
                // if the internal state was modified, we have to recalculate the self-correlation
                calc_self_correlation(i);
            } else if (externally_modified[i]) {
                // if the external state was modified, we have to update the coordinate representations
                this->coords_a[i] = detail::CompactCoordinates(this->protein->get_body(i));
            }
        }
    }

    // check if the hydration layer was modified
    if (this->statemanager->get_modified_hydration()) {
        this->coords_w = detail::CompactCoordinates(this->protein->get_waters()); // if so, first update the compact coordinate representation
        calc_ww(); // then update the partial histogram

        // iterate through the lower triangle
        for (unsigned int i = 0; i < this->body_size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) {
                    calc_aa(i, j);
                }
            }
            calc_aw(i); // we then update its partial histograms
        }
    }

    // if the hydration layer was not modified
    else {
        for (unsigned int i = 0; i < this->body_size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) { // if either of the two bodies were modified
                    calc_aa(i, j); // recalculate their partial histogram
                }
            }
            if (externally_modified[i]) { // if a body was modified
                calc_aw(i); // update its partial histogram with the hydration layer
            }
        }
    }
    this->statemanager.reset();

    // downsize our axes to only the relevant area
    GenericDistribution1D_t p_tot = this->master;
    int max_bin = 10; // minimum size is 10
    for (int i = (int) p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_tot.resize(max_bin);
    return std::make_unique<DistanceHistogram>(std::move(p_tot));
}

template<bool use_weighted_distribution> 
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManager<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;

    auto total = calculate();
    int bins = total->get_total_counts().size();

    // determine p_tot
    GenericDistribution1D_t p_tot(bins);
    for (int i = 0; i < bins; i++) {
        p_tot.index(i) = master.index(i);
    }

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    Distribution1D p_ww = Distribution1D(partials_ww);
    Distribution1D p_aa = Distribution1D(this->master.base);
    Distribution1D p_aw(bins);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; i++) {
        for (unsigned int j = 0; j < i; j++) {
            GenericDistribution1D_t& current = partials_aa.index(i, j);

            // iterate through each entry in the partial histogram
            for (int k = 0; k < bins; k++) {
                p_aa.add(k, current.get_content(k)); // add to p_aa
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; i++) {
        GenericDistribution1D_t& current = partials_aw.index(i);

        // iterate through each entry in the partial histogram
        for (int k = 0; k < bins; k++) {
            p_aw.add(k, current.get_content(k)); // add to p_aa
        }
    }

    // p_aw is already resized
    p_ww.resize(bins);
    p_aa.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(
        std::move(p_aa), 
        std::move(p_aw), 
        std::move(p_ww), 
        std::move(p_tot)
    );
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_self_correlation(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    detail::CompactCoordinates current(this->protein->get_body(index));

    // calculate internal distances between atoms
    GenericDistribution1D_t p_aa(this->master.axis.bins);
    for (unsigned int i = 0; i < current.size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < current.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_aa, current, current, i, j);
        }

        for (; j+3 < current.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_aa, current, current, i, j);
        }

        for (; j < current.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_aa, current, current, i, j);
        }
    }

    // calculate self-correlation
    p_aa.add(0, std::accumulate(current.get_data().begin(), current.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));

    // store the coordinates for later
    this->coords_a[index] = std::move(current);

    this->master.base -= partials_aa.index(index, index);
    this->master -= partials_aa.index(index, index);
    partials_aa.index(index, index) = std::move(p_aa);
    this->master += partials_aa.index(index, index);
    this->master.base += partials_aa.index(index, index);
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_aa(unsigned int n, unsigned int m) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    detail::CompactCoordinates& coords_n = this->coords_a[n];
    detail::CompactCoordinates& coords_m = this->coords_a[m];

    GenericDistribution1D_t p_aa(this->master.axis.bins);
    for (unsigned int i = 0; i < coords_n.size(); i++) {
        unsigned int j = 0;
        for (; j+7 < coords_m.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_aa, coords_n, coords_m, i, j);
        }

        for (; j+3 < coords_m.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_aa, coords_n, coords_m, i, j);
        }

        for (; j < coords_m.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_aa, coords_n, coords_m, i, j);
        }
    }

    this->master -= partials_aa.index(n, m);
    partials_aa.index(n, m) = std::move(p_aa);
    this->master += partials_aa.index(n, m);
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::initialize() {
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<use_weighted_distribution>(p_base, axis);

    partials_ww = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
    for (unsigned int n = 0; n < this->body_size; n++) {
        partials_aw.index(n) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        partials_aa.index(n, n) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        calc_self_correlation(n);

        for (unsigned int k = 0; k < n; k++) {
            partials_aa.index(n, k) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        }
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_aw(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    detail::CompactCoordinates& coords = this->coords_a[index];

    GenericDistribution1D_t p_aw(this->master.axis.bins);
    for (unsigned int i = 0; i < coords.size(); i++) {
        unsigned int j = 0;
        for (; j+7 < this->coords_w.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_aw, coords, this->coords_w, i, j);
        }

        for (; j+3 < this->coords_w.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_aw, coords, this->coords_w, i, j);
        }

        for (; j < this->coords_w.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_aw, coords, this->coords_w, i, j);
        }
    }

    this->master -= partials_aw.index(index); // subtract the previous hydration histogram
    partials_aw.index(index) = std::move(p_aw);
    this->master += partials_aw.index(index); // add the new hydration histogram
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_ww() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    GenericDistribution1D_t p_ww(this->master.axis.bins);

    // calculate internal distances for the hydration layer
    for (unsigned int i = 0; i < this->coords_w.size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < this->coords_w.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_ww, this->coords_w, this->coords_w, i, j);
        }

        for (; j+3 < this->coords_w.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_ww, this->coords_w, this->coords_w, i, j);
        }

        for (; j < this->coords_w.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_ww, this->coords_w, this->coords_w, i, j);
        }
    }

    // calculate self-correlation
    p_ww.add(0, std::accumulate(this->coords_w.get_data().begin(), this->coords_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));

    this->master -= partials_ww; // subtract the previous hydration histogram
    partials_ww = std::move(p_ww);
    this->master += partials_ww; // add the new hydration histogram
}

template class hist::PartialHistogramManager<true>;
template class hist::PartialHistogramManager<false>;