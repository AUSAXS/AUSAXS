#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <data/state/StateManager.h>
#include <settings/HistogramSettings.h>
#include <constants/Axes.h>

using namespace hist;

template<bool use_weighted_distribution> 
PartialHistogramManager<use_weighted_distribution>::PartialHistogramManager(observer_ptr<const data::Molecule> protein) 
    : HistogramManager<use_weighted_distribution>(protein), coords_p(this->body_size), partials_pp(this->body_size, this->body_size), partials_hp(this->body_size) {}

template<bool use_weighted_distribution> 
PartialHistogramManager<use_weighted_distribution>::~PartialHistogramManager() = default;

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> PartialHistogramManager<use_weighted_distribution>::calculate() {
    const std::vector<bool> externally_modified = this->statemanager->get_externally_modified_bodies();
    const std::vector<bool> internally_modified = this->statemanager->get_internally_modified_bodies();

    // check if the object has already been initialized
    if (master.get_counts().size() == 0) [[unlikely]] {
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
                coords_p[i] = detail::CompactCoordinates(this->protein->get_body(i));
            }
        }
    }

    // check if the hydration layer was modified
    if (this->statemanager->get_modified_hydration()) {
        coords_h = detail::CompactCoordinates(this->protein->get_waters()); // if so, first update the compact coordinate representation
        calc_hh(); // then update the partial histogram

        // iterate through the lower triangle
        for (unsigned int i = 0; i < this->body_size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) {
                    calc_pp(i, j);
                }
            }
            calc_hp(i); // we then update its partial histograms
        }
    }

    // if the hydration layer was not modified
    else {
        for (unsigned int i = 0; i < this->body_size; i++) {
            for (unsigned int j = 0; j < i; j++) {
                if (externally_modified[i] || externally_modified[j]) { // if either of the two bodies were modified
                    calc_pp(i, j); // recalculate their partial histogram
                }
            }
            if (externally_modified[i]) { // if a body was modified
                calc_hp(i); // update its partial histogram with the hydration layer
            }
        }
    }
    this->statemanager.reset();
    Distribution1D p = master.get_counts();
    return std::make_unique<DistanceHistogram>(std::move(p), master.get_axis());
}

template<bool use_weighted_distribution> 
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManager<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;

    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_hh = partials_hh.get_counts();
    GenericDistribution1D_t p_pp = master.base.get_counts();
    GenericDistribution1D_t p_hp(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; i++) {
        for (unsigned int j = 0; j <= i; j++) {
            detail::PartialHistogram& current = partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < bins; k++) {
                p_pp.index(k) += current.get_count(k); // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; i++) {
        detail::PartialHistogram& current = partials_hp.index(i);

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < bins; k++) {
            p_hp.index(k) += current.get_count(k); // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(bins);
    p_pp.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(
        std::move(p_pp), 
        std::move(p_hp), 
        std::move(p_hh), 
        std::move(total->get_counts()), 
        total->get_axis()
    );
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_self_correlation(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    detail::CompactCoordinates current(this->protein->get_body(index));

    // calculate internal distances between atoms
    GenericDistribution1D_t p_pp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < current.size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < current.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_pp, current, current, i, j);
        }

        for (; j+3 < current.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_pp, current, current, i, j);
        }

        for (; j < current.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_pp, current, current, i, j);
        }
    }

    // calculate self-correlation
    p_pp.add(0, std::accumulate(current.get_data().begin(), current.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));

    // store the coordinates for later
    coords_p[index] = std::move(current);

    master.base -= partials_pp.index(index, index);
    master -= partials_pp.index(index, index);
    partials_pp.index(index, index).get_counts() = std::move(p_pp.get_data());
    master += partials_pp.index(index, index);
    master.base += partials_pp.index(index, index);
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::initialize() {
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    master = detail::MasterHistogram(p_base, axis);

    partials_hh = detail::PartialHistogram(axis);
    for (unsigned int n = 0; n < this->body_size; n++) {
        partials_hp.index(n) = detail::PartialHistogram(axis);
        partials_pp.index(n, n) = detail::PartialHistogram(axis);
        calc_self_correlation(n);

        for (unsigned int k = 0; k < n; k++) {
            partials_pp.index(n, k) = detail::PartialHistogram(axis);
        }
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_hp(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    detail::CompactCoordinates& coords = coords_p[index];

    GenericDistribution1D_t p_hp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < coords.size(); i++) {
        unsigned int j = 0;
        for (; j+7 < coords_h.size(); j+=8) {
            evaluate8<use_weighted_distribution, 1>(p_hp, coords, coords_h, i, j);
        }

        for (; j+3 < coords_h.size(); j+=4) {
            evaluate4<use_weighted_distribution, 1>(p_hp, coords, coords_h, i, j);
        }

        for (; j < coords_h.size(); ++j) {
            evaluate1<use_weighted_distribution, 1>(p_hp, coords, coords_h, i, j);
        }
    }

    master -= partials_hp.index(index)*2; // subtract the previous hydration histogram
    partials_hp.index(index).get_counts() = std::move(p_hp.get_data());
    master += partials_hp.index(index)*2; // add the new hydration histogram
}

template<bool use_weighted_distribution> 
void PartialHistogramManager<use_weighted_distribution>::calc_hh() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    GenericDistribution1D_t p_hh(master.get_axis().bins, 0);

    // calculate internal distances for the hydration layer
    for (unsigned int i = 0; i < coords_h.size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < coords_h.size(); j+=8) {
            evaluate8<use_weighted_distribution, 2>(p_hh, coords_h, coords_h, i, j);
        }

        for (; j+3 < coords_h.size(); j+=4) {
            evaluate4<use_weighted_distribution, 2>(p_hh, coords_h, coords_h, i, j);
        }

        for (; j < coords_h.size(); ++j) {
            evaluate1<use_weighted_distribution, 2>(p_hh, coords_h, coords_h, i, j);
        }
    }

    // calculate self-correlation
    p_hh.add(0, std::accumulate(coords_h.get_data().begin(), coords_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.get_counts() = std::move(p_hh.get_data());
    master += partials_hh; // add the new hydration histogram
}

template class hist::PartialHistogramManager<true>;
template class hist::PartialHistogramManager<false>;