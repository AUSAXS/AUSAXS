#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/state/StateManager.h>
#include <settings/HistogramSettings.h>
#include <constants/Constants.h>

using namespace hist;

PartialHistogramManager::PartialHistogramManager(const data::Molecule* const protein) : HistogramManager(protein), coords_p(body_size), partials_pp(body_size, body_size), partials_hp(body_size) {}

PartialHistogramManager::~PartialHistogramManager() = default;

std::unique_ptr<DistanceHistogram>  PartialHistogramManager::calculate() {
    const std::vector<bool> externally_modified = statemanager->get_externally_modified_bodies();
    const std::vector<bool> internally_modified = statemanager->get_internally_modified_bodies();

    // check if the object has already been initialized
    if (master.get_counts().size() == 0) [[unlikely]] {
        initialize(); 
    } 
    
    // if not, we must first check if the coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < body_size; i++) {
            if (internally_modified[i]) {
                // if the internal state was modified, we have to recalculate the self-correlation
                calc_self_correlation(i);
            } else if (externally_modified[i]) {
                // if the external state was modified, we have to update the coordinate representations
                coords_p[i] = detail::CompactCoordinates(protein->get_body(i));
            }
        }
    }

    // check if the hydration layer was modified
    if (statemanager->get_modified_hydration()) {
        coords_h = detail::CompactCoordinates(protein->get_waters()); // if so, first update the compact coordinate representation
        calc_hh(); // then update the partial histogram

        // iterate through the lower triangle
        for (unsigned int i = 0; i < body_size; i++) {
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
        for (unsigned int i = 0; i < body_size; i++) {
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
    statemanager.reset();
    std::vector<double> p = master.get_counts();
    return std::make_unique<DistanceHistogram>(std::move(p), master.get_axis());
}

std::unique_ptr<CompositeDistanceHistogram> PartialHistogramManager::calculate_all() {
    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    std::vector<double> p_hh = partials_hh.get_counts();
    std::vector<double> p_pp = master.base.get_counts();
    std::vector<double> p_hp(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < body_size; i++) {
        for (unsigned int j = 0; j <= i; j++) {
            detail::PartialHistogram& current = partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < bins; k++) {
                p_pp[k] += current.get_count(k); // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < body_size; i++) {
        detail::PartialHistogram& current = partials_hp.index(i);

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < bins; k++) {
            p_hp[k] += current.get_count(k); // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(bins);
    p_pp.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(total->get_counts()), total->get_axis());
}

void PartialHistogramManager::calc_self_correlation(unsigned int index) {
    detail::CompactCoordinates current(protein->get_body(index));

    // calculate internal distances between atoms
    std::vector<double> p_pp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < current.get_size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < current.get_size(); j+=8) {
            auto res = current[i].evaluate(current[j], current[j+1], current[j+2], current[j+3], current[j+4], current[j+5], current[j+6], current[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < current.get_size(); j+=4) {
            auto res = current[i].evaluate(current[j], current[j+1], current[j+2], current[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < current.get_size(); ++j) {
            auto res = current[i].evaluate(current[j]);
            p_pp[res.distance] += 2*res.weight;
        }
    }

    // calculate self-correlation
    p_pp[0] = std::accumulate(current.get_data().begin(), current.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );

    // store the coordinates for later
    coords_p[index] = std::move(current);

    master.base -= partials_pp.index(index, index);
    master -= partials_pp.index(index, index);
    partials_pp.index(index, index).get_counts() = std::move(p_pp);
    master += partials_pp.index(index, index);
    master.base += partials_pp.index(index, index);
}

/**
 * @brief This initializes some necessary variables and precalculates the internal distances between atoms in each body.
 */
void PartialHistogramManager::initialize() {
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    master = detail::MasterHistogram(p_base, axis);

    partials_hh = detail::PartialHistogram(axis);
    for (unsigned int n = 0; n < body_size; n++) {
        partials_hp.index(n) = detail::PartialHistogram(axis);
        partials_pp.index(n, n) = detail::PartialHistogram(axis);
        calc_self_correlation(n);

        for (unsigned int k = 0; k < n; k++) {
            partials_pp.index(n, k) = detail::PartialHistogram(axis);
        }
    }
}

void PartialHistogramManager::calc_pp(unsigned int n, unsigned int m) {
    detail::CompactCoordinates& coords_n = coords_p[n];
    detail::CompactCoordinates& coords_m = coords_p[m];
    std::vector<double> p_pp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < coords_n.get_size(); i++) {
        unsigned int j = 0;
        for (; j+7 < coords_m.get_size(); j+=8) {
            auto res = coords_n[i].evaluate(coords_m[j], coords_m[j+1], coords_m[j+2], coords_m[j+3], coords_m[j+4], coords_m[j+5], coords_m[j+6], coords_m[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < coords_m.get_size(); j+=4) {
            auto res = coords_n[i].evaluate(coords_m[j], coords_m[j+1], coords_m[j+2], coords_m[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < coords_m.get_size(); ++j) {
            auto res = coords_n[i].evaluate(coords_m[j]);
            p_pp[res.distance] += 2*res.weight;
        }
    }
    master -= partials_pp.index(n, m);
    partials_pp.index(n, m).get_counts() = std::move(p_pp);
    master += partials_pp.index(n, m);
}

void PartialHistogramManager::calc_hp(unsigned int index) {
    std::vector<double> p_hp(master.get_axis().bins, 0);

    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.get_size(); i++) {
        unsigned int j = 0;
        for (; j+7 < coords_h.get_size(); j+=8) {
            auto res = coords[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3], coords_h[j+4], coords_h[j+5], coords_h[j+6], coords_h[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_hp[res.distance[k]] += res.weight[k];}
        }

        for (; j+3 < coords_h.get_size(); j+=4) {
            auto res = coords[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_hp[res.distance[k]] += res.weight[k];}
        }

        for (; j < coords_h.get_size(); ++j) {
            auto res = coords[i].evaluate(coords_h[j]);
            p_hp[res.distance] += res.weight;
        }
    }

    master -= partials_hp.index(index)*2; // subtract the previous hydration histogram
    partials_hp.index(index).get_counts() = std::move(p_hp);
    master += partials_hp.index(index)*2; // add the new hydration histogram
}

void PartialHistogramManager::calc_hh() {
    std::vector<double> p_hh(master.get_axis().bins, 0);

    // calculate internal distances for the hydration layer
    coords_h = detail::CompactCoordinates(protein->get_waters()); //! Remove?
    for (unsigned int i = 0; i < coords_h.get_size(); i++) {
        unsigned int j = i+1;
        for (; j+7 < coords_h.get_size(); j+=8) {
            auto res = coords_h[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3], coords_h[j+4], coords_h[j+5], coords_h[j+6], coords_h[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < coords_h.get_size(); j+=4) {
            auto res = coords_h[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_hh[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < coords_h.get_size(); ++j) {
            auto res = coords_h[i].evaluate(coords_h[j]);
            p_hh[res.distance] += 2*res.weight;
        }
    }

    // calculate self-correlation
    p_hh[0] = std::accumulate(coords_h.get_data().begin(), coords_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.get_counts() = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
}