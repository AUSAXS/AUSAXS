#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/state/StateManager.h>
#include <settings/HistogramSettings.h>

using namespace hist;

PartialHistogramManager::PartialHistogramManager(Protein* protein) : HistogramManager(protein), coords_p(body_size), partials_pp(body_size, body_size), partials_hp(body_size) {}

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
    double width = settings::axes::distance_bin_width;
    detail::CompactCoordinates current(protein->get_body(index));

    // calculate internal distances between atoms
    std::vector<double> p_pp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < current.get_size(); i++) {
        for (unsigned int j = i+1; j < current.get_size(); j++) {
            float weight = current[i].w*current[j].w;
            float dx = current[i].x - current[j].x;
            float dy = current[i].y - current[j].y;
            float dz = current[i].z - current[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < current.get_size(); i++) {p_pp[0] += current[i].w*current[i].w;}

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
    double width = settings::axes::distance_bin_width;
    Axis axis(0, settings::axes::max_distance, settings::axes::max_distance/width); 
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
    double width = settings::axes::distance_bin_width;

    detail::CompactCoordinates& coords_n = coords_p[n];
    detail::CompactCoordinates& coords_m = coords_p[m];
    std::vector<double> p_pp(master.get_axis().bins, 0);
    for (unsigned int i = 0; i < coords_n.get_size(); i++) {
        for (unsigned int j = 0; j < coords_m.get_size(); j++) {
            float weight = coords_n[i].w*coords_m[j].w;
            float dx = coords_n[i].x - coords_m[j].x;
            float dy = coords_n[i].y - coords_m[j].y;
            float dz = coords_n[i].z - coords_m[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }
    master -= partials_pp.index(n, m);
    partials_pp.index(n, m).get_counts() = std::move(p_pp);
    master += partials_pp.index(n, m);
}

void PartialHistogramManager::calc_pp(unsigned int index) {
    double width = settings::axes::distance_bin_width;
    detail::CompactCoordinates& coords_i = coords_p[index];

    // we do not want to calculate the self-correlation, so we have to skip entry 'index'
    for (unsigned int n = 0; n < index; n++) { // loop from (0, index]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.get_axis().bins, 0);
        for (unsigned int i = 0; i < coords_i.get_size(); i++) {
            for (unsigned int j = 0; j < coords_j.get_size(); j++) {
                float weight = coords_i[i].w*coords_j[j].w;
                float dx = coords_i[i].x - coords_j[j].x;
                float dy = coords_i[i].y - coords_j[j].y;
                float dz = coords_i[i].z - coords_j[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp.index(index, n);
        partials_pp.index(index, n).get_counts() = std::move(p_pp);
        master += partials_pp.index(index, n);
    }

    for (unsigned int n = index+1; n < body_size; n++) { // loop from (index, size]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.get_axis().bins, 0);
        for (unsigned int i = 0; i < coords_i.get_size(); i++) {
            for (unsigned int j = 0; j < coords_j.get_size(); j++) {
                float weight = coords_i[i].w*coords_j[j].w;
                float dx = coords_i[i].x - coords_j[j].x;
                float dy = coords_i[i].y - coords_j[j].y;
                float dz = coords_i[i].z - coords_j[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp.index(index, n);
        partials_pp.index(index, n).get_counts() = std::move(p_pp);
        master += partials_pp.index(index, n);
    }
}

void PartialHistogramManager::calc_hp(unsigned int index) {
    double width = settings::axes::distance_bin_width;
    std::vector<double> p_hp(master.get_axis().bins, 0);

    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.get_size(); i++) {
        for (unsigned int j = 0; j < coords_h.get_size(); j++) {
            float weight = coords[i].w*coords_h[j].w;
            float dx = coords[i].x - coords_h[j].x;
            float dy = coords[i].y - coords_h[j].y;
            float dz = coords[i].z - coords_h[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += weight;
        }
    }

    master -= partials_hp.index(index)*2; // subtract the previous hydration histogram
    partials_hp.index(index).get_counts() = std::move(p_hp);
    master += partials_hp.index(index)*2; // add the new hydration histogram
}

void PartialHistogramManager::calc_hh() {
    const double& width = settings::axes::distance_bin_width;
    std::vector<double> p_hh(master.get_axis().bins, 0);

    // calculate internal distances for the hydration layer
    coords_h = detail::CompactCoordinates(protein->get_waters()); //! Remove?
    for (unsigned int i = 0; i < coords_h.get_size(); i++) {
        for (unsigned int j = i+1; j < coords_h.get_size(); j++) {
            float weight = coords_h[i].w*coords_h[j].w;
            float dx = coords_h[i].x - coords_h[j].x;
            float dy = coords_h[i].y - coords_h[j].y;
            float dz = coords_h[i].z - coords_h[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < coords_h.get_size(); i++) {p_hh[0] += coords_h[i].w*coords_h[i].w;}

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.get_counts() = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
}