#include <hist/PartialHistogramManager.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/state/StateManager.h>
#include <settings/HistogramSettings.h>

using namespace hist;

PartialHistogramManager::PartialHistogramManager(Protein* protein) : HistogramManager(protein), coords_p(size), partials_pp(size, size), partials_hp(size) {}

PartialHistogramManager::~PartialHistogramManager() = default;

std::unique_ptr<DistanceHistogram>  PartialHistogramManager::calculate() {
    const std::vector<bool> externally_modified = statemanager->get_externally_modified_bodies();
    const std::vector<bool> internally_modified = statemanager->get_internally_modified_bodies();

    // check if the object has already been initialized
    if (master.p.size() == 0) [[unlikely]] {
        initialize(); 
    } 
    
    // if not, we must first check if the coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < size; i++) {
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
        for (unsigned int i = 0; i < size; i++) {
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
        for (unsigned int i = 0; i < size; i++) {
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
    return std::make_unique<DistanceHistogram>(master.p.copy(), master.axis);
}

std::unique_ptr<CompositeDistanceHistogram> PartialHistogramManager::calculate_all() {
    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    std::vector<double> p_hh = partials_hh.p;
    std::vector<double> p_pp = master.base.p;
    std::vector<double> p_hp(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < i; j++) {
            detail::PartialHistogram& current = partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < bins; k++) {
                p_pp[k] += current.p[k]; // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < size; i++) {
        detail::PartialHistogram& current = partials_hp.index(i);

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < bins; k++) {
            p_hp[k] += current.p[k]; // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(bins);
    p_pp.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(total->p), total->get_axis());
}

void PartialHistogramManager::calc_self_correlation(unsigned int index) {
    double width = settings::axes::distance_bin_width;
    detail::CompactCoordinates current(protein->get_body(index));

    // calculate internal distances between atoms
    std::vector<double> p_pp(master.axis.bins, 0);
    for (unsigned int i = 0; i < current.size; i++) {
        for (unsigned int j = i+1; j < current.size; j++) {
            float weight = current.data[i].w*current.data[j].w;
            float dx = current.data[i].x - current.data[j].x;
            float dy = current.data[i].y - current.data[j].y;
            float dz = current.data[i].z - current.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < current.size; i++) {p_pp[0] += current.data[i].w*current.data[i].w;}

    // store the coordinates for later
    coords_p[index] = std::move(current);

    master.base -= partials_pp.index(index, index);
    master -= partials_pp.index(index, index);
    partials_pp.index(index, index).p = std::move(p_pp);
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
    for (unsigned int n = 0; n < size; n++) {
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
    std::vector<double> p_pp(master.axis.bins, 0);
    for (unsigned int i = 0; i < coords_n.size; i++) {
        for (unsigned int j = 0; j < coords_m.size; j++) {
            float weight = coords_n.data[i].w*coords_m.data[j].w;
            float dx = coords_n.data[i].x - coords_m.data[j].x;
            float dy = coords_n.data[i].y - coords_m.data[j].y;
            float dz = coords_n.data[i].z - coords_m.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_pp[dist/width] += 2*weight;
        }
    }
    master -= partials_pp.index(n, m);
    partials_pp.index(n, m).p = std::move(p_pp);
    master += partials_pp.index(n, m);
}

void PartialHistogramManager::calc_pp(unsigned int index) {
    double width = settings::axes::distance_bin_width;
    detail::CompactCoordinates& coords_i = coords_p[index];

    // we do not want to calculate the self-correlation, so we have to skip entry 'index'
    for (unsigned int n = 0; n < index; n++) { // loop from (0, index]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = 0; i < coords_i.size; i++) {
            for (unsigned int j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[i].w*coords_j.data[j].w;
                float dx = coords_i.data[i].x - coords_j.data[j].x;
                float dy = coords_i.data[i].y - coords_j.data[j].y;
                float dz = coords_i.data[i].z - coords_j.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp.index(index, n);
        partials_pp.index(index, n).p = std::move(p_pp);
        master += partials_pp.index(index, n);
    }

    for (unsigned int n = index+1; n < size; n++) { // loop from (index, size]
        detail::CompactCoordinates& coords_j = coords_p[n];
        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = 0; i < coords_i.size; i++) {
            for (unsigned int j = 0; j < coords_j.size; j++) {
                float weight = coords_i.data[i].w*coords_j.data[j].w;
                float dx = coords_i.data[i].x - coords_j.data[j].x;
                float dy = coords_i.data[i].y - coords_j.data[j].y;
                float dz = coords_i.data[i].z - coords_j.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        master -= partials_pp.index(index, n);
        partials_pp.index(index, n).p = std::move(p_pp);
        master += partials_pp.index(index, n);
    }
}

void PartialHistogramManager::calc_hp(unsigned int index) {
    double width = settings::axes::distance_bin_width;
    std::vector<double> p_hp(master.axis.bins, 0);

    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.size; i++) {
        for (unsigned int j = 0; j < coords_h.size; j++) {
            float weight = coords.data[i].w*coords_h.data[j].w;
            float dx = coords.data[i].x - coords_h.data[j].x;
            float dy = coords.data[i].y - coords_h.data[j].y;
            float dz = coords.data[i].z - coords_h.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hp[dist/width] += 2*weight;
        }
    }

    master -= partials_hp.index(index); // subtract the previous hydration histogram
    partials_hp.index(index).p = std::move(p_hp);
    master += partials_hp.index(index); // add the new hydration histogram
}

void PartialHistogramManager::calc_hh() {
    const double& width = settings::axes::distance_bin_width;
    std::vector<double> p_hh(master.axis.bins, 0);

    // calculate internal distances for the hydration layer
    coords_h = detail::CompactCoordinates(protein->get_waters()); //! Remove?
    for (unsigned int i = 0; i < coords_h.size; i++) {
        for (unsigned int j = i+1; j < coords_h.size; j++) {
            float weight = coords_h.data[i].w*coords_h.data[j].w;
            float dx = coords_h.data[i].x - coords_h.data[j].x;
            float dy = coords_h.data[i].y - coords_h.data[j].y;
            float dz = coords_h.data[i].z - coords_h.data[j].z;
            float dist = sqrt(dx*dx + dy*dy + dz*dz);
            p_hh[dist/width] += 2*weight;
        }
    }

    // calculate self-correlation
    for (unsigned int i = 0; i < coords_h.size; i++) {p_hh[0] += coords_h.data[i].w*coords_h.data[i].w;}

    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.p = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
}