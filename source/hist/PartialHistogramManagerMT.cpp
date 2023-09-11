#include <hist/PartialHistogramManagerMT.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h> 
#include <data/Water.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>

#include <mutex>
#include <list>
#include <BS_thread_pool.hpp>

using namespace hist;

PartialHistogramManagerMT::PartialHistogramManagerMT(Protein* protein) : PartialHistogramManager(protein), pool(std::make_unique<BS::thread_pool>(settings::general::threads)) {}

PartialHistogramManagerMT::PartialHistogramManagerMT(PartialHistogramManager& phm) : PartialHistogramManager(phm), pool(std::make_unique<BS::thread_pool>(settings::general::threads)) {}

PartialHistogramManagerMT::~PartialHistogramManagerMT() = default;

Histogram PartialHistogramManagerMT::calculate() {
    std::vector<BS::multi_future<std::vector<double>>> futures_self_corr;
    std::vector<BS::multi_future<std::vector<double>>> futures_pp;
    std::vector<BS::multi_future<std::vector<double>>> futures_hp;
    BS::multi_future<std::vector<double>> futures_hh;
    futures_self_corr.reserve(size);
    futures_pp.reserve(size*(size-1)/2);
    futures_hp.reserve(size);

    const std::vector<bool> externally_modified = statemanager->get_externally_modified_bodies();
    const std::vector<bool> internally_modified = statemanager->get_internally_modified_bodies();
    const bool hydration_modified = statemanager->get_modified_hydration();

    // check if the object has already been initialized
    if (master.p.size() == 0) [[unlikely]] {
        initialize(); 
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < size; ++i) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[i]) {
                futures_self_corr.push_back(calc_self_correlation(i));
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations
            else if (externally_modified[i]) {
                pool->push_task(&PartialHistogramManagerMT::update_compact_representation_body, this, i);
            }
        }

        // merge the partial results from each thread and add it to the master histogram
        unsigned int counter = 0;
        for (unsigned int i = 0; i < size; ++i) {
            if (internally_modified[i]) {
                pool->push_task(&PartialHistogramManagerMT::combine_self_correlation, this, i, std::ref(futures_self_corr[counter++]));
            }
        }
    }

    // small efficiency improvement: if the hydration layer was modified, we can update the compact representations in parallel with the self-correlation
    if (hydration_modified) {
        pool->push_task(&PartialHistogramManagerMT::update_compact_representation_water, this);
    }
    pool->wait_for_tasks(); // ensure the compact representations have been updated before continuing

    // check if the hydration layer was modified
    if (hydration_modified) {
        futures_hh = calc_hh();
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (unsigned int i = 0; i < size; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            if (externally_modified[i] || externally_modified[j]) {
                // one of the bodies was modified, so we recalculate its partial histogram
                futures_pp.push_back(calc_pp(i, j));
            }
        }

        // we also have to remember to update the partial histograms with the hydration layer
        if (externally_modified[i] || hydration_modified) {
            futures_hp.push_back(calc_hp(i));
        }
    }

    // merge the partial results from each thread and add it to the master histogram
    {
        if (hydration_modified) {
            pool->push_task(&PartialHistogramManagerMT::combine_hh, this, std::ref(futures_hh));
        }

        unsigned int counter_pp = 0, counter_hp = 0;
        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                if (externally_modified[i] || externally_modified[j]) {
                    pool->push_task(&PartialHistogramManagerMT::combine_pp, this, i, j, std::ref(futures_pp[counter_pp++]));
                }
            }

            if (externally_modified[i] || hydration_modified) {
                pool->push_task(&PartialHistogramManagerMT::combine_hp, this, i, std::ref(futures_hp[counter_hp++]));
            }
        }
    }

    statemanager->reset();
    pool->wait_for_tasks();
    return Histogram(master.p, master.axis);
}

void PartialHistogramManagerMT::update_compact_representation_body(unsigned int index) {
    coords_p[index] = detail::CompactCoordinates(protein->get_body(index));
}

void PartialHistogramManagerMT::update_compact_representation_water() {
    coords_h = detail::CompactCoordinates(protein->get_waters());
}

ScatteringHistogram PartialHistogramManagerMT::calculate_all() {
    Histogram total = calculate();
    total.shorten_axis();

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    std::vector<double> p_hh = partials_hh.p;
    std::vector<double> p_pp = master.base.p;
    std::vector<double> p_hp(total.axis.bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < size; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            detail::PartialHistogram& current = partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < total.axis.bins; k++) {
                p_pp[k] += current.p[k]; // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < size; ++i) {
        detail::PartialHistogram& current = partials_hp.index(i);

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < total.axis.bins; k++) {
            p_hp[k] += current.p[k]; // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(total.axis.bins);
    p_pp.resize(total.axis.bins);

    return ScatteringHistogram(p_pp, p_hh, p_hp, std::move(total.p), total.axis);
}

/**
 * @brief This initializes some necessary variables and precalculates the internal distances between atoms in each body.
 */
void PartialHistogramManagerMT::initialize() {
    double width = settings::axes::distance_bin_width;
    Axis axis(0, settings::axes::max_distance, settings::axes::max_distance/width); 
    std::vector<double> p_base(axis.bins, 0);
    master = detail::MasterHistogram(p_base, axis);

    static std::vector<BS::multi_future<std::vector<double>>> futures;
    futures = std::vector<BS::multi_future<std::vector<double>>>(size);
    partials_hh = detail::PartialHistogram(axis);
    for (unsigned int i = 0; i < size; ++i) {
        partials_hp.index(i) = detail::PartialHistogram(axis);
        partials_pp.index(i, i) = detail::PartialHistogram(axis);
        futures[i] = calc_self_correlation(i);

        for (unsigned int j = 0; j < i; j++) {
            partials_pp.index(i, j) = detail::PartialHistogram(axis);
        }
    }

    for (unsigned int i = 0; i < size; ++i) {
        pool->push_task(&PartialHistogramManagerMT::combine_self_correlation, this, i, std::ref(futures[i]));
    }
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_self_correlation(unsigned int index) {
    update_compact_representation_body(index);

    // calculate internal distances between atoms
    static auto calc_internal = [this] (unsigned int index, unsigned int imin, unsigned int imax) {
        double width = settings::axes::distance_bin_width;
        detail::CompactCoordinates& current = coords_p[index];

        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = i+1; j < current.size; ++j) {
                float weight = current.data[i].w*current.data[j].w;
                float dx = current.data[i].x - current.data[j].x;
                float dy = current.data[i].y - current.data[j].y;
                float dz = current.data[i].z - current.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        return p_pp;
    };

    // calculate self correlation
    static auto calc_self = [this] (unsigned int index) {
        detail::CompactCoordinates& current = coords_p[index];
        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = 0; i < current.size; ++i) {
            p_pp[0] += current.data[i].w*current.data[i].w;
        }
        return p_pp;
    };

    BS::multi_future<std::vector<double>> futures_self_corr;
    for (unsigned int i = 0; i < size; i += settings::general::detail::job_size) {
        futures_self_corr.push_back(pool->submit(calc_internal, index, i, std::min(i+settings::general::detail::job_size, size)));
    }
    futures_self_corr.push_back(pool->submit(calc_self, index));
    return futures_self_corr;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_pp(unsigned int n, unsigned int m) {
    static auto calc_pp = [this] (unsigned int n, unsigned int m, unsigned int imin, unsigned int imax) {
        double width = settings::axes::distance_bin_width;
        detail::CompactCoordinates& coords_n = coords_p[n];
        detail::CompactCoordinates& coords_m = coords_p[m];

        std::vector<double> p_pp(master.axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = 0; j < coords_m.size; ++j) {
                float weight = coords_n.data[i].w*coords_m.data[j].w;
                float dx = coords_n.data[i].x - coords_m.data[j].x;
                float dy = coords_n.data[i].y - coords_m.data[j].y;
                float dz = coords_n.data[i].z - coords_m.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_pp[dist/width] += 2*weight;
            }
        }
        return p_pp;
    };

    BS::multi_future<std::vector<double>> futures_pp;
    detail::CompactCoordinates& coords_n = coords_p[n];
    for (unsigned int i = 0; i < coords_n.size; i += settings::general::detail::job_size) {
        futures_pp.push_back(pool->submit(calc_pp, n, m, i, std::min(i+settings::general::detail::job_size, coords_n.size)));
    }
    return futures_pp;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_hp(unsigned int index) {
    static auto calc_hp = [this] (unsigned int index, unsigned int imin, unsigned int imax) {
        double width = settings::axes::distance_bin_width;
        detail::CompactCoordinates& coords = coords_p[index];

        std::vector<double> p_hp(master.axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = 0; j < coords_h.size; ++j) {
                float weight = coords.data[i].w*coords_h.data[j].w;
                float dx = coords.data[i].x - coords_h.data[j].x;
                float dy = coords.data[i].y - coords_h.data[j].y;
                float dz = coords.data[i].z - coords_h.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_hp[dist/width] += 2*weight;
            }
        }
        return p_hp;
    };

    BS::multi_future<std::vector<double>> futures_hp;
    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.size; i += settings::general::detail::job_size) {
        futures_hp.push_back(pool->submit(calc_hp, index, i, std::min(i+settings::general::detail::job_size, coords.size)));
    }
    return futures_hp;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_hh() {
    // calculate internal distances for the hydration layer
    static auto calc_hh = [this] (unsigned int imin, unsigned int imax) {
        double width = settings::axes::distance_bin_width;

        std::vector<double> p_hh(master.axis.bins, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            for (unsigned int j = i+1; j < protein->get_waters().size(); ++j) {
                float weight = coords_h.data[i].w*coords_h.data[j].w;
                float dx = coords_h.data[i].x - coords_h.data[j].x;
                float dy = coords_h.data[i].y - coords_h.data[j].y;
                float dz = coords_h.data[i].z - coords_h.data[j].z;
                float dist = sqrt(dx*dx + dy*dy + dz*dz);
                p_hh[dist/width] += 2*weight;
            }
        }
        return p_hh;
    };

    // calculate self correlation
    static auto calc_self = [this] () {
        std::vector<double> p_hh(master.axis.bins, 0);
        for (unsigned int i = 0; i < protein->get_waters().size(); ++i) {
            p_hh[0] += coords_h.data[i].w*coords_h.data[i].w;
        }
        return p_hh;
    };

    BS::multi_future<std::vector<double>> futures_hh;
    for (unsigned int i = 0; i < coords_h.size; i += settings::general::detail::job_size) {
        futures_hh.push_back(pool->submit(calc_hh, i, std::min(i+settings::general::detail::job_size, coords_h.size)));
    }
    futures_hh.push_back(pool->submit(calc_self));
    return futures_hh;
}

void PartialHistogramManagerMT::combine_self_correlation(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(master.axis.bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master.base -= partials_pp.index(index, index);
    master -= partials_pp.index(index, index);
    partials_pp.index(index, index).p = std::move(p_pp);
    master += partials_pp.index(index, index);
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_pp(unsigned int n, unsigned int m, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(master.axis.bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_pp.index(n, m);
    partials_pp.index(n, m).p = std::move(p_pp);
    master += partials_pp.index(n, m);
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_hp(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hp(master.axis.bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_hp.index(index); // subtract the previous hydration histogram
    partials_hp.index(index).p = std::move(p_hp);
    master += partials_hp.index(index); // add the new hydration histogram
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_hh(BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hh(master.axis.bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.begin(), p_hh.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.p = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
    master_hist_mutex.unlock();
}