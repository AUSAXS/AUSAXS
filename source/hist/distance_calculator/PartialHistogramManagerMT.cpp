#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <hist/DistanceHistogram.h>
#include <hist/CompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/record/Atom.h> 
#include <data/record/Water.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>
#include <constants/Constants.h>

#include <mutex>
#include <list>
#include <BS_thread_pool.hpp>

using namespace hist;

PartialHistogramManagerMT::PartialHistogramManagerMT(view_ptr<const data::Molecule> protein) : PartialHistogramManager(protein), pool(std::make_unique<BS::thread_pool>(settings::general::threads)) {}

PartialHistogramManagerMT::PartialHistogramManagerMT(PartialHistogramManager& phm) : PartialHistogramManager(phm), pool(std::make_unique<BS::thread_pool>(settings::general::threads)) {}

PartialHistogramManagerMT::~PartialHistogramManagerMT() = default;

std::unique_ptr<DistanceHistogram> PartialHistogramManagerMT::calculate() {
    std::vector<BS::multi_future<std::vector<double>>> futures_self_corr;
    std::vector<BS::multi_future<std::vector<double>>> futures_pp;
    std::vector<BS::multi_future<std::vector<double>>> futures_hp;
    BS::multi_future<std::vector<double>> futures_hh;
    futures_self_corr.reserve(body_size);
    futures_pp.reserve(body_size*(body_size-1)/2);
    futures_hp.reserve(body_size);

    const std::vector<bool>& externally_modified = statemanager->get_externally_modified_bodies();
    const std::vector<bool>& internally_modified = statemanager->get_internally_modified_bodies();
    const bool hydration_modified = statemanager->get_modified_hydration();

    // check if the object has already been initialized
    if (master.get_counts().size() == 0) [[unlikely]] {
        initialize(); 
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < body_size; ++i) {

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
        for (unsigned int i = 0; i < body_size; ++i) {
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
    for (unsigned int i = 0; i < body_size; ++i) {
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
        for (unsigned int i = 0; i < body_size; ++i) {
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
    std::vector<double> p = master.get_counts();
    return std::make_unique<DistanceHistogram>(std::move(p), master.get_axis());
}

void PartialHistogramManagerMT::update_compact_representation_body(unsigned int index) {
    coords_p[index] = detail::CompactCoordinates(protein->get_body(index));
}

void PartialHistogramManagerMT::update_compact_representation_water() {
    coords_h = detail::CompactCoordinates(protein->get_waters());
}

std::unique_ptr<CompositeDistanceHistogram> PartialHistogramManagerMT::calculate_all() {
    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    std::vector<double> p_hh = partials_hh.get_counts();
    std::vector<double> p_pp = master.base.get_counts();
    std::vector<double> p_hp(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < body_size; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            detail::PartialHistogram& current = partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            for (unsigned int k = 0; k < bins; ++k) {
                p_pp[k] += current.get_count(k); // add to p_pp
            }
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < body_size; ++i) {
        detail::PartialHistogram& current = partials_hp.index(i);

        // iterate through each entry in the partial histogram
        for (unsigned int k = 0; k < bins; ++k) {
            p_hp[k] += current.get_count(k); // add to p_pp
        }
    }

    // p_hp is already resized
    p_hh.resize(bins);
    p_pp.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(std::move(p_pp), std::move(p_hp), std::move(p_hh), std::move(total->get_counts()), total->get_axis());
}

/**
 * @brief This initializes some necessary variables and precalculates the internal distances between atoms in each body.
 */
void PartialHistogramManagerMT::initialize() {
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    master = detail::MasterHistogram(p_base, axis);

    static std::vector<BS::multi_future<std::vector<double>>> futures;
    futures = std::vector<BS::multi_future<std::vector<double>>>(body_size);
    partials_hh = detail::PartialHistogram(axis);
    for (unsigned int i = 0; i < body_size; ++i) {
        partials_hp.index(i) = detail::PartialHistogram(axis);
        partials_pp.index(i, i) = detail::PartialHistogram(axis);
        futures[i] = calc_self_correlation(i);

        for (unsigned int j = 0; j < i; ++j) {
            partials_pp.index(i, j) = detail::PartialHistogram(axis);
        }
    }

    for (unsigned int i = 0; i < body_size; ++i) {
        pool->push_task(&PartialHistogramManagerMT::combine_self_correlation, this, i, std::ref(futures[i]));
    }
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_self_correlation(unsigned int index) {
    update_compact_representation_body(index);

    // calculate internal distances between atoms
    static auto calc_internal = [] (const detail::CompactCoordinates& coords, unsigned int pp_size, unsigned int imin, unsigned int imax) {
        std::vector<double> p_pp(pp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
        unsigned int j = i+1;
        for (; j+7 < coords.get_size(); j+=8) {
            auto res = coords[i].evaluate(coords[j], coords[j+1], coords[j+2], coords[j+3], coords[j+4], coords[j+5], coords[j+6], coords[j+7]);
            for (unsigned int k = 0; k < 8; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j+3 < coords.get_size(); j+=4) {
            auto res = coords[i].evaluate(coords[j], coords[j+1], coords[j+2], coords[j+3]);
            for (unsigned int k = 0; k < 4; ++k) {p_pp[res.distance[k]] += 2*res.weight[k];}
        }

        for (; j < coords.get_size(); ++j) {
            auto res = coords[i].evaluate(coords[j]);
            p_pp[res.distance] += 2*res.weight;
        }
        }
        return p_pp;
    };

    // calculate self correlation
    static auto calc_self = [] (const detail::CompactCoordinates& coords, unsigned int pp_size) {
        std::vector<double> p_pp(pp_size, 0);
        p_pp[0] = std::accumulate(coords.get_data().begin(), coords.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );
        return p_pp;
    };

    BS::multi_future<std::vector<double>> futures_self_corr;
    unsigned int atom_size = protein->atom_size();
    for (unsigned int i = 0; i < atom_size; i += settings::general::detail::job_size) {
        futures_self_corr.push_back(pool->submit(calc_internal, std::cref(coords_p[index]), master.get_axis().bins, i, std::min(i+settings::general::detail::job_size, atom_size)));
    }
    futures_self_corr.push_back(pool->submit(calc_self, std::cref(coords_p[index]), master.get_axis().bins));
    return futures_self_corr;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_pp(unsigned int n, unsigned int m) {
    static auto calc_pp = [] (const detail::CompactCoordinates& coords_n, const detail::CompactCoordinates& coords_m, unsigned int pp_size, unsigned int imin, unsigned int imax) {
        std::vector<double> p_pp(pp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
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
        return p_pp;
    };

    BS::multi_future<std::vector<double>> futures_pp;
    detail::CompactCoordinates& coords_n = coords_p[n];
    for (unsigned int i = 0; i < coords_n.get_size(); i += settings::general::detail::job_size) {
        futures_pp.push_back(pool->submit(calc_pp, std::cref(coords_p[n]), std::cref(coords_p[m]), master.get_axis().bins, i, std::min(i+settings::general::detail::job_size, coords_n.get_size())));
    }
    return futures_pp;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_hp(unsigned int index) {
    static auto calc_hp = [] (const detail::CompactCoordinates& coords_i, const detail::CompactCoordinates& coords_h, unsigned int hp_size, unsigned int imin, unsigned int imax) {
        std::vector<double> p_hp(hp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = 0;
            for (; j+7 < coords_h.get_size(); j+=8) {
                auto res = coords_i[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3], coords_h[j+4], coords_h[j+5], coords_h[j+6], coords_h[j+7]);
                for (unsigned int k = 0; k < 8; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j+3 < coords_h.get_size(); j+=4) {
                auto res = coords_i[i].evaluate(coords_h[j], coords_h[j+1], coords_h[j+2], coords_h[j+3]);
                for (unsigned int k = 0; k < 4; ++k) {p_hp[res.distance[k]] += res.weight[k];}
            }

            for (; j < coords_h.get_size(); ++j) {
                auto res = coords_i[i].evaluate(coords_h[j]);
                p_hp[res.distance] += res.weight;
            }
        }
        return p_hp;
    };

    BS::multi_future<std::vector<double>> futures_hp;
    detail::CompactCoordinates& coords = coords_p[index];
    for (unsigned int i = 0; i < coords.get_size(); i += settings::general::detail::job_size) {
        futures_hp.push_back(pool->submit(calc_hp, std::cref(coords_p[index]), std::cref(coords_h), master.get_axis().bins, i, std::min(i+settings::general::detail::job_size, coords.get_size())));
    }
    return futures_hp;
}

BS::multi_future<std::vector<double>> PartialHistogramManagerMT::calc_hh() {
    // calculate internal distances for the hydration layer
    static auto calc_hh = [] (const detail::CompactCoordinates& coords_h, unsigned int hh_size, unsigned int imin, unsigned int imax) {
        std::vector<double> p_hh(hh_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
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
        return p_hh;
    };

    // calculate self correlation
    static auto calc_self = [] (const detail::CompactCoordinates& coords_h, unsigned int hh_size) {
        std::vector<double> p_hh(hh_size, 0);
        p_hh[0] = std::accumulate(coords_h.get_data().begin(), coords_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;} );
        return p_hh;
    };

    BS::multi_future<std::vector<double>> futures_hh;
    for (unsigned int i = 0; i < coords_h.get_size(); i += settings::general::detail::job_size) {
        futures_hh.push_back(pool->submit(calc_hh, std::cref(coords_h), master.get_axis().bins, i, std::min(i+settings::general::detail::job_size, coords_h.get_size())));
    }
    futures_hh.push_back(pool->submit(calc_self, std::cref(coords_h), master.get_axis().bins));
    return futures_hh;
}

void PartialHistogramManagerMT::combine_self_correlation(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_pp.index(index, index);
    partials_pp.index(index, index).get_counts() = std::move(p_pp);
    master += partials_pp.index(index, index);
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_pp(unsigned int n, unsigned int m, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_pp.index(n, m);
    partials_pp.index(n, m).get_counts() = std::move(p_pp);
    master += partials_pp.index(n, m);
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_hp(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hp(master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_hp.index(index)*2; // subtract the previous hydration histogram
    partials_hp.index(index).get_counts() = std::move(p_hp);
    master += partials_hp.index(index)*2; // add the new hydration histogram
    master_hist_mutex.unlock();
}

void PartialHistogramManagerMT::combine_hh(BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hh(master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.begin(), p_hh.begin(), std::plus<double>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    master -= partials_hh; // subtract the previous hydration histogram
    partials_hh.get_counts() = std::move(p_hh);
    master += partials_hh; // add the new hydration histogram
    master_hist_mutex.unlock();
}