#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>
#include <constants/Axes.h>
#include <utility/MultiThreading.h>
#include <container/ThreadLocalWrapper.h>

#include <mutex>
#include <list>

using namespace hist;

template<bool use_weighted_distribution> 
PartialHistogramManagerMT<use_weighted_distribution>::PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein) : PartialHistogramManager<use_weighted_distribution>(protein) {}

template<bool use_weighted_distribution> 
PartialHistogramManagerMT<use_weighted_distribution>::~PartialHistogramManagerMT() = default;

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> PartialHistogramManagerMT<use_weighted_distribution>::calculate() {
    std::vector<BS::multi_future<std::vector<double>>> futures_self_corr;
    std::vector<BS::multi_future<std::vector<double>>> futures_pp;
    std::vector<BS::multi_future<std::vector<double>>> futures_hp;
    BS::multi_future<std::vector<double>> futures_hh;
    futures_self_corr.reserve(this->body_size);
    futures_pp.reserve(this->body_size*(this->body_size-1)/2);
    futures_hp.reserve(this->body_size);

    const std::vector<bool>& externally_modified = this->statemanager->get_externally_modified_bodies();
    const std::vector<bool>& internally_modified = this->statemanager->get_internally_modified_bodies();
    const bool hydration_modified = this->statemanager->get_modified_hydration();
    auto pool = utility::multi_threading::get_global_pool();

    // check if the object has already been initialized
    if (this->master.get_counts().size() == 0) [[unlikely]] {
        initialize(); 
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < this->body_size; ++i) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[i]) {
                futures_self_corr.push_back(calc_self_correlation(i));
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations
            else if (externally_modified[i]) {
                pool->detach_task(
                    [&update_compact_representation_body, i] () {update_compact_representation_body(i);}
                );
            }
        }

        // merge the partial results from each thread and add it to the master histogram
        unsigned int counter = 0;
        for (unsigned int i = 0; i < this->body_size; ++i) {
            if (internally_modified[i]) {
                pool->detach_task(
                    [&combine_self_correlation, &futures_self_corr, i] () {combine_self_correlation(i, futures_self_corr);}
                );
            }
        }
    }

    // small efficiency improvement: if the hydration layer was modified, we can update the compact representations in parallel with the self-correlation
    if (hydration_modified) {
        pool->detach_task(
            update_compact_representation_water()
        );
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    // check if the hydration layer was modified
    if (hydration_modified) {
        futures_hh = calc_hh();
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (unsigned int i = 0; i < this->body_size; ++i) {
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
        for (unsigned int i = 0; i < this->body_size; ++i) {
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

    this->statemanager->reset();
    pool->wait();
    Distribution1D p = this->master.get_counts();
    return std::make_unique<DistanceHistogram>(std::move(p), this->master.get_axis());
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::update_compact_representation_body(unsigned int index) {
    this->coords_p[index] = detail::CompactCoordinates(this->protein->get_body(index));
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::update_compact_representation_water() {
    this->coords_h = detail::CompactCoordinates(this->protein->get_waters());
}

template<bool use_weighted_distribution> 
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManagerMT<use_weighted_distribution>::calculate_all() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;

    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_hh = this->partials_hh.get_counts();
    GenericDistribution1D_t p_pp = this->master.base.get_counts();
    GenericDistribution1D_t p_hp(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            detail::PartialHistogram& current = this->partials_pp.index(i, j);

            // iterate through each entry in the partial histogram
            std::transform(p_pp.begin(), p_pp.end(), current.get_counts().begin(), p_pp.begin(), std::plus<>());
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; ++i) {
        detail::PartialHistogram& current = this->partials_hp.index(i);

        // iterate through each entry in the partial histogram
        std::transform(p_hp.begin(), p_hp.end(), current.get_counts().begin(), p_hp.begin(), std::plus<>());
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
void PartialHistogramManagerMT<use_weighted_distribution>::initialize() {
    auto pool = utility::multi_threading::get_global_pool();
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram(p_base, axis);

    static std::vector<BS::multi_future<std::vector<double>>> futures;
    futures = std::vector<BS::multi_future<std::vector<double>>>(this->body_size);
    this->partials_hh = detail::PartialHistogram(axis);
    for (unsigned int i = 0; i < this->body_size; ++i) {
        this->partials_hp.index(i) = detail::PartialHistogram(axis);
        this->partials_pp.index(i, i) = detail::PartialHistogram(axis);
        futures[i] = calc_self_correlation(i);

        for (unsigned int j = 0; j < i; ++j) {
            this->partials_pp.index(i, j) = detail::PartialHistogram(axis);
        }
    }

    for (unsigned int i = 0; i < this->body_size; ++i) {
        pool->push_task(&PartialHistogramManagerMT::combine_self_correlation, this, i, std::ref(futures[i]));
    }
}

template<bool use_weighted_distribution> 
BS::multi_future<std::vector<double>> PartialHistogramManagerMT<use_weighted_distribution>::calc_self_correlation(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();
    update_compact_representation_body(index);

    // calculate internal distances between atoms
    static auto calc_internal = [] (const detail::CompactCoordinates& coords, unsigned int pp_size, unsigned int imin, unsigned int imax) {
        GenericDistribution1D_t p_pp(pp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = i+1;
            for (; j+7 < coords.size(); j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_pp, coords, coords, i, j);
            }

            for (; j+3 < coords.size(); j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_pp, coords, coords, i, j);
            }

            for (; j < coords.size(); ++j) {
                evaluate1<use_weighted_distribution, 2>(p_pp, coords, coords, i, j);
            }
        }
        return p_pp.get_data();
    };

    // calculate self correlation
    static auto calc_self = [] (const detail::CompactCoordinates& coords, unsigned int pp_size) {
        GenericDistribution1D_t p_pp(pp_size, 0);
        p_pp.add(0, std::accumulate(coords.get_data().begin(), coords.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));
        return p_pp.get_data();
    };

    BS::multi_future<std::vector<double>> futures_self_corr;
    unsigned int atom_size = this->protein->atom_size();
    for (unsigned int i = 0; i < atom_size; i += settings::general::detail::job_size) {
        futures_self_corr.push_back(pool->submit(calc_internal, std::cref(this->coords_p[index]), this->master.get_axis().bins, i, std::min(i+settings::general::detail::job_size, atom_size)));
    }
    futures_self_corr.push_back(pool->submit(calc_self, std::cref(this->coords_p[index]), this->master.get_axis().bins));
    return futures_self_corr;
}

template<bool use_weighted_distribution> 
BS::multi_future<std::vector<double>> PartialHistogramManagerMT<use_weighted_distribution>::calc_pp(unsigned int n, unsigned int m) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_pp = [] (const detail::CompactCoordinates& coords_n, const detail::CompactCoordinates& coords_m, unsigned int pp_size, unsigned int imin, unsigned int imax) {
        GenericDistribution1D_t p_pp(pp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = 0;
            for (; j+7 < coords_m.size(); j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_pp, coords_n, coords_m, i, j);
            }

            for (; j+3 < coords_m.size(); j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_pp, coords_n, coords_m, i, j);
            }

            for (; j < coords_m.size(); ++j) {
                evaluate1<use_weighted_distribution, 2>(p_pp, coords_n, coords_m, i, j);
            }
        }
        return p_pp.get_data();
    };

    BS::multi_future<std::vector<double>> futures_pp;
    detail::CompactCoordinates& coords_n = this->coords_p[n];
    for (unsigned int i = 0; i < coords_n.size(); i += settings::general::detail::job_size) {
        futures_pp.push_back(pool->submit(calc_pp, std::cref(this->coords_p[n]), std::cref(this->coords_p[m]), this->master.get_axis().bins, i, std::min<int>(i+settings::general::detail::job_size, coords_n.size())));
    }
    return futures_pp;
}

template<bool use_weighted_distribution> 
BS::multi_future<std::vector<double>> PartialHistogramManagerMT<use_weighted_distribution>::calc_hp(unsigned int index) {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_hp = [] (const detail::CompactCoordinates& coords_i, const detail::CompactCoordinates& coords_h, unsigned int hp_size, unsigned int imin, unsigned int imax) {
        GenericDistribution1D_t p_hp(hp_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = 0;
            for (; j+7 < coords_h.size(); j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_hp, coords_i, coords_h, i, j);
            }

            for (; j+3 < coords_h.size(); j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_hp, coords_i, coords_h, i, j);
            }

            for (; j < coords_h.size(); ++j) {
                evaluate1<use_weighted_distribution, 1>(p_hp, coords_i, coords_h, i, j);
            }
        }
        return p_hp.get_data();
    };

    BS::multi_future<std::vector<double>> futures_hp;
    detail::CompactCoordinates& coords = this->coords_p[index];
    for (unsigned int i = 0; i < coords.size(); i += settings::general::detail::job_size) {
        futures_hp.push_back(pool->submit(calc_hp, std::cref(this->coords_p[index]), std::cref(this->coords_h), this->master.get_axis().bins, i, std::min<int>(i+settings::general::detail::job_size, coords.size())));
    }
    return futures_hp;
}

template<bool use_weighted_distribution> 
BS::multi_future<std::vector<double>> PartialHistogramManagerMT<use_weighted_distribution>::calc_hh() {
    using GenericDistribution1D_t = typename hist::GenericDistribution1D<use_weighted_distribution>::type;
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_hh = [] (const detail::CompactCoordinates& coords_h, unsigned int hh_size, unsigned int imin, unsigned int imax) {
        GenericDistribution1D_t p_hh(hh_size, 0);
        for (unsigned int i = imin; i < imax; ++i) {
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
        return p_hh.get_data();
    };

    // calculate self correlation
    static auto calc_self = [] (const detail::CompactCoordinates& coords_h, unsigned int hh_size) {
        GenericDistribution1D_t p_hh(hh_size, 0);
        p_hh.add(0, std::accumulate(coords_h.get_data().begin(), coords_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));
        return p_hh.get_data();
    };

    BS::multi_future<std::vector<double>> futures_hh;
    for (unsigned int i = 0; i < this->coords_h.size(); i += settings::general::detail::job_size) {
        futures_hh.push_back(pool->submit(calc_hh, std::cref(this->coords_h), this->master.get_axis().bins, i, std::min<int>(i+settings::general::detail::job_size, this->coords_h.size())));
    }
    futures_hh.push_back(pool->submit(calc_self, std::cref(this->coords_h), this->master.get_axis().bins));
    return futures_hh;
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_self_correlation(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(this->master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    this->master -= this->partials_pp.index(index, index);
    this->partials_pp.index(index, index).get_counts() = std::move(p_pp);
    this->master += this->partials_pp.index(index, index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_pp(unsigned int n, unsigned int m, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_pp(this->master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_pp.begin(), p_pp.end(), tmp.begin(), p_pp.begin(), std::plus<>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    this->master -= this->partials_pp.index(n, m);
    this->partials_pp.index(n, m).get_counts() = std::move(p_pp);
    this->master += this->partials_pp.index(n, m);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_hp(unsigned int index, BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hp(this->master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hp.begin(), p_hp.end(), tmp.begin(), p_hp.begin(), std::plus<>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    this->master -= this->partials_hp.index(index)*2; // subtract the previous hydration histogram
    this->partials_hp.index(index).get_counts() = std::move(p_hp);
    this->master += this->partials_hp.index(index)*2; // add the new hydration histogram
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_hh(BS::multi_future<std::vector<double>>& futures) {
    std::vector<double> p_hh(this->master.get_axis().bins, 0);

    // iterate through all partial results. Each partial result is from a different thread calculation
    for (auto& tmp : futures.get()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.begin(), p_hh.begin(), std::plus<>());
    }

    // update the master histogram
    master_hist_mutex.lock();
    this->master -= this->partials_hh; // subtract the previous hydration histogram
    this->partials_hh.get_counts() = std::move(p_hh);
    this->master += this->partials_hh; // add the new hydration histogram
    master_hist_mutex.unlock();
}

template class hist::PartialHistogramManagerMT<true>;
template class hist::PartialHistogramManagerMT<false>;