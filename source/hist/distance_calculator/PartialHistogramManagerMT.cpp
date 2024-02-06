#include "data/Molecule.h"
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

using namespace hist;

template<bool use_weighted_distribution> 
PartialHistogramManagerMT<use_weighted_distribution>::PartialHistogramManagerMT(observer_ptr<const data::Molecule> protein) 
    : PartialHistogramManager<use_weighted_distribution>(protein), partials_pp_all(this->body_size, this->body_size), partials_hp_all(this->body_size) {
        std::cout << "body size: " << this->body_size << std::endl;
        std::cout << partials_pp_all.get().size_x() << " " << partials_pp_all.get().size_y() << std::endl;
    }

template<bool use_weighted_distribution> 
PartialHistogramManagerMT<use_weighted_distribution>::~PartialHistogramManagerMT() = default;

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> PartialHistogramManagerMT<use_weighted_distribution>::calculate() {
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
                calc_self_correlation(i);
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations
            else if (externally_modified[i]) {
                pool->detach_task(
                    [this, i] () {update_compact_representation_body(i);}
                );
            }
        }

        // merge the partial results from each thread and add it to the master histogram
        for (unsigned int i = 0; i < this->body_size; ++i) {
            if (internally_modified[i]) {
                pool->detach_task(
                    [this, i] () {combine_self_correlation(i);}
                );
            }
        }
    }

    // small efficiency improvement: if the hydration layer was modified, we can update the compact representations in parallel with the self-correlation
    if (hydration_modified) {
        pool->detach_task(
            [this] () {update_compact_representation_water();}
        );
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    // check if the hydration layer was modified
    if (hydration_modified) {
        calc_hh();
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            if (externally_modified[i] || externally_modified[j]) {
                // one of the bodies was modified, so we recalculate its partial histogram
                calc_pp(i, j);
            }
        }

        // we also have to remember to update the partial histograms with the hydration layer
        if (externally_modified[i] || hydration_modified) {
            calc_hp(i);
        }
    }

    // merge the partial results from each thread and add it to the master histogram
    {
        if (hydration_modified) {
            pool->detach_task(
                [this] () {combine_hh();}
            );
        }

        for (unsigned int i = 0; i < this->body_size; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                if (externally_modified[i] || externally_modified[j]) {
                    pool->detach_task(
                        [this, i, j] () {combine_pp(i, j);}
                    );
                }
            }

            if (externally_modified[i] || hydration_modified) {
                pool->detach_task(
                    [this, i] () {combine_hp(i);}
                );
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

    this->partials_hh = detail::PartialHistogram(axis);
    for (unsigned int i = 0; i < this->body_size; ++i) {
        this->partials_hp.index(i) = detail::PartialHistogram(axis);
        this->partials_pp.index(i, i) = detail::PartialHistogram(axis);
        
        calc_self_correlation(i);

        for (unsigned int j = 0; j < i; ++j) {
            this->partials_pp.index(i, j) = detail::PartialHistogram(axis);
        }
    }

    for (unsigned int i = 0; i < this->body_size; ++i) {
        pool->detach_task(
            [this, i] () {combine_self_correlation(i);}
        );
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_self_correlation(unsigned int index) {
    auto pool = utility::multi_threading::get_global_pool();
    update_compact_representation_body(index);

    // calculate internal distances between atoms
    static auto calc_internal = [] (
        container::ThreadLocalWrapper<container::Container2D<GenericDistribution1D_t>>& p_pp_all, 
        unsigned int index, 
        const detail::CompactCoordinates& coords, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_pp = p_pp_all.get().index(index, index);
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
    };

    // calculate self correlation
    static auto calc_self = [] (
        container::ThreadLocalWrapper<container::Container2D<GenericDistribution1D_t>>& p_pp_all,
        unsigned int index,
        const detail::CompactCoordinates& coords
    ) -> void {
        auto& p_pp = p_pp_all.get().index(index, index);
        p_pp.add(0, std::accumulate(coords.get_data().begin(), coords.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));
    };

    unsigned int atom_size = this->protein->atom_size();
    for (unsigned int i = 0; i < atom_size; i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, i, index, atom_size] () {calc_internal(this->partials_pp_all, index, this->coords_p[index], i, std::min(i+settings::general::detail::job_size, atom_size));}
        );
    }
    pool->detach_task(
        [this, index] () {calc_self(this->partials_pp_all, index, this->coords_p[index]);}
    );
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_pp(unsigned int n, unsigned int m) {
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_pp = [] (
        container::ThreadLocalWrapper<container::Container2D<GenericDistribution1D_t>>& p_pp_all,
        unsigned int n,
        unsigned int m,
        const detail::CompactCoordinates& coords_n, 
        const detail::CompactCoordinates& coords_m, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_pp = p_pp_all.get().index(n, m);
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
    };

    detail::CompactCoordinates& coords_n = this->coords_p[n];
    for (unsigned int i = 0; i < coords_n.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, n, m, i, &coords_n] () {calc_pp(this->partials_pp_all, n, m, coords_n, this->coords_p[m], i, std::min<int>(i+settings::general::detail::job_size, coords_n.size()));}
        );
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_hp(unsigned int index) {
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_hp = [] (
        container::ThreadLocalWrapper<container::Container1D<GenericDistribution1D_t>>& p_hp_all,
        unsigned int index,
        const detail::CompactCoordinates& coords_i, 
        const detail::CompactCoordinates& coords_h, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_hp = p_hp_all.get().index(index);
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
    };

    detail::CompactCoordinates& coords = this->coords_p[index];
    for (unsigned int i = 0; i < coords.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, index, i, &coords] () {
                calc_hp(this->partials_hp_all, index, coords, this->coords_h, i, std::min<int>(i+settings::general::detail::job_size, coords.size()));
            }
        );
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_hh() {
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_hh = [] (
        container::ThreadLocalWrapper<GenericDistribution1D_t>& p_hh_all,
        const detail::CompactCoordinates& coords_h, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_hh = p_hh_all.get();
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
    };

    // calculate self correlation
    static auto calc_self = [] (
        container::ThreadLocalWrapper<GenericDistribution1D_t>& p_hh_all,
        const detail::CompactCoordinates& coords_h
    ) -> void {
        auto& p_hh = p_hh_all.get();
        p_hh.add(0, std::accumulate(coords_h.get_data().begin(), coords_h.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));
    };

    for (unsigned int i = 0; i < this->coords_h.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, i] () {calc_hh(this->partials_hh_all, this->coords_h, i, std::min<int>(i+settings::general::detail::job_size, this->coords_h.size()));}
        );
    }
    pool->detach_task(
        [this] () {calc_self(this->partials_hh_all, this->coords_h);}
    );
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_self_correlation(unsigned int index) {
    std::cout << "combine self correlation(" << index << ")" << std::endl;
    GenericDistribution1D_t p_pp(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_pp_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(index, index).begin(), p_pp.begin(), std::plus<>());
    }
    auto v = p_pp.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_pp.index(index, index);
    this->partials_pp.index(index, index).get_counts() = std::move(v);
    this->master += this->partials_pp.index(index, index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_pp(unsigned int n, unsigned int m) {
    GenericDistribution1D_t p_pp(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_pp_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(n, m).begin(), p_pp.begin(), std::plus<>());
    }
    auto v = p_pp.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_pp.index(n, m);
    this->partials_pp.index(n, m).get_counts() = std::move(v);
    this->master += this->partials_pp.index(n, m);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_hp(unsigned int index) {
    GenericDistribution1D_t p_hp(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_hp_all.get_all()) {
        std::transform(p_hp.begin(), p_hp.end(), tmp.get().index(index).begin(), p_hp.begin(), std::plus<>());
    }
    auto v = p_hp.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_hp.index(index)*2;
    this->partials_hp.index(index).get_counts() = std::move(v);
    this->master += this->partials_hp.index(index)*2;
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_hh() {
    GenericDistribution1D_t p_hh(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_hh_all.get_all()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.get().begin(), p_hh.begin(), std::plus<>());
    }
    auto v = p_hh.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_hh;
    this->partials_hh.get_counts() = std::move(v);
    this->master += this->partials_hh;
    master_hist_mutex.unlock();
}

template class hist::PartialHistogramManagerMT<true>;
template class hist::PartialHistogramManagerMT<false>;