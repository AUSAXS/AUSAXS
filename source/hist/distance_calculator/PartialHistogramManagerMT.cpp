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
    : PartialHistogramManager<use_weighted_distribution>(protein) {}

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
            calc_aw(i);
        }
    }

    // merge the partial results from each thread and add it to the master histogram
    pool->wait(); // we have to wait for all calculations to finish before we can merge them
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
                    [this, i] () {combine_aw(i);}
                );
            }
        }
    }

    pool->wait();
    this->statemanager->reset();
    Distribution1D p = this->master.get_counts();
    return std::make_unique<DistanceHistogram>(std::move(p), this->master.get_axis());
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::update_compact_representation_body(unsigned int index) {
    this->coords_a[index] = detail::CompactCoordinates(this->protein->get_body(index));
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::update_compact_representation_water() {
    this->coords_w = detail::CompactCoordinates(this->protein->get_waters());
}

template<bool use_weighted_distribution> 
std::unique_ptr<ICompositeDistanceHistogram> PartialHistogramManagerMT<use_weighted_distribution>::calculate_all() {
    auto total = calculate();
    total->shorten_axis();
    unsigned int bins = total->get_axis().bins;

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_hh = this->partials_ww.get_counts();
    GenericDistribution1D_t p_pp = this->master.base.get_counts();
    GenericDistribution1D_t p_aw(bins, 0);
    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            detail::PartialHistogram& current = this->partials_aa.index(i, j);

            // iterate through each entry in the partial histogram
            std::transform(p_pp.begin(), p_pp.end(), current.get_counts().begin(), p_pp.begin(), std::plus<>());
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; ++i) {
        detail::PartialHistogram& current = this->partials_aw.index(i);

        // iterate through each entry in the partial histogram
        std::transform(p_aw.begin(), p_aw.end(), current.get_counts().begin(), p_aw.begin(), std::plus<>());
    }

    // p_aw is already resized
    p_hh.resize(bins);
    p_pp.resize(bins);

    return std::make_unique<CompositeDistanceHistogram>(
        std::move(p_pp), 
        std::move(p_aw), 
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
    this->partials_aa_all = container::Container2D<GenericDistribution1D_t>(this->body_size, this->body_size, axis.bins);
    this->partials_aw_all = container::Container1D<GenericDistribution1D_t>(this->body_size, axis.bins);
    this->partials_ww_all = GenericDistribution1D_t(axis.bins, 0);

    this->partials_ww = detail::PartialHistogram(axis);
    for (unsigned int i = 0; i < this->body_size; ++i) {
        this->partials_aw.index(i) = detail::PartialHistogram(axis);
        this->partials_aa.index(i, i) = detail::PartialHistogram(axis);
        
        calc_self_correlation(i);

        for (unsigned int j = 0; j < i; ++j) {
            this->partials_aa.index(i, j) = detail::PartialHistogram(axis);
        }
    }

    pool->wait();
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
            [this, i, index, atom_size] () {calc_internal(this->partials_aa_all, index, this->coords_a[index], i, std::min(i+settings::general::detail::job_size, atom_size));}
        );
    }
    pool->detach_task(
        [this, index] () {calc_self(this->partials_aa_all, index, this->coords_a[index]);}
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

    detail::CompactCoordinates& coords_n = this->coords_a[n];
    for (unsigned int i = 0; i < coords_n.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, n, m, i, &coords_n] () {calc_pp(this->partials_aa_all, n, m, coords_n, this->coords_a[m], i, std::min<int>(i+settings::general::detail::job_size, coords_n.size()));}
        );
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_aw(unsigned int index) {
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_aw = [] (
        container::ThreadLocalWrapper<container::Container1D<GenericDistribution1D_t>>& p_aw_all,
        unsigned int index,
        const detail::CompactCoordinates& coords_i, 
        const detail::CompactCoordinates& coords_w, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_aw = p_aw_all.get().index(index);
        for (unsigned int i = imin; i < imax; ++i) { // atom
            unsigned int j = 0;                      // water
            for (; j+7 < coords_w.size(); j+=8) {
                evaluate8<use_weighted_distribution, 1>(p_aw, coords_i, coords_w, i, j);
            }

            for (; j+3 < coords_w.size(); j+=4) {
                evaluate4<use_weighted_distribution, 1>(p_aw, coords_i, coords_w, i, j);
            }

            for (; j < coords_w.size(); ++j) {
                evaluate1<use_weighted_distribution, 1>(p_aw, coords_i, coords_w, i, j);
            }
        }
    };

    detail::CompactCoordinates& coords = this->coords_a[index];
    for (unsigned int i = 0; i < coords.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, index, i, &coords] () {
                calc_aw(this->partials_aw_all, index, coords, this->coords_w, i, std::min<int>(i+settings::general::detail::job_size, coords.size()));
            }
        );
    }
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::calc_hh() {
    auto pool = utility::multi_threading::get_global_pool();

    static auto calc_hh = [] (
        container::ThreadLocalWrapper<GenericDistribution1D_t>& p_hh_all,
        const detail::CompactCoordinates& coords_w, 
        unsigned int imin, 
        unsigned int imax
    ) -> void {
        auto& p_hh = p_hh_all.get();
        for (unsigned int i = imin; i < imax; ++i) {
            unsigned int j = i+1;
            for (; j+7 < coords_w.size(); j+=8) {
                evaluate8<use_weighted_distribution, 2>(p_hh, coords_w, coords_w, i, j);
            }

            for (; j+3 < coords_w.size(); j+=4) {
                evaluate4<use_weighted_distribution, 2>(p_hh, coords_w, coords_w, i, j);
            }

            for (; j < coords_w.size(); ++j) {
                evaluate1<use_weighted_distribution, 2>(p_hh, coords_w, coords_w, i, j);
            }
        }
    };

    // calculate self correlation
    static auto calc_self = [] (
        container::ThreadLocalWrapper<GenericDistribution1D_t>& p_hh_all,
        const detail::CompactCoordinates& coords_w
    ) -> void {
        auto& p_hh = p_hh_all.get();
        p_hh.add(0, std::accumulate(coords_w.get_data().begin(), coords_w.get_data().end(), 0.0, [](double sum, const hist::detail::CompactCoordinatesData& val) {return sum + val.value.w*val.value.w;}));
    };

    for (unsigned int i = 0; i < this->coords_w.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, i] () {calc_hh(this->partials_ww_all, this->coords_w, i, std::min<int>(i+settings::general::detail::job_size, this->coords_w.size()));}
        );
    }
    pool->detach_task(
        [this] () {calc_self(this->partials_ww_all, this->coords_w);}
    );
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_self_correlation(unsigned int index) {
    GenericDistribution1D_t p_pp(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_aa_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(index, index).begin(), p_pp.begin(), std::plus<>());
    }
    auto v = p_pp.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(index, index);
    this->partials_aa.index(index, index).get_counts() = std::move(v);
    this->master += this->partials_aa.index(index, index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_pp(unsigned int n, unsigned int m) {
    GenericDistribution1D_t p_pp(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_aa_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(n, m).begin(), p_pp.begin(), std::plus<>());
    }
    auto v = p_pp.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(n, m);
    this->partials_aa.index(n, m).get_counts() = std::move(v);
    this->master += this->partials_aa.index(n, m);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_aw(unsigned int index) {
    GenericDistribution1D_t p_aw(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_aw_all.get_all()) {
        std::transform(p_aw.begin(), p_aw.end(), tmp.get().index(index).begin(), p_aw.begin(), std::plus<>());
    }
    auto v = p_aw.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_aw.index(index)*2;
    this->partials_aw.index(index).get_counts() = std::move(v);
    this->master += this->partials_aw.index(index)*2;
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialHistogramManagerMT<use_weighted_distribution>::combine_hh() {
    GenericDistribution1D_t p_hh(this->master.get_axis().bins, 0);
    for (auto& tmp : this->partials_ww_all.get_all()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.get().begin(), p_hh.begin(), std::plus<>());
    }
    auto v = p_hh.as_vector();

    master_hist_mutex.lock();
    this->master -= this->partials_ww;
    this->partials_ww.get_counts() = std::move(v);
    this->master += this->partials_ww;
    master_hist_mutex.unlock();
}

template class hist::PartialHistogramManagerMT<true>;
template class hist::PartialHistogramManagerMT<false>;