/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/symmetry/PartialSymmetryManagerMT.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/MultiThreading.h>
#include <container/ThreadLocalWrapper.h>

#include <numeric>

using namespace ausaxs;
using namespace ausaxs::hist;

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_body(unsigned int index) {
    this->coords_a[index] = detail::CompactCoordinates(this->protein->get_body(index).get_atoms());
    hist::detail::SimpleExvModel::apply_simple_excluded_volume(this->coords_a[index], this->protein);
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_water() {
    this->coords_w = detail::CompactCoordinates(this->protein->get_waters());
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_self_correlation(unsigned int index) {
    auto pool = utility::multi_threading::get_global_pool();
    update_compact_representation_body(index);
    for (auto& tmp : this->partials_aa_all.get_all()) {
        tmp.get().index(index, index) = GenericDistribution1D_t(this->master.axis.bins);
    }

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

    unsigned int atom_size = this->protein->size_atom();
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
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aa(unsigned int n, unsigned int m) {
    auto pool = utility::multi_threading::get_global_pool();
    for (auto& tmp : this->partials_aa_all.get_all()) {
        tmp.get().index(n, m) = GenericDistribution1D_t(this->master.axis.bins);
    }

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

    auto& coords_n = this->coords_a[n];
    for (unsigned int i = 0; i < coords_n.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, n, m, i, &coords_n] () {calc_pp(this->partials_aa_all, n, m, coords_n, this->coords_a[m], i, std::min<int>(i+settings::general::detail::job_size, coords_n.size()));}
        );
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aw(unsigned int index) {
    auto pool = utility::multi_threading::get_global_pool();
    for (auto& tmp : this->partials_aw_all.get_all()) {
        tmp.get().index(index) = GenericDistribution1D_t(this->master.axis.bins);
    }

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

    auto& coords = this->coords_a[index];
    for (unsigned int i = 0; i < coords.size(); i += settings::general::detail::job_size) {
        pool->detach_task(
            [this, index, i, &coords] () {
                calc_aw(this->partials_aw_all, index, coords, this->coords_w, i, std::min<int>(i+settings::general::detail::job_size, coords.size()));
            }
        );
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_ww() {
    auto pool = utility::multi_threading::get_global_pool();
    for (auto& tmp : this->partials_ww_all.get_all()) {
        tmp.get() = GenericDistribution1D_t(this->master.axis.bins);
    }

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
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_self_correlation(unsigned int index) {
    GenericDistribution1D_t p_pp(this->master.axis.bins);
    for (auto& tmp : this->partials_aa_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(index, index).begin(), p_pp.begin(), std::plus<>());
    }

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(index, index);
    this->partials_aa.index(index, index) = std::move(p_pp);
    this->master += this->partials_aa.index(index, index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa(unsigned int n, unsigned int m) {
    GenericDistribution1D_t p_pp(this->master.axis.bins);
    for (auto& tmp : this->partials_aa_all.get_all()) { // std::reference_wrapper<container::Container2D<GenericDistribution1D_t>>
        std::transform(p_pp.begin(), p_pp.end(), tmp.get().index(n, m).begin(), p_pp.begin(), std::plus<>());
    }

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(n, m);
    this->partials_aa.index(n, m) = std::move(p_pp);
    this->master += this->partials_aa.index(n, m);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aw(unsigned int index) {
    GenericDistribution1D_t p_aw(this->master.axis.bins);
    for (auto& tmp : this->partials_aw_all.get_all()) {
        std::transform(p_aw.begin(), p_aw.end(), tmp.get().index(index).begin(), p_aw.begin(), std::plus<>());
    }

    master_hist_mutex.lock();
    this->master -= 2*this->partials_aw.index(index);
    this->partials_aw.index(index) = std::move(p_aw);
    this->master += 2*this->partials_aw.index(index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_ww() {
    GenericDistribution1D_t p_hh(this->master.axis.bins);
    for (auto& tmp : this->partials_ww_all.get_all()) {
        std::transform(p_hh.begin(), p_hh.end(), tmp.get().begin(), p_hh.begin(), std::plus<>());
    }

    master_hist_mutex.lock();
    this->master -= this->partials_ww;
    this->partials_ww = std::move(p_hh);
    this->master += this->partials_ww;
    master_hist_mutex.unlock();
}

template class hist::PartialSymmetryManagerMT<true>;
template class hist::PartialSymmetryManagerMT<false>;