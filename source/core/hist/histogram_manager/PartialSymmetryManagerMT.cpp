/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <hist/histogram_manager/detail/SymmetryHelpers.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <hist/distance_calculator/SimpleCalculator.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/DistanceHistogram.h>
#include <hist/intensity_calculator/CompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <data/state/StateManager.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <utility/MultiThreading.h>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;

template<bool use_weighted_distribution> 
PartialSymmetryManagerMT<use_weighted_distribution>::PartialSymmetryManagerMT(observer_ptr<const data::Molecule> protein) 
    : IPartialHistogramManager(protein), 
      protein(protein),
      coords(this->body_size), 
      partials_aa(this->body_size, this->body_size), 
      partials_aw(this->body_size)
{}

template<bool use_weighted_distribution> 
PartialSymmetryManagerMT<use_weighted_distribution>::~PartialSymmetryManagerMT() = default;

int water_res_index = 1.31e8;
int to_res_index(int body, int symmetry) {
    return (body+1)*100 + symmetry;
}

int to_res_index(int body1, int symmetry1, int body2, int symmetry2) {
    return to_res_index(body1, symmetry1)*100 + to_res_index(body2, symmetry2);
}

template<bool use_weighted_distribution> 
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMT<use_weighted_distribution>::calculate() {
    std::vector<bool> externally_modified = this->statemanager->get_externally_modified_bodies();
    std::vector<bool> internally_modified = this->statemanager->get_internally_modified_bodies();
    bool hydration_modified = this->statemanager->is_modified_hydration();
    auto pool = utility::multi_threading::get_global_pool();
    auto calculator = std::make_unique<distance_calculator::SimpleCalculator<use_weighted_distribution>>();

    // check if the object has already been initialized
    if (this->master.empty()) [[unlikely]] {
        initialize(calculator.get()); 

        // since the initialization also calculates the self-correlation, mark it as unmodified to avoid desyncing its state
        internally_modified = std::vector<bool>(this->body_size, false);
    }

    // if not, we must first check if the atom coordinates have been changed in any of the bodies
    else {
        for (unsigned int i = 0; i < this->body_size; ++i) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[i]) {
                calc_aa_self(calculator.get(), i);
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations (implicitly done in calc_self_correlation)
            else if (externally_modified[i]) {
                pool->detach_task(
                    [this, i] () {update_compact_representation_body(i);}
                );
            }
        }
    }

    // small efficiency improvement: if the hydration layer was modified, we can update the compact representations in parallel with the self-correlation
    if (hydration_modified) {
        for (unsigned int i = 0; i < this->body_size; ++i) {
            pool->detach_task(
                [this, i] () {update_compact_representation_water(i);}
            );
        }
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    // check if the hydration layer was modified
    if (hydration_modified) {
        calc_ww(calculator.get());
    }

    // iterate through the lower triangle and check if either of each pair of bodies was modified
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            if (externally_modified[i] || externally_modified[j]) {
                // one of the bodies was modified, so we recalculate its partial histogram
                calc_aa(calculator.get(), i, j);
            }
        }

        // we also have to remember to update the partial histograms with the hydration layer
        if (externally_modified[i] || hydration_modified) {
            calc_aw(calculator.get(), i);
        }
    }

    // merge the partial results from each thread and add it to the master histogram
    // for this process, we first have to wait for all threads to finish
    // then we extract the results in the same order they were submitted to ensure correctness
    auto res = calculator->run();
    {
        if (hydration_modified) {
            std::cout << "searching for water result at index " << water_res_index << std::endl;
            assert(res.self.contains(water_res_index) && "SymmetryManager::calculate: water result not found");
            pool->detach_task(
                [this, r = std::move(res.self[water_res_index])] () mutable {combine_ww(std::move(r));}
            );
        }

        for (unsigned int i = 0; i < this->body_size; ++i) {
            if (internally_modified[i]) {
                std::cout << "searching for aa self " << i << " result at index " << to_res_index(i, 0) << std::endl;
                assert(res.self.contains(to_res_index(i, 0)) && "SymmetryManager::calculate: aa self result not found");
                pool->detach_task(
                    [this, i, r = std::move(res.self[to_res_index(i, 0)])] () mutable {combine_aa_self(i, std::move(r));}
                );
            }

            for (unsigned int j = 0; j < i; ++j) {
                if (externally_modified[i] || externally_modified[j]) {
                    std::cout << "searching for aa cross " << i << ", " << j << " result at index " << to_res_index(i, 0, j, 0) << std::endl;
                    assert(res.cross.contains(to_res_index(i, 0, j, 0)) && "SymmetryManager::calculate: aa cross result not found");
                    pool->detach_task(
                        [this, i, j, r = std::move(res.cross[to_res_index(i, 0, j, 0)])] () mutable {combine_aa(i, j, std::move(r));}
                    );
                }
            }

            if (externally_modified[i] || hydration_modified) {
                std::cout << "searching for aw cross " << i << " result at index " << to_res_index(i, 0) + water_res_index << std::endl;
                assert(res.cross.contains(to_res_index(i, 0) + water_res_index) && "SymmetryManager::calculate: aw cross result not found");
                pool->detach_task(
                    [this, i, r = std::move(res.cross[to_res_index(i, 0) + water_res_index])] () mutable {combine_aw(i, std::move(r));}
                );
            }
        }
    }

    this->statemanager->reset_to_false();
    pool->wait();

    // downsize our axes to only the relevant area
    GenericDistribution1D_t p_tot = this->master;
    int max_bin = 10; // minimum size is 10
    for (int i = (int) p_tot.size()-1; i >= 10; i--) {
        if (p_tot.index(i) != 0) {
            max_bin = i+1; // +1 since we usually use this for looping (i.e. i < max_bin)
            break;
        }
    }
    p_tot.resize(max_bin);

    return std::make_unique<DistanceHistogram>(std::move(p_tot));
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_body(int index) {
    coords[index] = symmetry::detail::generate_transformed_data(this->protein->get_body(index));
    for (auto& c : coords[index].atomic) {
        for (auto& sym : c) {
            hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
        }
    }
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_water(int index) {
    if (this->protein->get_body(index).size_water() == 0) {return;}
    coords[index].waters = CompactCoordinates(this->protein->get_body(index).get_waters());
}

template<bool use_weighted_distribution>
std::unique_ptr<ICompositeDistanceHistogram> PartialSymmetryManagerMT<use_weighted_distribution>::calculate_all() {
    auto total = calculate();
    int bins = total->get_total_counts().size();

    // determine p_tot
    GenericDistribution1D_t p_tot(bins);
    for (int i = 0; i < bins; ++i) {
        p_tot.index(i) = this->master.index(i);
    }

    // after calling calculate(), everything is already calculated, and we only have to extract the individual contributions
    GenericDistribution1D_t p_ww = this->partials_ww;
    GenericDistribution1D_t p_aa = this->master.base;
    GenericDistribution1D_t p_aw(bins);
    p_ww.resize(bins);
    p_aa.resize(bins);

    // iterate through all partial histograms in the upper triangle
    for (unsigned int i = 0; i < this->body_size; ++i) {
        for (unsigned int j = 0; j <= i; ++j) {
            // iterate through each entry in the partial histogram
            std::transform(p_aa.begin(), p_aa.end(), this->partials_aa.index(i, j).begin(), p_aa.begin(), std::plus<>());
        }
    }

    // iterate through all partial hydration-protein histograms
    for (unsigned int i = 0; i < this->body_size; ++i) {
        // iterate through each entry in the partial histogram
        std::transform(p_aw.begin(), p_aw.end(), this->partials_aw.index(i).begin(), p_aw.begin(), std::plus<>());
    }

    if constexpr (use_weighted_distribution) {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(Distribution1D(p_aa)), 
            std::move(Distribution1D(p_aw)), 
            std::move(Distribution1D(p_ww)), 
            std::move(p_tot)
        );
    } else {
        return std::make_unique<CompositeDistanceHistogram>(
            std::move(p_aa), 
            std::move(p_aw), 
            std::move(p_ww), 
            std::move(p_tot)
        );
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::initialize(calculator_t calculator) {
    auto pool = utility::multi_threading::get_global_pool();
    const Axis& axis = constants::axes::d_axis; 
    std::vector<double> p_base(axis.bins, 0);
    this->master = detail::MasterHistogram<use_weighted_distribution>(p_base, axis);
    this->partials_ww = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
    for (unsigned int i = 0; i < this->body_size; ++i) {
        this->partials_aw.index(i) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        this->partials_aa.index(i, i) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        calc_aa_self(calculator, i);

        for (unsigned int j = 0; j < i; ++j) {
            this->partials_aa.index(i, j) = detail::PartialHistogram<use_weighted_distribution>(axis.bins);
        }
    }

    auto res = calculator->run();
    assert(res.self.size() == this->body_size && "The number of self-correlation results does not match the number of bodies.");
    for (int i = 0; i < static_cast<int>(body_size); ++i) {
        assert(res.self.contains(to_res_index(i, 0)) && "The self-correlation result does not contain the expected index.");
        pool->detach_task(
            [this, i, r = std::move(res.self[to_res_index(i, 0)])] () mutable {combine_aa_self(i, std::move(r));}
        );
    }
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aa_self(calculator_t calculator, int index) {
    std::cout << "storing aa self " << index << " at index " << to_res_index(index, 0) << std::endl;
    update_compact_representation_body(index);
    const auto& body = protein->get_body(index);
    calculator->enqueue_calculate_self(coords[index].atomic[0][0], 1 + body.size_symmetry_total(), to_res_index(index, 0));
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_ww(calculator_t calculator) {
    std::cout << "storing ww at index " << water_res_index << std::endl;
    for (int i_body1 = 0; i_body1 < static_cast<int>(body_size); ++i_body1) {
        const auto& body1_waters = coords[i_body1].waters;
        calculator->enqueue_calculate_self(body1_waters, 1, water_res_index);
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aa(calculator_t calculator, int n, int m) {
    const auto& body1 = protein->get_body(n);
    const auto& body2 = protein->get_body(m);
    int res_index = to_res_index(n, 0, m, 0);
    std::cout << "storing aa cross " << n << ", " << m << " at index " << res_index << std::endl;

    // for every symmetry i_sym1 of body 1
    for (int i_sym1 = 0; i_sym1 < static_cast<int>(body1.size_symmetry()); ++i_sym1) {
        const auto& sym1 = body1.symmetry().get(i_sym1);

        // for every replication i_repeat1 in i_sym1
        for (int i_repeat1 = 0; i_repeat1 < sym1.repeat; ++i_repeat1) {
            const auto& body1_sym_atomic = coords[n].atomic[1+i_sym1][i_repeat1];

            // for every symmetry j_sym1 of body 2
            for (int j_sym1 = 0; j_sym1 < static_cast<int>(body2.size_symmetry()); ++j_sym1) {
                const auto& sym2 = body2.symmetry().get(j_sym1);

                // for every replication j_repeat1 in j_sym1
                for (int j_repeat1 = 0; j_repeat1 < sym2.repeat; ++j_repeat1) {
                    const auto& body2_sym_atomic = coords[m].atomic[1+j_sym1][j_repeat1];

                    // calculate the cross-correlation between the two replicates
                    calculator->enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, res_index);
                }
            }
        }
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aw(calculator_t calculator, int index) {
    std::cout << "storing aw " << index << " at index " << to_res_index(index, 0) + water_res_index << std::endl;
    const auto& body = protein->get_body(index);
    const auto& body1_atomic = coords[index].atomic[0][0];
    const auto& waters = coords[index].waters;
    int res_index = to_res_index(index, 0) + water_res_index;

    calculator->enqueue_calculate_cross(body1_atomic, waters, 1, res_index);
    for (int i_sym1 = 0; i_sym1 < static_cast<int>(body.size_symmetry()); ++i_sym1) {
        const auto& sym1 = body.symmetry().get(i_sym1);
        for (int i_repeat1 = 0; i_repeat1 < sym1.repeat; ++i_repeat1) {
            const auto& body1_sym_atomic = coords[index].atomic[1+i_sym1][i_repeat1];
            calculator->enqueue_calculate_cross(body1_sym_atomic, waters, 1, res_index);
        }
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa_self(int index, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(index, index);
    this->partials_aa.index(index, index) = std::move(res);
    this->master += this->partials_aa.index(index, index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa(int n, int m, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(n, m);
    this->partials_aa.index(n, m) = std::move(res);
    this->master += this->partials_aa.index(n, m);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aw(int index, GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_aw.index(index);
    this->partials_aw.index(index) = std::move(res);
    this->master += this->partials_aw.index(index);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_ww(GenericDistribution1D_t&& res) {
    master_hist_mutex.lock();
    this->master -= this->partials_ww;
    this->partials_ww = std::move(res);
    this->master += this->partials_ww;
    master_hist_mutex.unlock();
}

template class hist::PartialSymmetryManagerMT<true>;
template class hist::PartialSymmetryManagerMT<false>;