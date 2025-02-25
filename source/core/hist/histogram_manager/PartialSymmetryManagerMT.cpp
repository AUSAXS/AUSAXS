/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include "constants/ConstantsAxes.h"
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

#include <list>
#include <functional>

/**
The indexing in this file is a bit tricky. 
The body symmetries (body.symmetry.get) contains the first actual symmetry at index 0
The coordinates are indexed such that the main body is at index 0, with the symmetries starting from index 1
Thus, we have a lot of +1s and -1s in the indexing to account for this.
**/

#define DEBUG_INFO_PSMMT true
#define DEBUG_INFO_PSMMT_EXTENDED false

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
int to_res_index(int body1, int symmetry1, int body2, int symmetry2) {
    auto f = [] (int body, int symmetry) {return (body+1)*100 + symmetry;};
    return f(body1, symmetry1)*100 + f(body2, symmetry2);
}

int to_res_index_self(int body, int symmetry) {
    return to_res_index(body, symmetry, body, symmetry);
}

int to_res_index_water(int body, int symmetry) {
    return to_res_index_self(body, symmetry) + water_res_index;
}

template<bool use_weighted_distribution>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMT<use_weighted_distribution>::calculate() {
    if (protein->size_water() == 0) {
        return _calculate<false>();
    } else {
        return _calculate<true>();
    }
}

template<bool use_weighted_distribution> template<bool hydration_enabled>
std::unique_ptr<DistanceHistogram> PartialSymmetryManagerMT<use_weighted_distribution>::_calculate() {
    auto externally_modified = this->statemanager->get_externally_modified_bodies();
    auto internally_modified = this->statemanager->get_internally_modified_bodies();
    auto symmetry_modified = this->statemanager->get_symmetry_modified_bodies();
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
        for (int ibody = 0; ibody < static_cast<int>(this->body_size); ++ibody) {

            // if the internal state was modified, we have to recalculate the self-correlation
            if (internally_modified[ibody]) {
                update_compact_representation_body(ibody); //? unnecessary to update whole body; enough to update main body
                calc_aa_self(calculator.get(), ibody);
            }

            // if the external state was modified, we have to update the coordinate representations for later calculations 
            // (implicitly done in calc_self_correlation)
            else if (externally_modified[ibody]) {
                pool->detach_task(
                    [this, ibody] () {update_compact_representation_body(ibody);}
                );
            }

            for (int isym = 0; isym < static_cast<int>(this->protein->get_body(ibody).size_symmetry()); ++isym) {
                if (symmetry_modified[ibody][isym]) {
                    pool->detach_task(
                        [this, ibody, isym] () {update_compact_representation_symmetry(ibody, isym+1);}
                    );
                }
            }
        }
    }

    if constexpr (hydration_enabled) {
        // small efficiency improvement: if the hydration layer was modified, 
        // we can update the compact representations in parallel with the self-correlation
        if (hydration_modified) {
            pool->detach_task(
                [this] () {update_compact_representation_water();}
            );
        }
    }
    pool->wait(); // ensure the compact representations have been updated before continuing

    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "atomic setup: " << std::endl;
        for (int ibody = 0; ibody < static_cast<int>(this->body_size); ++ibody) {
            for (int iatom = 0; iatom < static_cast<int>(this->protein->get_body(ibody).get_atoms().size()); ++iatom) {
                std::cout << "[" << ibody << 0 << 0 << iatom << "]: ";
                for (int i = 0; i < 3; ++i) {
                    std::cout << coords[ibody].atomic[0][0].get_data()[iatom].data[i] << " ";
                }
                std::cout << std::endl;
            }

            for (int isym = 0; isym < static_cast<int>(this->protein->get_body(ibody).size_symmetry()); ++isym) {
                for (int irepeat = 0; irepeat < static_cast<int>(this->protein->get_body(ibody).symmetry().get(isym).repeat); ++irepeat) {
                    for (int iatom = 0; iatom < static_cast<int>(this->protein->get_body(ibody).get_atoms().size()); ++iatom) {
                        std::cout << "[" << ibody << isym+1 << irepeat << iatom << "]: ";
                        for (int i = 0; i < 3; ++i) {
                            std::cout << coords[ibody].atomic[isym+1][irepeat].get_data()[iatom].data[i] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }
    #endif

    // prepare a list of tasks to be run after the calculations are done
    // this way we avoid having to maintain two identical but separate loops for calculations and combining
    using res_t = distance_calculator::SimpleCalculator<use_weighted_distribution>::run_result;
    std::list<std::function<void()>> combine_tasks;
    res_t res; // placeholder for future results

    auto enqueue_combine_ww = [this, pool, &res, &combine_tasks] () {
        combine_tasks.emplace_back(
            [this, pool, &res] () {
                #if DEBUG_INFO_PSMMT
                    std::cout << "accessing self index " << water_res_index << std::endl;
                #endif
                assert(res.self.contains(water_res_index) && "SymmetryManager::calculate: water result not found");
                pool->detach_task(
                    [this, r = std::move(res.self[water_res_index])] () mutable {combine_ww(std::move(r));}
                );
            }
        );
    };

    auto enqueue_combine_aw = [this, pool, &res, &combine_tasks] (int ibody, int isym) {
        combine_tasks.emplace_back(
            [this, pool, &res, ibody, isym] () {
                #if DEBUG_INFO_PSMMT
                    std::cout << "accessing cross index " << to_res_index_water(ibody, isym) << std::endl;
                #endif
                assert(res.cross.contains(to_res_index_water(ibody, isym)) && "SymmetryManager::calculate: aw cross result not found");
                pool->detach_task(
                    [this, ibody, isym, r = std::move(res.cross[to_res_index_water(ibody, isym)])] 
                    () mutable {combine_aw(ibody, isym, std::move(r));}
                );
            }
        );
    };

    auto enqueue_combine_aa = [this, pool, &res, &combine_tasks] (int ibody1, int isym1, int ibody2, int isym2) {
        combine_tasks.emplace_back(
            [this, pool, &res, ibody1, isym1, ibody2, isym2] () {
                #if DEBUG_INFO_PSMMT
                    std::cout << "accessing cross index " << to_res_index(ibody1, isym1, ibody2, isym2) << std::endl;
                #endif
                assert(res.cross.contains(to_res_index(ibody1, isym1, ibody2, isym2)) && "SymmetryManager::calculate: aa cross result not found");
                pool->detach_task(
                    [this, ibody1, isym1, ibody2, isym2, r = std::move(res.cross[to_res_index(ibody1, isym1, ibody2, isym2)])] 
                    () mutable {combine_aa(ibody1, isym1, ibody2, isym2, std::move(r));}
                );    
            }
        );
    };

    if constexpr (hydration_enabled) {
        // check if the hydration layer was modified
        if (hydration_modified) {
            calc_ww(calculator.get());
            enqueue_combine_ww();
        }
    }

    // iterate through the lower triangle and check if either of each pair of different bodies were modified
    for (int ibody1 = 0; ibody1 < static_cast<int>(this->body_size); ++ibody1) {
        // note: off-diagonal elements only
        for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
            if (externally_modified[ibody1] || externally_modified[ibody2]) {

                // external modification requires recalculation of all affected symmetries
                for (int isym1 = 0; isym1 < 1+static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
                    for (int isym2 = 0; isym2 < 1+static_cast<int>(this->protein->get_body(ibody2).size_symmetry()); ++isym2) {
                        calc_aa(calculator.get(), ibody1, isym1, ibody2, isym2);
                        enqueue_combine_aa(ibody1, isym1, ibody2, isym2);
                    }
                }
            }
        }

        if (!externally_modified[ibody1]) {

            // correlations between this main body and symmetries of other main bodies
            for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {
                if (externally_modified[ibody2]) {continue;}
                for (int isym2 = 0; isym2 < static_cast<int>(this->protein->get_body(ibody2).size_symmetry()); ++isym2) {
                    if (symmetry_modified[ibody2][isym2]) {
                        calc_aa(calculator.get(), ibody1, 0, ibody2, isym2+1);
                        enqueue_combine_aa(ibody1, 0, ibody2, isym2+1);
                    }
                }
            }

            // correlations between symmetries of this body and other bodies
            for (int isym1 = 0; isym1 < static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
                // cross-correlations with other bodies
                for (int ibody2 = 0; ibody2 < ibody1; ++ibody2) {

                    // cross-correlation with other main body
                    if (symmetry_modified[ibody1][isym1]) {
                        calc_aa(calculator.get(), ibody1, isym1+1, ibody2, 0);
                        enqueue_combine_aa(ibody1, isym1+1, ibody2, 0);
                    }
                    
                    // cross-correlations with symmetries in other main body
                    for (int isym2 = 0; isym2 < static_cast<int>(this->protein->get_body(ibody2).size_symmetry()); ++isym2) {
                        if (!(symmetry_modified[ibody1][isym1] || symmetry_modified[ibody2][isym2])) {continue;}
                        calc_aa(calculator.get(), ibody1, isym1+1, ibody2, isym2+1);
                        enqueue_combine_aa(ibody1, isym1+1, ibody2, isym2+1);
                    }
                }
            }
        }

        // diagonal elements only have to be recalculated if the symmetry was modified
        {
            for (int isym1 = 0; isym1 < static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
                if (!symmetry_modified[ibody1][isym1]) {continue;}

                // cross-correlation with main body
                calc_aa(calculator.get(), ibody1, isym1+1, ibody1, 0);
                enqueue_combine_aa(ibody1, isym1+1, ibody1, 0);

                // cross-correlations with other symmetries
                for (int isym2 = 0; isym2 < isym1; ++isym2) {
                    calc_aa(calculator.get(), ibody1, isym1+1, ibody1, isym2+1);
                    enqueue_combine_aa(ibody1, isym1+1, ibody1, isym2+1);
                }
            }
        }

        if constexpr (hydration_enabled) {
            // we also have to remember to update the partial histograms with the hydration layer
            if (externally_modified[ibody1] || hydration_modified) {

                // loop from 0 (main body) to 1+size_symmetry (last symmetry)
                for (int isym1 = 0; isym1 < 1+static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
                    calc_aw(calculator.get(), ibody1, isym1);
                    enqueue_combine_aw(ibody1, isym1);
                }
            } else {

                // hydration layer not modified, check for symmetry modifications
                for (int isym1 = 0; isym1 < static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
                    if (symmetry_modified[ibody1][isym1]) {
                        calc_aw(calculator.get(), ibody1, isym1+1);
                        enqueue_combine_aw(ibody1, isym1+1);
                    }
                }
            }
        }
    }

    // wait for all calculations to finish
    res = calculator->run();

    // start all queued combine tasks
    // these will update all partial histograms with the results of the calculations
    for (auto& task : combine_tasks) {
        task();
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
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_body(int ibody) {
    coords[ibody] = symmetry::detail::generate_transformed_data(this->protein->get_body(ibody));
    for (auto& c : coords[ibody].atomic) {
        for (auto& sym : c) {
            hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
        }
    }
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_symmetry(int ibody, int isym) {
    coords[ibody].atomic[isym] = symmetry::detail::generate_transformed_data(this->protein->get_body(ibody), isym-1).data;
    for (auto& sym : coords[ibody].atomic[isym]) {
        hist::detail::SimpleExvModel::apply_simple_excluded_volume(sym, this->protein);
    }
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::update_compact_representation_water() {
    coords_w = CompactCoordinates(this->protein->get_waters());
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
    for (int ibody1 = 0; ibody1 < static_cast<int>(this->body_size); ++ibody1) {
        for (int isym1 = 0; isym1 < 1+static_cast<int>(this->protein->get_body(ibody1).size_symmetry()); ++isym1) {
            for (int ibody2 = 0; ibody2 <= ibody1; ++ibody2) {
                for (int isym2 = 0; isym2 < 1+static_cast<int>(this->protein->get_body(ibody2).size_symmetry()); ++isym2) {
                    // iterate through each entry in the partial histogram
                    std::transform(p_aa.begin(), p_aa.end(), this->partials_aa.index(ibody1, ibody2).index(isym1, isym2).begin(), p_aa.begin(), std::plus<>());
                }
            }
            // iterate through each entry in the partial hydration-atom histograms
            std::transform(p_aw.begin(), p_aw.end(), this->partials_aw.index(ibody1).index(isym1).begin(), p_aw.begin(), std::plus<>());
        }
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
    for (int ibody = 0; ibody < static_cast<int>(this->body_size); ++ibody) {
        auto& body = this->protein->get_body(ibody);
        this->partials_aw.index(ibody) = SymmetryIndexer1D<detail::PartialHistogram<use_weighted_distribution>>(
            body.size_symmetry()+1,
            detail::PartialHistogram<use_weighted_distribution>(axis.bins)
        );

        // initialize body diagonal
        this->partials_aa.index(ibody, ibody) = 
            SymmetryIndexer2D<detail::PartialHistogram<use_weighted_distribution>>(
                body.size_symmetry()+1, 
                axis.bins
            )
        ;
        
        // initialize body cross indexes
        for (int ibody2 = 0; ibody2 < ibody; ++ibody2) {
            auto& body2 = this->protein->get_body(ibody2);
            this->partials_aa.index(ibody, ibody2) = 
                SymmetryIndexer2D<detail::PartialHistogram<use_weighted_distribution>>(
                    body.size_symmetry()+1, 
                    body2.size_symmetry()+1, 
                    axis.bins
                )
            ;
        }

        // calculate the self-correlation of each body
        update_compact_representation_body(ibody); //? unnecessary to update whole body; enough to update main body
        calc_aa_self(calculator, ibody);
    }

    auto res = calculator->run();
    for (int ibody = 0; ibody < static_cast<int>(body_size); ++ibody) {
        #if DEBUG_INFO_PSMMT
            std::cout << "accessing index " << to_res_index_self(ibody, 0) << std::endl;
        #endif
        assert(res.self.contains(to_res_index_self(ibody, 0)) && "The self-correlation result does not contain the expected index.");
        pool->detach_task(
            [this, ibody, r = std::move(res.self[to_res_index_self(ibody, 0)])] () mutable {combine_aa_self(ibody, std::move(r));}
        );
    }
}

template<bool use_weighted_distribution>
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aa_self(calculator_t calculator, int ibody) const {
    const auto& body = protein->get_body(ibody);
    #if DEBUG_INFO_PSMMT
        std::cout << "calc_aa_self[" << ibody << "]" << std::endl;
    #endif

    #if DEBUG_INFO_PSMMT
        std::cout << "\t[" << ibody << "][0]" << std::endl;
        std::cout << "\t\tstored at self index " << to_res_index_self(ibody, 0) << std::endl;
    #endif
    // calculate the self correlation within each body and symmetry, equal to (N_sym+1) * (main body self corr)
    calculator->enqueue_calculate_self(coords[ibody].atomic[0][0], 1+body.size_symmetry_total(), to_res_index_self(ibody, 0));
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_ww(calculator_t calculator) const {
    #if DEBUG_INFO_PSMMT
        std::cout << "calc_ww" << std::endl;
        std::cout << "\tstored at self index " << water_res_index << std::endl;
    #endif
    calculator->enqueue_calculate_self(coords_w, 1, water_res_index);
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aa(calculator_t calculator, int ibody1, int isym1, int ibody2, int isym2) const {
    const auto& body1 = protein->get_body(ibody1);
    const auto& body2 = protein->get_body(ibody2);
    int res_index = to_res_index(ibody1, isym1, ibody2, isym2);
    assert(ibody2 <= ibody1 && "ibody2 must be less than or equal to ibody1 to avoid double-counting");

    #if DEBUG_INFO_PSMMT
        std::cout << "calc_aa[" << ibody1 << isym1 << "][" << ibody2 << isym2 << "]" << std::endl;
    #endif

    // internal correlations within the same body
    if (ibody1 == ibody2) {
        assert(isym1 != isym2 && "This method is unsuitable for calculating self-correlations");

        // correlations between a symmetry and its host body
        if (isym2 == 0) {
            assert(isym1 < 1+static_cast<int>(body1.size_symmetry()) && "symmetry index out of bounds");
            const auto& sym1 = body1.symmetry().get(isym1-1);
            bool closed = sym1.is_closed();
            for (int irepeat1 = 0; irepeat1 < (sym1.repeat - closed); ++irepeat1) {
                const auto& body1_sym_atomic = coords[ibody1].atomic[isym1][irepeat1];
                int scale = sym1.repeat - irepeat1;
                if (irepeat1 == 0 && closed) {scale += 1;}
                calculator->enqueue_calculate_cross(coords[ibody2].atomic[0][0], body1_sym_atomic, scale, res_index);

                #if DEBUG_INFO_PSMMT
                    std::cout << "\t[" << ibody1 << isym1 << irepeat1 << "][" << ibody2 << isym2 << "0] x" << scale << std::endl;
                    std::cout << "\t\tstored at cross index " << res_index << std::endl;
                #endif
            }
            return;
        }
        assert(isym1 != 0 && "Attempting to calculate cross-correlations outside the lower triangle");
    }

    // symmetry 0 is the main body, so we have to treat it separately
    if (isym1 == 0 && isym2 == 0) {
        #if DEBUG_INFO_PSMMT
            std::cout << "\t[" << ibody1 << "00][" << ibody2 << "00]" << std::endl;
            std::cout << "\t\tstored at cross index " << res_index << std::endl;
        #endif
        calculator->enqueue_calculate_cross(coords[ibody1].atomic[0][0], coords[ibody2].atomic[0][0], 1, res_index);
        return;
    } else if (isym1 == 0) {
        assert(isym2 < 1+static_cast<int>(body2.size_symmetry()) && "symmetry index out of bounds");
        const auto& sym2 = body2.symmetry().get(isym2-1);
        for (int irepeat2 = 0; irepeat2 < sym2.repeat; ++irepeat2) {
            const auto& body2_sym_atomic = coords[ibody2].atomic[isym2][irepeat2];
            calculator->enqueue_calculate_cross(coords[ibody1].atomic[0][0], body2_sym_atomic, 1, res_index);

            #if DEBUG_INFO_PSMMT
                std::cout << "\t[" << ibody1 << "00][" << ibody2 << isym2 << irepeat2 << "]" << std::endl;
                std::cout << "\t\tstored at cross index " << res_index << std::endl;
            #endif
        }
        return;
    } else if (isym2 == 0) {
        assert(isym1 < 1+static_cast<int>(body1.size_symmetry()) && "symmetry index out of bounds");
        const auto& sym1 = body1.symmetry().get(isym1-1);
        for (int irepeat1 = 0; irepeat1 < sym1.repeat; ++irepeat1) {
            const auto& body1_sym_atomic = coords[ibody1].atomic[isym1][irepeat1];
            calculator->enqueue_calculate_cross(body1_sym_atomic, coords[ibody2].atomic[0][0], 1, res_index);

            #if DEBUG_INFO_PSMMT
                std::cout << "\t[" << ibody1 << isym1 << irepeat1 << "][" << ibody2 << "00]" << std::endl;
                std::cout << "\t\tstored at cross index " << res_index << std::endl;
            #endif
        }
        return;
    }

    // else iterate over the replications of both symmetries
    assert(isym1 < 1+static_cast<int>(body1.size_symmetry()) && "symmetry index out of bounds");
    assert(isym2 < 1+static_cast<int>(body2.size_symmetry()) && "symmetry index out of bounds");
    const auto& sym1 = body1.symmetry().get(isym1-1);
    const auto& sym2 = body2.symmetry().get(isym2-1);

    for (int irepeat1 = 0; irepeat1 < sym1.repeat; ++irepeat1) {
        const auto& body1_sym_atomic = coords[ibody1].atomic[isym1][irepeat1];
        for (int irepeat2 = 0; irepeat2 < sym2.repeat; ++irepeat2) {
            const auto& body2_sym_atomic = coords[ibody2].atomic[isym2][irepeat2];
            calculator->enqueue_calculate_cross(body1_sym_atomic, body2_sym_atomic, 1, res_index);

            #if DEBUG_INFO_PSMMT
                std::cout << "\t[" << ibody1 << isym1 << irepeat1 << "][" << ibody2 << isym2 << irepeat2 << "]" << std::endl;
                std::cout << "\t\tstored at cross index " << res_index << std::endl;
            #endif
        }
    }
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::calc_aw(calculator_t calculator, int ibody, int isym) const {
    const auto& body = protein->get_body(ibody);
    const auto& waters = coords_w;
    int res_index = to_res_index_water(ibody, isym);

    #if DEBUG_INFO_PSMMT
        std::cout << "calc_aw[" << ibody << isym << "]" << std::endl;
    #endif

    // symmetry 0 is the main body, so we have to treat it separately
    if (isym == 0) {
        #if DEBUG_INFO_PSMMT
            std::cout << "\t[" << ibody << "0]" << std::endl;
            std::cout << "\t\tstored at cross index " << res_index << std::endl;
        #endif
        calculator->enqueue_calculate_cross(coords[ibody].atomic[0][0], waters, 1, res_index);
        return;
    }

    // else iterate over its repititions
    assert(isym < 1+static_cast<int>(body.size_symmetry()) && "symmetry index out of bounds");
    const auto& sym = body.symmetry().get(isym-1);
    for (int irepeat = 0; irepeat < sym.repeat; ++irepeat) {
        const auto& body1_sym_atomic = coords[ibody].atomic[isym][irepeat];
        #if DEBUG_INFO_PSMMT
            std::cout << "\t[" << ibody << isym << irepeat << "]" << std::endl;
            std::cout << "\t\tstored at cross index " << res_index << std::endl;
        #endif
        calculator->enqueue_calculate_cross(body1_sym_atomic, waters, 1, res_index);
    }
}

#if DEBUG_INFO_PSMMT_EXTENDED
    #include <iomanip>
#endif
template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa_self(int ibody, int isym, GenericDistribution1D_t&& res) {
    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "combine_aa_self[" << ibody << isym << "]" << std::endl;
        std::cout << "\tremoving " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
        }
        std::cout << "\n\t\t" << std::flush;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << this->partials_aa.index(ibody, ibody).index(isym, isym).get_content(i) << " ";
        }
        std::cout << "\n\tadding " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << res.get_content(i) << " ";
        }
        std::cout << std::endl;
    #endif

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(ibody, ibody).index(isym, isym);
    this->partials_aa.index(ibody, ibody).index(isym, isym) = std::move(res);
    this->master += this->partials_aa.index(ibody, ibody).index(isym, isym);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa_self(int index, GenericDistribution1D_t&& res) {
    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "combine_aa_self[" << index << "]" << std::endl;
        std::cout << "\tremoving " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
        }
        std::cout << "\n\t\t" << std::flush;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< this->partials_aa.index(index, index).index(0, 0).get_content(i) << " ";
        }
        std::cout << "\n\tadding " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< res.get_content(i) << " ";
        }
        std::cout << std::endl;
    #endif

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(index, index).index(0, 0);
    this->partials_aa.index(index, index).index(0, 0) = std::move(res);
    this->master += this->partials_aa.index(index, index).index(0, 0);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aa(int ibody1, int isym1, int ibody2, int isym2, GenericDistribution1D_t&& res) {
    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "combine_aa[" << ibody1 << isym1 << ", " << ibody2 << isym2 << "]" << std::endl;
        std::cout << "\tremoving " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
        }
        std::cout << "\n\t\t" << std::flush;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< this->partials_aa.index(ibody1, ibody2).index(isym1, isym2).get_content(i) << " ";
        }
        std::cout << "\n\tadding " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< res.get_content(i) << " ";
        }
        std::cout << std::endl;
    #endif

    master_hist_mutex.lock();
    this->master -= this->partials_aa.index(ibody1, ibody2).index(isym1, isym2);
    this->partials_aa.index(ibody1, ibody2).index(isym1, isym2) = std::move(res);
    this->master += this->partials_aa.index(ibody1, ibody2).index(isym1, isym2);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_aw(int ibody, int isym, GenericDistribution1D_t&& res) {
    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "combine_aw[" << ibody << isym << "]" << std::endl;
        std::cout << "\tremoving " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
        }
        std::cout << "\n\t\t" << std::flush;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< this->partials_aw.index(ibody).index(isym).get_content(i) << " ";
        }
        std::cout << "\n\tadding " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< res.get_content(i) << " ";
        }
        std::cout << std::endl;
    #endif

    master_hist_mutex.lock();
    this->master -= this->partials_aw.index(ibody).index(isym);
    this->partials_aw.index(ibody).index(isym) = std::move(res);
    this->master += this->partials_aw.index(ibody).index(isym);
    master_hist_mutex.unlock();
}

template<bool use_weighted_distribution> 
void PartialSymmetryManagerMT<use_weighted_distribution>::combine_ww(GenericDistribution1D_t&& res) {
    #if DEBUG_INFO_PSMMT_EXTENDED
        std::cout << "combine_ww" << std::endl;
        std::cout << "\tremoving " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
        }
        std::cout << "\n\t\t" << std::flush;
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< this->partials_ww.get_content(i) << " ";
        }
        std::cout << "\n\tadding " << std::endl << "\t\t";
        for (int i = 0; i < 20; ++i) {
            std::cout << std::setw(4)<< this->partials_ww.get_content(i) << " ";
        }
        std::cout << std::endl;
    #endif

    master_hist_mutex.lock();
    this->master -= this->partials_ww;
    this->partials_ww = std::move(res);
    this->master += this->partials_ww;
    master_hist_mutex.unlock();
}

template std::unique_ptr<hist::DistanceHistogram> hist::PartialSymmetryManagerMT<true>::_calculate<true>();
template std::unique_ptr<hist::DistanceHistogram> hist::PartialSymmetryManagerMT<true>::_calculate<false>();
template std::unique_ptr<hist::DistanceHistogram> hist::PartialSymmetryManagerMT<false>::_calculate<true>();
template std::unique_ptr<hist::DistanceHistogram> hist::PartialSymmetryManagerMT<false>::_calculate<false>();
template class hist::PartialSymmetryManagerMT<true>;
template class hist::PartialSymmetryManagerMT<false>;