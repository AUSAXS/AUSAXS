// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/TemplateHelpers.h>
#include <container/ThreadLocalWrapper.h>
#include <utility/MultiThreading.h>
#include <settings/HistogramSettings.h>

#include <vector>
#include <unordered_map>

#define DEBUG_INFO false

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Simple interface to queue histogram calculations. 
     *        Submit the data and then call calculate to get the result.
     *        The caller must guarantee the lifetime of all submitted data.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SimpleCalculator {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        public:
            struct run_result {
                std::unordered_map<int, GenericDistribution1D_t> self;
                std::unordered_map<int, GenericDistribution1D_t> cross;
            };

            /**
             * @brief Queue a self-correlation calculation. 
             *        This is faster than calling the cross-correlation method with the same data, as some optimizations can be made. 
             *
             * @param a The data to calculate the self-correlation for. The reference must be valid until calculate is called.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources. 
             * @param scaling The scaling factor to apply to the result.
             *
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_self(const hist::detail::CompactCoordinates<variable_bin_width>& a, int scaling = 1, int merge_id = -1);

            /**
             * @brief Queue a cross-correlation calculation. 
             *
             * @param a1 The first set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param a2 The second set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources.
             * @param scaling The scaling factor to apply to the result.
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates<variable_bin_width>& a1, const hist::detail::CompactCoordinates<variable_bin_width>& a2, int scaling = 1, int merge_id = -1);

            /**
             * @brief Get the current size of the result vector. 
             */
            int size_self_result() const; 
            int size_cross_result() const; //< @copydoc size_self_result

            /**
             * @brief Calculate the queued histograms. 
             *        This will block until all calculations are done.
             *
             * @return The calculated histograms. 
             */
            run_result run();

        private:
            std::vector<std::unique_ptr<container::ThreadLocalWrapper<GenericDistribution1D_t>>> self_results, cross_results;
            std::unordered_map<int, int> self_merge_ids, cross_merge_ids;

            template<int scaling>
            int enqueue_calculate_self(const hist::detail::CompactCoordinates<variable_bin_width>& data, int merge_id);

            template<int scaling>
            int enqueue_calculate_cross(const hist::detail::CompactCoordinates<variable_bin_width>& data_1, const hist::detail::CompactCoordinates<variable_bin_width>& data_2, int merge_id);
    };
}

template<bool weighted_bins, bool variable_bin_width> template<int scaling>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::enqueue_calculate_self(
    const hist::detail::CompactCoordinates<variable_bin_width>& data, 
    int merge_id
) {
    auto pool = utility::multi_threading::get_global_pool();

    int res_idx;
    if (!self_merge_ids.contains(merge_id)) {
        res_idx = self_results.size();
        merge_id = merge_id == -1 ? res_idx : merge_id;
        self_merge_ids[merge_id] = res_idx;
        self_results.emplace_back(std::make_unique<container::ThreadLocalWrapper<GenericDistribution1D_t>>(settings::axes::bin_count));
    } else {
        res_idx = self_merge_ids[merge_id];
        assert(self_results[res_idx]->get().size() == settings::axes::bin_count && "The result vector has the wrong size.");
    }

    auto res_ptr = self_results[res_idx].get();
    int data_size = static_cast<int>(data.size());
    int job_size = settings::general::detail::job_size;

    // calculate upper triangle
    for (int i = 0; i < data_size; i+=job_size) {
        pool->detach_task(
            [&data, res_ptr, data_size, imin = i, imax = std::min(i+job_size, data_size)] () {
                auto& p_aa = res_ptr->get();
                for (int i = imin; i < imax; ++i) { // atom
                    int j = i+1;                    // atom
                    for (; j+7 < data_size; j+=8) {
                        evaluate8<weighted_bins, variable_bin_width, 2*scaling>(p_aa, data, data, i, j);
                    }

                    for (; j+3 < data_size; j+=4) {
                        evaluate4<weighted_bins, variable_bin_width, 2*scaling>(p_aa, data, data, i, j);
                    }

                    for (; j < data_size; ++j) {
                        evaluate1<weighted_bins, variable_bin_width, 2*scaling>(p_aa, data, data, i, j);
                    }
                }
            }
        );
    }

    // calculate skipped diagonal
    pool->detach_task(
        [&data, res_ptr] () {
            auto& p_aa = res_ptr->get();
            p_aa.add(0, scaling*std::accumulate(
                data.get_data().begin(), 
                data.get_data().end(), 
                0.0, 
                [] (double sum, const hist::detail::CompactCoordinatesXYZW<variable_bin_width>& val) {return sum + val.value.w*val.value.w;}
            ));
        }
    );

    return res_idx;
}

template<bool weighted_bins, bool variable_bin_width> template<int scaling>
int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::enqueue_calculate_cross(
    const hist::detail::CompactCoordinates<variable_bin_width>& data_1, 
    const hist::detail::CompactCoordinates<variable_bin_width>& data_2, 
    int merge_id
) {
    auto pool = utility::multi_threading::get_global_pool();

    int res_idx;
    if (!cross_merge_ids.contains(merge_id)) {
        res_idx = cross_results.size();
        merge_id = merge_id == -1 ? res_idx : merge_id;
        cross_merge_ids[merge_id] = res_idx;
        cross_results.emplace_back(std::make_unique<container::ThreadLocalWrapper<GenericDistribution1D_t>>(settings::axes::bin_count));
    } else {
        res_idx = cross_merge_ids[merge_id];
        assert(cross_results[res_idx]->get().size() == settings::axes::bin_count && "The result vector has the wrong size.");
    }

    auto res_ptr = cross_results[res_idx].get();
    int data_1_size = static_cast<int>(data_1.size());
    int data_2_size = static_cast<int>(data_2.size());
    int job_size = settings::general::detail::job_size;

    for (int i = 0; i < data_2_size; i+=job_size) {
        pool->detach_task(
            [&data_1, &data_2, res_ptr, data_1_size, imin = i, imax = std::min(i+job_size, data_2_size)] () {
                auto& p_ab = res_ptr->get();
                for (int i = imin; i < imax; ++i) { // b
                    int j = 0;                      // a
                    for (; j+7 < data_1_size; j+=8) {
                        evaluate8<weighted_bins, variable_bin_width, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j+3 < data_1_size; j+=4) {
                        evaluate4<weighted_bins, variable_bin_width, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }

                    for (; j < data_1_size; ++j) {
                        evaluate1<weighted_bins, variable_bin_width, 2*scaling>(p_ab, data_2, data_1, i, j);
                    }
                }
            }
        );
    }

    return res_idx;
}

template<bool weighted_bins, bool variable_bin_width>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::size_self_result() const {
    return self_results.size();
}

template<bool weighted_bins, bool variable_bin_width>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::size_cross_result() const {
    return cross_results.size();
}

template<bool weighted_bins, bool variable_bin_width>
inline typename ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::run_result ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::run() {
    auto pool = utility::multi_threading::get_global_pool();
    pool->wait();
    run_result result;

    #if DEBUG_INFO
        if (!self_merge_ids.empty()) {
            std::cout << "self results:" << std::endl;
            std::cout << "\t" << std::flush;
            for (int i = 0; i < 20; ++i) {
                std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
            }
            std::cout << std::endl;
        }
        for (auto[i, j] : self_merge_ids) {
            result.self[i] = self_results[j]->merge();
            std::cout << "\t";
            for (int k = 0; k < 20; ++k) {
                std::cout << std::setw(4) << result.self[i].get_content(k) << " ";
            }
            std::cout << std::endl;
        }
        if (!cross_merge_ids.empty()) {
            std::cout << "cross results:" << std::endl;
            std::cout << "\t" << std::flush;
            for (int i = 0; i < 20; ++i) {
                std::cout << std::setw(4) << constants::axes::d_vals[i] << " ";
            }
            std::cout << std::endl;
        }
        for (auto[i, j] : cross_merge_ids) {
            result.cross[i] = cross_results[j]->merge();
            std::cout << "\t";
            for (int k = 0; k < 20; ++k) {
                std::cout << std::setw(4) << result.cross[i].get_content(k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    #endif

    for (auto[i, j] : self_merge_ids) {
        result.self[i] = self_results[j]->merge();
    }

    for (auto[i, j] : cross_merge_ids) {
        result.cross[i] = cross_results[j]->merge();
    }

    // cleanup
    self_results.clear();
    cross_results.clear();
    self_merge_ids.clear();
    cross_merge_ids.clear();

    return result;
}

// this approach is not sustainable
// should higher scaling factors be needed, new add1, add4, and add8 functions should be created which accepts the scaling factor as a parameter
// for now, this is primarily meant for rigidbody optimizations, where larger symmetries are not expected
template<bool weighted_bins, bool variable_bin_width>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::enqueue_calculate_self(
    const hist::detail::CompactCoordinates<variable_bin_width>& data, 
    int scaling,
    int merge_id
) {
    switch (scaling) {
        case 1:  return enqueue_calculate_self<1>(data, merge_id);
        case 2:  return enqueue_calculate_self<2>(data, merge_id);
        case 3:  return enqueue_calculate_self<3>(data, merge_id);
        case 4:  return enqueue_calculate_self<4>(data, merge_id);
        case 5:  return enqueue_calculate_self<5>(data, merge_id);
        case 6:  return enqueue_calculate_self<6>(data, merge_id);
        case 7:  return enqueue_calculate_self<7>(data, merge_id);
        case 8:  return enqueue_calculate_self<8>(data, merge_id);
        case 9:  return enqueue_calculate_self<9>(data, merge_id);
        case 10: return enqueue_calculate_self<10>(data, merge_id);
        case 11: return enqueue_calculate_self<11>(data, merge_id);
        case 12: return enqueue_calculate_self<12>(data, merge_id);
        case 13: return enqueue_calculate_self<13>(data, merge_id);
        case 14: return enqueue_calculate_self<14>(data, merge_id);
        case 15: return enqueue_calculate_self<15>(data, merge_id);
        case 16: return enqueue_calculate_self<16>(data, merge_id);
        case 17: return enqueue_calculate_self<17>(data, merge_id);
        case 18: return enqueue_calculate_self<18>(data, merge_id);
        case 19: return enqueue_calculate_self<19>(data, merge_id);
        case 20: return enqueue_calculate_self<20>(data, merge_id);
        case 21: return enqueue_calculate_self<21>(data, merge_id);
        case 22: return enqueue_calculate_self<22>(data, merge_id);
        case 23: return enqueue_calculate_self<23>(data, merge_id);
        case 24: return enqueue_calculate_self<24>(data, merge_id);
        case 25: return enqueue_calculate_self<25>(data, merge_id);
        case 26: return enqueue_calculate_self<26>(data, merge_id);
        case 27: return enqueue_calculate_self<27>(data, merge_id);
        case 28: return enqueue_calculate_self<28>(data, merge_id);
        case 29: return enqueue_calculate_self<29>(data, merge_id);
        case 30: return enqueue_calculate_self<30>(data, merge_id);
        default: throw std::runtime_error("SimpleCalculator::enqueue_calculate_self: too large scaling factor (" + std::to_string(scaling) + ")");
    }
}

// this approach is not sustainable
// should higher scaling factors be needed, new add1, add4, and add8 functions should be created which accepts the scaling factor as a parameter
// for now, this is primarily meant for rigidbody optimizations, where larger symmetries are not expected
template<bool weighted_bins, bool variable_bin_width>
inline int ausaxs::hist::distance_calculator::SimpleCalculator<weighted_bins, variable_bin_width>::enqueue_calculate_cross(
    const hist::detail::CompactCoordinates<variable_bin_width>& data_1, 
    const hist::detail::CompactCoordinates<variable_bin_width>& data_2, 
    int scaling,
    int merge_id
) {
    switch (scaling) {
        case 1:  return enqueue_calculate_cross<1>(data_1, data_2, merge_id);
        case 2:  return enqueue_calculate_cross<2>(data_1, data_2, merge_id);
        case 3:  return enqueue_calculate_cross<3>(data_1, data_2, merge_id);
        case 4:  return enqueue_calculate_cross<4>(data_1, data_2, merge_id);
        case 5:  return enqueue_calculate_cross<5>(data_1, data_2, merge_id);
        case 6:  return enqueue_calculate_cross<6>(data_1, data_2, merge_id);
        case 7:  return enqueue_calculate_cross<7>(data_1, data_2, merge_id);
        case 8:  return enqueue_calculate_cross<8>(data_1, data_2, merge_id);
        case 9:  return enqueue_calculate_cross<9>(data_1, data_2, merge_id);
        case 10: return enqueue_calculate_cross<10>(data_1, data_2, merge_id);
        case 11: return enqueue_calculate_cross<11>(data_1, data_2, merge_id);
        case 12: return enqueue_calculate_cross<12>(data_1, data_2, merge_id);
        case 13: return enqueue_calculate_cross<13>(data_1, data_2, merge_id);
        case 14: return enqueue_calculate_cross<14>(data_1, data_2, merge_id);
        case 15: return enqueue_calculate_cross<15>(data_1, data_2, merge_id);
        case 16: return enqueue_calculate_cross<16>(data_1, data_2, merge_id);
        case 17: return enqueue_calculate_cross<17>(data_1, data_2, merge_id);
        case 18: return enqueue_calculate_cross<18>(data_1, data_2, merge_id);
        case 19: return enqueue_calculate_cross<19>(data_1, data_2, merge_id);
        case 20: return enqueue_calculate_cross<20>(data_1, data_2, merge_id);
        case 21: return enqueue_calculate_cross<21>(data_1, data_2, merge_id);
        case 22: return enqueue_calculate_cross<22>(data_1, data_2, merge_id);
        case 23: return enqueue_calculate_cross<23>(data_1, data_2, merge_id);
        case 24: return enqueue_calculate_cross<24>(data_1, data_2, merge_id);
        case 25: return enqueue_calculate_cross<25>(data_1, data_2, merge_id);
        case 26: return enqueue_calculate_cross<26>(data_1, data_2, merge_id);
        case 27: return enqueue_calculate_cross<27>(data_1, data_2, merge_id);
        case 28: return enqueue_calculate_cross<28>(data_1, data_2, merge_id);
        case 29: return enqueue_calculate_cross<29>(data_1, data_2, merge_id);
        case 30: return enqueue_calculate_cross<30>(data_1, data_2, merge_id);
        default: throw std::runtime_error("SimpleCalculator::enqueue_calculate_cross: too large scaling factor (" + std::to_string(scaling) + ")");
    }
}