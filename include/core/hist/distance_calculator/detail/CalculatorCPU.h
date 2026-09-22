// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/CPUKernel.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <settings/GeneralSettings.h>
#include <settings/HistogramSettings.h>
#include <utility/MultiThreading.h>
#include <utility/observer_ptr.h>

#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

#define DEBUG_INFO false

namespace ausaxs::hist::distance_calculator::detail {
    /**
     * @brief Queues and evaluates pairwise distance histograms on the global thread pool.
     *
     * Each enqueue call immediately dispatches its work: the calculation is split into chunks that
     * are submitted to the thread pool right away, rather than being deferred until run(). run()
     * therefore does not start any work — it only blocks until all previously enqueued tasks have
     * finished, then merges the per-thread partial histograms and returns them. If the tasks have
     * already completed by the time run() is called, it does not block at all.
     *
     * Jobs that share a @c merge_id accumulate into the same result histogram, which saves memory
     * when many calculations contribute to a single histogram (e.g. symmetry copies). The optional
     * integer @c scaling factor multiplies a job's contribution; it is dispatched to a templated
     * implementation, so only a bounded set of values is supported (see dispatch_scaling).
     *
     * The caller must keep all submitted data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class CalculatorCPU {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        using ThreadLocalResult = container::ThreadLocalWrapper<GenericDistribution1D_t>;
        public:
            struct run_result {
                std::unordered_map<int, GenericDistribution1D_t> self;
                std::unordered_map<int, GenericDistribution1D_t> cross;
            };

            /**
             * @brief Construct a calculator whose result histograms span @a bin_count bins.
             */
            explicit CalculatorCPU(int bin_count) : bin_count(bin_count) {}

            /**
             * @brief Drain the thread pool before any of this object's state is released.
             *
             * Every enqueue call dispatches tasks that capture this object's result vectors, so those tasks
             * must not outlive it. run() normally waits for them, but an exception between enqueue and run
             * would otherwise unwind past the results while the pool is still writing to them.
             */
            ~CalculatorCPU() {
                auto* pool = utility::multi_threading::get_global_pool();
                pool->purge();
                pool->wait();
            }

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
            int enqueue_calculate_self(const hist::detail::CompactCoordinates<variable_bin_width>& a, int scaling = 1, int merge_id = -1) {
                auto [target, index] = resolve(self_results, self_merge_ids, merge_id);
                // every unordered pair is counted twice by this convention, but the diagonal only once
                dispatch_scaling(scaling, [&a, target] (auto s) {
                    enqueue_self<weighted_bins, variable_bin_width, 2*decltype(s)::value, decltype(s)::value>(a, target);
                });
                return index;
            }

            /**
             * @brief Queue a cross-correlation calculation.
             *
             * @param a1 The first set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param a2 The second set of data to calculate the cross-correlation for. The reference must be valid until calculate is called.
             * @param merge_id The result vector id this calculation can be merged into. Supplying this can save significant memory resources.
             * @param scaling The scaling factor to apply to the result.
             * @return The index of the data in the result vector.
             */
            int enqueue_calculate_cross(
                const hist::detail::CompactCoordinates<variable_bin_width>& a1,
                const hist::detail::CompactCoordinates<variable_bin_width>& a2,
                int scaling = 1, int merge_id = -1
            ) {
                auto [target, index] = resolve(cross_results, cross_merge_ids, merge_id);
                dispatch_scaling(scaling, [&a1, &a2, target] (auto s) {
                    enqueue_cross<variable_bin_width, 2*decltype(s)::value>(a1, a2, target);
                });
                return index;
            }

            /**
             * @brief Get the current size of the result vector.
             */
            int size_self_result() const {return static_cast<int>(self_results.size());}
            int size_cross_result() const {return static_cast<int>(cross_results.size());} //< @copydoc size_self_result

            /**
             * @brief Calculate the queued histograms.
             *        This will block until all calculations are done.
             *
             * @return The calculated histograms.
             */
            run_result run();

        private:
            /**
             * @brief The handle the queued tasks accumulate through; see Target.
             */
            struct Target {
                using entry_type = typename GenericDistribution1D_t::value_type;
                observer_ptr<ThreadLocalResult> results;
                std::span<entry_type> get() const {
                    auto& histogram = results->get();
                    return {&*histogram.begin(), static_cast<std::size_t>(histogram.size())};
                }
            };

            struct Resolved {
                Target target;
                int index; // position in the result vector, which is what callers index by
            };

            int bin_count;
            std::vector<std::unique_ptr<ThreadLocalResult>> self_results, cross_results;
            std::unordered_map<int, int> self_merge_ids, cross_merge_ids;

            /**
             * @brief Assign a merge id its result histogram, allocating one the first time it is seen.
             *        An id of -1 always allocates a fresh one.
             */
            Resolved resolve(
                std::vector<std::unique_ptr<ThreadLocalResult>>& results,
                std::unordered_map<int, int>& merge_ids,
                int merge_id
            ) {
                int res_idx;
                if (!merge_ids.contains(merge_id)) {
                    res_idx = static_cast<int>(results.size());
                    merge_id = merge_id == -1 ? res_idx : merge_id;
                    merge_ids[merge_id] = res_idx;
                    results.emplace_back(std::make_unique<ThreadLocalResult>(bin_count));
                } else {
                    res_idx = merge_ids[merge_id];
                    assert(results[res_idx]->get().size() == bin_count && "The result vector has the wrong size.");
                }
                return Resolved{.target=Target{results[res_idx].get()}, .index=res_idx};
            }
    };
}

template<bool weighted_bins, bool variable_bin_width>
inline typename ausaxs::hist::distance_calculator::detail::CalculatorCPU<weighted_bins, variable_bin_width>::run_result ausaxs::hist::distance_calculator::detail::CalculatorCPU<weighted_bins, variable_bin_width>::run() {
    auto* pool = utility::multi_threading::get_global_pool();
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
