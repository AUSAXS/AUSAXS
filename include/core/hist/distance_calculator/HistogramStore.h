// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <span>
#include <unordered_map>
#include <variant>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Storage for the results of the pairwise distance histograms queued on a Calculator.
     *
     * The caller allocates one result per quantity it needs, before queueing anything into it. The shape of a result
     * follows from the coordinate sets it is calculated from, where a set is either flat or partitioned into classes():
     *   - allocate_1d(): a single histogram, for a flat set with itself or with another flat set.
     *   - allocate_2d(): one histogram per class, for a partitioned set with a flat set.
     *   - allocate_3d(): one histogram per unordered class pair (Shape::Triangular), for a partitioned set with itself or with another partitioned set.
     *
     * The store owns the results. They are only written when the calculator runs, so until then they hold what the
     * previous run left, and they are zero before the first. Read them with get_*() or move them out with export_*().
     */
    template<bool weighted_bins>
    class HistogramStore {
        public:
            using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
            using GenericDistribution2D_t = typename hist::GenericDistribution2D<weighted_bins>::type;
            using GenericDistribution3D_t = typename hist::GenericDistribution3D<weighted_bins, hist::Shape::Triangular>::type;
            using entry_type = typename GenericDistribution1D_t::value_type;

            /**
             * @brief Construct an empty store of histograms spanning @a bins bins.
             * @param classes The number of classes the partitioned coordinate sets are split into.
             */
            explicit HistogramStore(int bins, int classes = 1) : n_bins(bins), n_classes(classes) {
                assert(0 <= bins && 0 < classes && "HistogramStore: invalid size.");
            }

            HistogramStore(const HistogramStore&) = delete;
            HistogramStore& operator=(const HistogramStore&) = delete;
            HistogramStore(HistogramStore&&) = default;
            HistogramStore& operator=(HistogramStore&&) = default;

            int bins() const {return n_bins;}
            int classes() const {return n_classes;}

            /**
             * @brief Allocate a zeroed result of a single histogram, and return its id.
             */
            int allocate_1d() {return allocate(GenericDistribution1D_t(n_bins));}

            /**
             * @brief Allocate a zeroed result of one histogram per class, and return its id.
             */
            int allocate_2d() {return allocate(GenericDistribution2D_t(n_classes, n_bins));}

            /**
             * @brief Allocate a zeroed result of one histogram per class pair, and return its id.
             */
            int allocate_3d() {return allocate(GenericDistribution3D_t(n_classes, n_classes, n_bins));}

            /**
             * @brief The result @a id, as the last run of the calculator left it.
             */
            const GenericDistribution1D_t& get_1d(int id) const {return std::get<GenericDistribution1D_t>(results[check_valid_id(id)]);}
            const GenericDistribution2D_t& get_2d(int id) const {return std::get<GenericDistribution2D_t>(results[check_valid_id(id)]);} //< @copydoc get_1d
            const GenericDistribution3D_t& get_3d(int id) const {return std::get<GenericDistribution3D_t>(results[check_valid_id(id)]);} //< @copydoc get_1d

            /**
             * @brief Move the result @a id out of the store. It must not be used again afterwards.
             */
            GenericDistribution1D_t export_1d(int id) {return take<GenericDistribution1D_t>(id);}
            GenericDistribution2D_t export_2d(int id) {return take<GenericDistribution2D_t>(id);} //< @copydoc export_1d
            GenericDistribution3D_t export_3d(int id) {return take<GenericDistribution3D_t>(id);} //< @copydoc export_1d

        private:
            template<bool, bool, bool> friend class Calculator;
            template<bool, bool, bool> friend class detail::CalculatorCPU;
            template<bool, bool, bool> friend class detail::GPUKernel;
            using Scratch = container::ThreadLocalWrapper<std::vector<entry_type>>;

            /**
             * @brief Where the tasks of one row accumulate, as handed to the kernels in CPUKernel.h.
             */
            class Target {
                public:
                    using entry_type = typename HistogramStore::entry_type;

                    /**
                     * @brief The calling thread's copy of the row.
                     */
                    std::span<entry_type> get() const {return scratch->get();}

                private:
                    friend class HistogramStore;
                    explicit Target(observer_ptr<Scratch> scratch) : scratch(scratch) {}

                    observer_ptr<Scratch> scratch;
            };

            int n_bins, n_classes;
            std::vector<std::variant<GenericDistribution1D_t, GenericDistribution2D_t, GenericDistribution3D_t>> results;
            std::unordered_map<entry_type*, Scratch> scratch; // per row with queued calculations, keyed by its first bin
            std::vector<std::span<entry_type>> resets;        // rows calculated from an empty set, zeroed by the next fold

            /**
             * @brief The histogram of the result @a id of allocate_1d(), the one of class @a i of the result @a id of
             *        allocate_2d(), and the one of the class pair (@a i, @a j) of the result @a id of allocate_3d().
             */
            std::span<entry_type> row(int id) {
                auto& result = std::get<GenericDistribution1D_t>(results[check_valid_id(id)]);
                return check_live_result({result.begin(), result.end()});
            }
            std::span<entry_type> row(int id, int i) {return check_live_result(std::get<GenericDistribution2D_t>(results[check_valid_id(id)]).row(i));}           //< @copydoc row(int)
            std::span<entry_type> row(int id, int i, int j) {return check_live_result(std::get<GenericDistribution3D_t>(results[check_valid_id(id)]).row(i, j));} //< @copydoc row(int)

            /**
             * @brief The target a CPU calculation into @a row accumulates through. Must be called on the thread that enqueues.
             */
            Target target(std::span<entry_type> row) {
                auto [it, _] = scratch.try_emplace(row.data(), n_bins);
                return Target(&it->second);
            }

            /**
             * @brief Zero @a row on the next fold, since the calculation queued into it has nothing to calculate.
             */
            void reset(std::span<entry_type> row) {resets.push_back(row);}

            /**
             * @brief Zero the rows that were reset, and assign every row given a target the sum of its per-thread copies.
             *        This must only be done after all calculations have finished.
             */
            void fold() {
                for (auto row : resets) {std::fill(row.begin(), row.end(), entry_type{});}
                resets.clear();
                for (auto& [start, local] : scratch) {
                    std::span<entry_type> row(start, static_cast<std::size_t>(n_bins));
                    std::fill(row.begin(), row.end(), entry_type{});
                    for (const auto& copy : local.get_all()) {
                        std::transform(row.begin(), row.end(), copy.get().begin(), row.begin(), std::plus<>());
                    }
                }
                scratch.clear();
            }

            template<typename T>
            int allocate(T&& result) {
                // a queued calculation holds a pointer into its result, which must not be moved while it is queued
                assert(scratch.empty() && resets.empty() && "HistogramStore: results must be allocated before anything is queued.");
                results.emplace_back(std::forward<T>(result));
                return static_cast<int>(results.size())-1;
            }

            template<typename T>
            T take(int id) {
                assert(scratch.empty() && resets.empty() && "HistogramStore: the calculator must run before its results are exported.");
                auto& result = std::get<T>(results[check_valid_id(id)]);
                T taken = std::move(result);
                result = T{};
                return taken;
            }

            int check_valid_id(int id) const {
                assert(0 <= id && id < static_cast<int>(results.size()) && "HistogramStore: unknown result id.");
                return id;
            }

            // an exported result is left empty, and fold() would write n_bins entries past its end
            std::span<entry_type> check_live_result(std::span<entry_type> row) const {
                assert(static_cast<int>(row.size()) == n_bins && "HistogramStore: calculation queued into an exported result.");
                return row;
            }
    };
}
