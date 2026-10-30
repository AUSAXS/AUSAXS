// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <hist/distance_calculator/DistanceCalculatorFwd.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <utility/Exceptions.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
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
     *
     * The first bin of a histogram holds both the zero distance of each point with itself and the distinct pairs closer
     * than a bin width. The store records the former separately for every histogram, see get_self_1d() and get_self_3d().
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
             * @brief The result @a id of type @a T, i.e. one of GenericDistribution1D_t, GenericDistribution2D_t and GenericDistribution3D_t.
             */
            template<typename T>
            const T& get(int id) const {return std::get<T>(results[check_valid_id(id)]);}

            /**
             * @brief Call @a f with the result @a id, whichever shape it was allocated with.
             *        Like get_*(), this throws except::invalid_operation if the result was exported.
             */
            template<typename F>
            decltype(auto) visit(int id, F&& f) const {
                return std::visit(
                    [&f, id] <typename T> (const T& result) {
                        if constexpr (std::is_same_v<T, std::monostate>) {HistogramStore::throw_exported(id);}
                        else {return f(result);}
                    },
                    results[check_valid_id(id)]
                );
            }

            /**
             * @brief The self-correlation weight of the result @a id of allocate_1d(): how much of its first bin is the zero distance
             *        of points with themselves, i.e. their summed squared weights. It is zero unless a self-correlation was calculated into it.
             *        It must be read before the result is exported.
             */
            double get_self_1d(int id) const {
                const auto& result = get_1d(id);
                return get_self({result.begin(), result.end()});
            }

            /**
             * @brief The self-correlation weight of each class with itself of the result @a id of allocate_3d(). See get_self_1d().
             */
            std::vector<double> get_self_3d(int id) const {
                const auto& result = get_3d(id);
                std::vector<double> self(n_classes);
                for (int i = 0; i < n_classes; ++i) {self[i] = get_self(result.row(i, i));}
                return self;
            }

            /**
             * @brief Move the result @a id out of the store. Any later use of @a id throws except::invalid_operation.
             */
            GenericDistribution1D_t export_1d(int id) {return take<GenericDistribution1D_t>(id);}
            GenericDistribution2D_t export_2d(int id) {return take<GenericDistribution2D_t>(id);} //< @copydoc export_1d
            GenericDistribution3D_t export_3d(int id) {return take<GenericDistribution3D_t>(id);} //< @copydoc export_1d

        private:
            template<bool, bool> friend class Calculator;
            template<bool, bool> friend class detail::CalculatorCPU;
            template<bool, bool> friend class detail::GPUKernel;
            using Scratch = container::ThreadLocalWrapper<std::vector<entry_type>>;

            /**
             * @brief Where the tasks of one row accumulate at one factor, as handed to the kernels in CPUKernel.h.
             *        Everything accumulated through it is multiplied by that factor when it is folded into the row.
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
            std::vector<std::variant<std::monostate, GenericDistribution1D_t, GenericDistribution2D_t, GenericDistribution3D_t>> results; // monostate once exported
            std::map<std::pair<entry_type*, int>, Scratch> scratch; // per (row, factor) with queued calculations, keyed by the first bin of the row
            std::vector<std::span<entry_type>> resets;        // rows calculated from an empty set, zeroed by the next fold
            std::map<const entry_type*, double> scratch_diagonal; // per row, the self-correlation weight added since the last fold, see add_diagonal()
            std::map<const entry_type*, double> self_weights; // per row, the self-correlation weight its first bin holds, see assign_diagonal()

            /**
             * @brief The histogram of the result @a id of allocate_1d(), the one of class @a i of the result @a id of
             *        allocate_2d(), and the one of the class pair (@a i, @a j) of the result @a id of allocate_3d().
             */
            std::span<entry_type> row(int id) {
                auto& result = std::get<GenericDistribution1D_t>(results[check_valid_id(id)]);
                return {result.begin(), result.end()};
            }
            std::span<entry_type> row(int id, int i) {return std::get<GenericDistribution2D_t>(results[check_valid_id(id)]).row(i);}           //< @copydoc row(int)
            std::span<entry_type> row(int id, int i, int j) {return std::get<GenericDistribution3D_t>(results[check_valid_id(id)]).row(i, j);} //< @copydoc row(int)

            /**
             * @brief The target a CPU calculation into @a row accumulates through, scaled by @a factor on the next fold.
             *        Calculations at different factors into the same row get separate targets, so each can be scaled on its own.
             *        Must be called on the thread that enqueues.
             */
            Target target(std::span<entry_type> row, int factor) {
                assert(0 < factor && "HistogramStore: the factor must be positive.");
                auto [it, _] = scratch.try_emplace(std::pair(row.data(), factor), n_bins);
                return Target(&it->second);
            }

            /**
             * @brief Add the zero distance of the points of a self-correlation with themselves, of total @a weight, to @a row on the
             *        next fold. Only rows with a target are folded, so the distinct pairs of the same calculation must be given a
             *        target into @a row. Must be called on the thread that enqueues.
             */
            void add_diagonal(std::span<entry_type> row, double weight) {scratch_diagonal[row.data()] += weight;}

            /**
             * @brief Add the zero distance of the points of @a row with themselves, of total @a weight, to its first bin, and record it
             *        as the self-correlation weight of the row. The record is replaced rather than added to, since a row is always
             *        calculated in full.
             */
            void assign_diagonal(std::span<entry_type> row, double weight) {
                if constexpr (weighted_bins) {row[0] += hist::detail::WeightedEntry(weight, static_cast<std::int64_t>(weight), 0);}
                else {row[0] += weight;}
                self_weights[row.data()] = weight;
            }

            /**
             * @brief Zero @a row on the next fold, since the calculation queued into it has nothing to calculate.
             */
            void reset(std::span<entry_type> row) {resets.push_back(row);}

            /**
             * @brief Zero the rows that were reset, and assign every row given a target the sum over its targets of the
             *        factor times the sum of its per-thread copies, plus the diagonal added to it.
             *        This must only be done after all calculations have finished.
             */
            void fold() {
                for (auto row : resets) {
                    std::fill(row.begin(), row.end(), entry_type{});
                    self_weights.erase(row.data());
                }
                resets.clear();
                entry_type* previous = nullptr;
                for (auto& [key, local] : scratch) {
                    auto [start, factor] = key;
                    std::span<entry_type> row(start, static_cast<std::size_t>(n_bins));

                    // the targets of a row are adjacent, since the map is ordered by row first
                    if (start != previous) {
                        std::fill(row.begin(), row.end(), entry_type{});
                        auto diagonal = scratch_diagonal.extract(start);
                        assign_diagonal(row, diagonal.empty() ? 0 : diagonal.mapped());
                    }
                    previous = start;

                    for (const auto& copy : local.get_all()) {
                        std::transform(
                            row.begin(), row.end(), copy.get().begin(), row.begin(),
                            [factor] (const entry_type& sum, const entry_type& bin) {return sum + bin*factor;}
                        );
                    }
                }
                assert(scratch_diagonal.empty() && "HistogramStore: a diagonal was added to a row without a target.");
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

                // the rows of a result are contiguous, so their self-correlation weights are adjacent in the map
                const entry_type* first = std::to_address(result.begin());
                const entry_type* last = std::to_address(result.end());
                self_weights.erase(self_weights.lower_bound(first), self_weights.lower_bound(last));

                T taken = std::move(result);
                results[id] = std::monostate{};
                return taken;
            }

            // @a id, after checking that it is a result of this store which has not been exported
            int check_valid_id(int id) const {
                assert(0 <= id && id < static_cast<int>(results.size()) && "HistogramStore: unknown result id.");
                if (std::holds_alternative<std::monostate>(results[id])) [[unlikely]] {throw_exported(id);}
                return id;
            }

            [[noreturn]] static void throw_exported(int id) {
                throw except::invalid_operation("HistogramStore: result " + std::to_string(id) + " was exported and can no longer be used.");
            }

            // the self-correlation weight held by the first bin of @a row, see assign_diagonal()
            double get_self(std::span<const entry_type> row) const {
                auto it = self_weights.find(row.data());
                return it == self_weights.end() ? 0 : it->second;
            }
    };
}
