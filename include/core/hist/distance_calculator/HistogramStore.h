// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <optional>
#include <span>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Temporary storage space for the pairwise distance histograms calculated by the CPU or GPU kernel.
     */
    template<bool weighted_bins>
    class HistogramStore {
        public:
            using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
            using GenericDistribution2D_t = typename hist::GenericDistribution2D<weighted_bins>::type;
            using GenericDistribution3D_t = typename hist::GenericDistribution3D<weighted_bins>::type;
            using entry_type = typename GenericDistribution1D_t::value_type;

            /**
             * @brief Where the tasks of one handle accumulate, as handed to the kernels in CPUKernel.h.
             *        Must be obtained from target(), on the thread that enqueues.
             */
            class Target {
                public:
                    using entry_type = typename HistogramStore::entry_type;

                    /**
                     * @brief The calling thread's copy of the row.
                     */
                    std::span<entry_type> get() const {return store->scratch[handle]->get();}

                private:
                    friend class HistogramStore;
                    Target(observer_ptr<HistogramStore> store, int handle) : store(store), handle(handle) {}

                    observer_ptr<HistogramStore> store;
                    int handle;
            };

            /**
             * @brief Construct a store of @a rows zeroed rows, each spanning @a bins bins.
             */
            HistogramStore(int rows, int bins)
                : n_rows(rows), n_bins(bins), primary(static_cast<std::size_t>(rows)*bins),
                  scratch(rows)
            {assert(0 <= rows && 0 <= bins && "HistogramStore: negative size.");}

            HistogramStore(const HistogramStore&) = delete;
            HistogramStore& operator=(const HistogramStore&) = delete;
            HistogramStore(HistogramStore&&) = default;
            HistogramStore& operator=(HistogramStore&&) = default;

            int rows() const {return n_rows;}
            int bins() const {return n_bins;}

            /**
             * @brief The bins of the row @a h.
             */
            std::span<entry_type> row(int h) {
                assert(0 <= h && h < n_rows && "HistogramStore::row: handle out of bounds.");
                return {primary.data() + static_cast<std::size_t>(h)*n_bins, static_cast<std::size_t>(n_bins)};
            }

            std::span<const entry_type> row(int h) const {
                assert(0 <= h && h < n_rows && "HistogramStore::row: handle out of bounds.");
                return {primary.data() + static_cast<std::size_t>(h)*n_bins, static_cast<std::size_t>(n_bins)};
            } //< @copydoc row(int)

            /**
             * @brief Zero the row @a h.
             */
            void clear(int h) {
                auto r = row(h);
                std::fill(r.begin(), r.end(), entry_type{});
            }

            /**
             * @brief The target a CPU calculation into @a h accumulates through.
             */
            Target target(int h) {
                assert(0 <= h && h < n_rows && "HistogramStore::target: handle out of bounds.");
                if (!scratch[h].has_value()) {
                    scratch[h].emplace(n_bins);
                    touched.push_back(h);
                }
                return Target(this, h);
            }

            /**
             * @brief Assign the sum of the per-thread copies of every handle given a target since the last call to its
             *        row, and free them. This must only be done after all calculations have finished. 
             */
            void fold() {
                for (int h : touched) {
                    auto r = row(h);
                    std::fill(r.begin(), r.end(), entry_type{});
                    for (const auto& local : scratch[h]->get_all()) {
                        std::transform(r.begin(), r.end(), local.get().begin(), r.begin(), std::plus<>());
                    }
                    scratch[h].reset();
                }
                touched.clear();
            }

            /**
             * @brief Get the result as a 1D distribution.
             */
            GenericDistribution1D_t export_1d(int h) const {
                auto r = row(h);
                return GenericDistribution1D_t(std::vector<entry_type>(r.begin(), r.end()));
            }

            /**
             * @brief Get the results as a 2D distribution. 
             */
            GenericDistribution2D_t export_2d(int first, int n) const {
                GenericDistribution2D_t result(n, n_bins);
                for (int i = 0; i < n; ++i) {
                    auto r = row(first + i);
                    std::copy(r.begin(), r.end(), result.begin(i));
                }
                return result;
            }

            /**
             * @brief Get the results as a 3D distribution. 
             */
            GenericDistribution3D_t export_3d(int first, int n1, int n2) const {
                GenericDistribution3D_t result(n1, n2, n_bins);
                for (int i = 0; i < n1; ++i) {
                    for (int j = 0; j < n2; ++j) {
                        auto r = row(first + i*n2 + j);
                        std::copy(r.begin(), r.end(), result.begin(i, j));
                    }
                }
                return result;
            }

        private:
            using Scratch = container::ThreadLocalWrapper<std::vector<entry_type>>;

            int n_rows = 0, n_bins = 0;
            std::vector<entry_type> primary; // the final, folded storage
            std::vector<std::optional<Scratch>> scratch; // one per row, lazily allocated
            std::vector<int> touched; // updated results since the last fold, required for the persistent (partial) managers
    };
}
