// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/distance_calculator/detail/CPUKernel.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <utility/MultiThreading.h>

#include <cassert>
#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Evaluates the form factor-resolved pairwise distance histograms on the global thread pool.
     *
     * The form factor of an atom only decides *which* histogram its pairs land in, never the distance itself, so a
     * calculation restricted to two form factor types is an ordinary pairwise distance calculation whose result
     * happens to be one row of a larger distribution. This calculator therefore splits each input by form factor and
     * hands the resulting subsets to the same detail::enqueue_self and detail::enqueue_cross the weight-based calculator
     * uses, with a target naming the row to accumulate into.
     *
     * The rows are rows of one distribution per thread, not histograms of their own, so a worker thread keeps
     * accumulating into the single array it owns and everything is merged once, at the end.
     *
     * Contributions are counts: the subsets carry unit weights, and the form factor amplitudes are applied later,
     * by the intensity calculator. The caller must keep all submitted data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class CalculatorFF {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        using GenericDistribution2D_t = typename hist::GenericDistribution2D<weighted_bins>::type;
        using GenericDistribution3D_t = typename hist::GenericDistribution3D<weighted_bins>::type;
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using CompactCoordinatesFF_t = hist::detail::CompactCoordinatesFF<variable_bin_width>;
        public:
            struct run_result {
                GenericDistribution3D_t aa; // ff_type1, ff_type2, distance
                GenericDistribution2D_t aw; // ff_type, distance
                GenericDistribution1D_t ww; // distance
            };

            /**
             * @brief Construct a calculator whose result histograms span @a bin_count bins.
             */
            explicit CalculatorFF(int bin_count)
                : bin_count(bin_count), n_ff(form_factor::get_active_count()),
                  aa(n_ff, n_ff, bin_count), aw(n_ff, bin_count), ww(bin_count) {}

            /**
             * @brief Queue the self-correlation of @a a, resolved by form factor on both sides, into the aa result.
             *
             * Every unordered pair is counted twice, once in each of the two form factor orders it could be read in;
             * a pair of unlike types lands wholly in one of them rather than being split between the two. Only the
             * sum over both form factor indices is ever read, so this is the same distribution.
             */
            void enqueue_self_by_ff(const CompactCoordinatesFF_t& a) {
                const auto& parts = partition(a);
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    if (parts[ff1].empty()) {continue;}
                    detail::enqueue_self<weighted_bins, variable_bin_width, 2, 1>(parts[ff1], target_aa(ff1, ff1));
                    for (int ff2 = ff1+1; ff2 < n_ff; ++ff2) {
                        detail::enqueue_balanced_cross<variable_bin_width, 2>(parts[ff1], parts[ff2], target_aa(ff1, ff2));
                    }
                }
            }

            /**
             * @brief Queue the cross-correlation of @a a and @a b, resolved by form factor on the @a a side only,
             *        into the aw result. Each pair is counted once, as that convention expects.
             */
            void enqueue_cross_by_ff(const CompactCoordinatesFF_t& a, const CompactCoordinatesFF_t& b) {
                const auto& parts = partition(a);
                const auto& whole = flatten(b);
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    detail::enqueue_balanced_cross<variable_bin_width, 1>(parts[ff1], whole, target_aw(ff1));
                }
            }

            /**
             * @brief Queue the self-correlation of @a b, ignoring form factors, into the ww result.
             *        Every unordered pair is counted twice, as that convention expects.
             */
            void enqueue_self_flat(const CompactCoordinatesFF_t& b) {
                const auto& whole = flatten(b);
                if (whole.empty()) {return;}
                detail::enqueue_self<weighted_bins, variable_bin_width, 2, 1>(whole, target_ww());
            }

            /**
             * @brief Calculate the queued histograms.
             *        This will block until all calculations are done.
             */
            run_result run() {
                auto* pool = utility::multi_threading::get_global_pool();
                pool->wait();

                run_result result{.aa=aa.merge(), .aw=aw.merge(), .ww=ww.merge()};

                // cleanup
                aa.reinitialize_all(n_ff, n_ff, bin_count);
                aw.reinitialize_all(n_ff, bin_count);
                ww.reinitialize_all(bin_count);
                partitions.clear();
                flattened.clear();

                return result;
            }

        private:
            int bin_count;
            int n_ff;
            container::ThreadLocalWrapper<GenericDistribution3D_t> aa;
            container::ThreadLocalWrapper<GenericDistribution2D_t> aw;
            container::ThreadLocalWrapper<GenericDistribution1D_t> ww;

            // the unit-weight subsets the tasks read, owned here so that they outlive the work they were built for
            std::unordered_map<const void*, std::vector<CompactCoordinates_t>> partitions;
            std::unordered_map<const void*, CompactCoordinates_t> flattened;

            // the handles the queued tasks accumulate through; see detail::RowTarget
            detail::RowTarget<GenericDistribution3D_t, 2> target_aa(int ff1, int ff2) {
                assert(0 <= ff1 && ff1 < n_ff && 0 <= ff2 && ff2 < n_ff && "CalculatorFF::target_aa: form factor index out of bounds.");
                return detail::row_target(aa, bin_count, ff1, ff2);
            }

            detail::RowTarget<GenericDistribution2D_t, 1> target_aw(int ff1) {
                assert(0 <= ff1 && ff1 < n_ff && "CalculatorFF::target_aw: form factor index out of bounds.");
                return detail::row_target(aw, bin_count, ff1);
            }

            detail::RowTarget<GenericDistribution1D_t, 0> target_ww() {return detail::row_target(ww, bin_count);}

            /**
             * @brief detail::partition_by_ff over the active form factors, memoized on the address of @a source.
             */
            const std::vector<CompactCoordinates_t>& partition(const CompactCoordinatesFF_t& source) {
                if (auto it = partitions.find(&source); it != partitions.end()) {return it->second;}
                return partitions.emplace(&source, detail::partition_by_ff(source, n_ff)).first->second;
            }

            /**
             * @brief detail::flatten, memoized on the address of @a source.
             */
            const CompactCoordinates_t& flatten(const CompactCoordinatesFF_t& source) {
                if (auto it = flattened.find(&source); it != flattened.end()) {return it->second;}
                return flattened.emplace(&source, detail::flatten(source)).first->second;
            }
    };
}
