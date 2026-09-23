// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/CPUKernel.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <utility/MultiThreading.h>

#include <cassert>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Evaluates the form factor-resolved pairwise distance histograms on the global thread pool.
     *
     * The form factor of an atom only decides *which* histogram its pairs land in, never the distance itself, so a
     * calculation restricted to two form factor types is an ordinary pairwise distance calculation whose result
     * happens to be one row of a larger distribution. This calculator therefore takes the atoms already split by form
     * factor, and hands the subsets to the same detail::enqueue_self and detail::enqueue_cross the weight-based
     * calculator uses, with a target naming the row to accumulate into.
     *
     * The rows are rows of one distribution per thread, not histograms of their own, so a worker thread keeps
     * accumulating into the single array it owns and everything is merged once, at the end.
     *
     * Contributions are counts: the inputs are expected to carry unit weights (see hist::detail::factory::construct_unit_weight),
     * and the form factor amplitudes are applied later, by the intensity calculator. The caller must keep all submitted
     * data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class CalculatorFF {
        using GenericDistribution1D_t = typename hist::GenericDistribution1D<weighted_bins>::type;
        using GenericDistribution2D_t = typename hist::GenericDistribution2D<weighted_bins>::type;
        using GenericDistribution3D_t = typename hist::GenericDistribution3D<weighted_bins>::type;
        using CompactCoordinates_t = hist::detail::CompactCoordinates<variable_bin_width>;
        using PartitionedCoordinates_t = std::vector<CompactCoordinates_t>;
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
             *        @a a holds one set per active form factor, indexed by type.
             *
             * Every unordered pair is counted twice, once in each of the two form factor orders it could be read in;
             * a pair of unlike types lands wholly in one of them rather than being split between the two. Only the
             * sum over both form factor indices is ever read, so this is the same distribution.
             */
            void enqueue_self_by_ff(const PartitionedCoordinates_t& a) {
                assert(static_cast<int>(a.size()) == n_ff && "CalculatorFF::enqueue_self_by_ff: expected one set per active form factor.");
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    if (a[ff1].empty()) {continue;}
                    detail::enqueue_self<weighted_bins, variable_bin_width, 2, 1>(a[ff1], target_aa(ff1, ff1));
                    for (int ff2 = ff1+1; ff2 < n_ff; ++ff2) {
                        detail::enqueue_balanced_cross<variable_bin_width, 2>(a[ff1], a[ff2], target_aa(ff1, ff2));
                    }
                }
            }

            /**
             * @brief Queue the cross-correlation of @a a and @a b, resolved by form factor on the @a a side only,
             *        into the aw result. Each pair is counted once, as that convention expects.
             */
            void enqueue_cross_by_ff(const PartitionedCoordinates_t& a, const CompactCoordinates_t& b) {
                assert(static_cast<int>(a.size()) == n_ff && "CalculatorFF::enqueue_cross_by_ff: expected one set per active form factor.");
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    detail::enqueue_balanced_cross<variable_bin_width, 1>(a[ff1], b, target_aw(ff1));
                }
            }

            /**
             * @brief Queue the self-correlation of @a b, ignoring form factors, into the ww result.
             *        Every unordered pair is counted twice, as that convention expects.
             */
            void enqueue_self_flat(const CompactCoordinates_t& b) {
                if (b.empty()) {return;}
                detail::enqueue_self<weighted_bins, variable_bin_width, 2, 1>(b, target_ww());
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

                return result;
            }

        private:
            int bin_count;
            int n_ff;
            container::ThreadLocalWrapper<GenericDistribution3D_t> aa;
            container::ThreadLocalWrapper<GenericDistribution2D_t> aw;
            container::ThreadLocalWrapper<GenericDistribution1D_t> ww;

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
    };
}
