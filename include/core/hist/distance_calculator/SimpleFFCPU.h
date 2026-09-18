// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <form_factor/FormFactorType.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/distance_calculator/detail/AccumulationTasks.h>
#include <hist/distribution/GenericDistribution1D.h>
#include <hist/distribution/GenericDistribution2D.h>
#include <hist/distribution/GenericDistribution3D.h>
#include <utility/MultiThreading.h>
#include <utility/observer_ptr.h>

#include <cassert>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

namespace ausaxs::hist::distance_calculator {
    /**
     * @brief Evaluates the form factor-resolved pairwise distance histograms on the global thread pool.
     *
     * The form factor of an atom only decides *which* histogram its pairs land in, never the distance itself, so a
     * calculation restricted to two form factor types is an ordinary pairwise distance calculation whose result
     * happens to be one row of a larger distribution. This kernel therefore splits each input by form factor and
     * hands the resulting subsets to the same detail::enqueue_self and detail::enqueue_cross the weight-based kernel
     * uses, with a target naming the row to accumulate into.
     *
     * The rows are rows of one distribution per thread, not histograms of their own, so a worker thread keeps
     * accumulating into the single array it owns and everything is merged once, at the end.
     *
     * Contributions are counts: the subsets carry unit weights, and the form factor amplitudes are applied later,
     * by the intensity calculator. The caller must keep all submitted data alive until run() returns.
     */
    template<bool weighted_bins, bool variable_bin_width>
    class SimpleFFCPU {
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
             * @brief Construct a kernel whose result histograms span @a bin_count bins.
             */
            explicit SimpleFFCPU(int bin_count)
                : bin_count(bin_count), n_ff(form_factor::get_active_count()),
                  aa(n_ff, n_ff, bin_count), aw(n_ff, bin_count), ww(bin_count) {}

            /**
             * @brief Drain the thread pool before any of this object's state is released.
             *        See SimpleCPU's destructor; the same reasoning applies.
             */
            ~SimpleFFCPU() {
                auto* pool = utility::multi_threading::get_global_pool();
                pool->purge();
                pool->wait();
            }

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
                        if (parts[ff2].empty()) {continue;}
                        enqueue_balanced_cross<2>(parts[ff1], parts[ff2], target_aa(ff1, ff2));
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
                if (whole.empty()) {return;}
                for (int ff1 = 0; ff1 < n_ff; ++ff1) {
                    if (parts[ff1].empty()) {continue;}
                    enqueue_balanced_cross<1>(parts[ff1], whole, target_aw(ff1));
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
            /**
             * @brief The handles the queued tasks accumulate through; see detail::Target.
             *
             * The distance axis varies fastest in all three results, so a fixed form factor index names a contiguous
             * run of bins that the evaluators can write into as if it were a histogram of its own.
             */
            struct TargetAA {
                using entry_type = typename GenericDistribution3D_t::value_type;
                observer_ptr<container::ThreadLocalWrapper<GenericDistribution3D_t>> results;
                int ff1, ff2, bins;
                std::span<entry_type> get() const {
                    return {&*results->get().begin(ff1, ff2), static_cast<std::size_t>(bins)};
                }
            };

            struct TargetAW {
                using entry_type = typename GenericDistribution2D_t::value_type;
                observer_ptr<container::ThreadLocalWrapper<GenericDistribution2D_t>> results;
                int ff1, bins;
                std::span<entry_type> get() const {
                    return {&*results->get().begin(ff1), static_cast<std::size_t>(bins)};
                }
            };

            struct TargetWW {
                using entry_type = typename GenericDistribution1D_t::value_type;
                observer_ptr<container::ThreadLocalWrapper<GenericDistribution1D_t>> results;
                int bins;
                std::span<entry_type> get() const {
                    return {&*results->get().begin(), static_cast<std::size_t>(bins)};
                }
            };

            int bin_count;
            int n_ff;
            container::ThreadLocalWrapper<GenericDistribution3D_t> aa;
            container::ThreadLocalWrapper<GenericDistribution2D_t> aw;
            container::ThreadLocalWrapper<GenericDistribution1D_t> ww;

            // the unit-weight subsets the tasks read, owned here so that they outlive the work they were built for
            std::unordered_map<const void*, std::vector<CompactCoordinates_t>> partitions;
            std::unordered_map<const void*, CompactCoordinates_t> flattened;

            /**
             * @brief Cross-correlate two sets, chunked over the larger of the two.
             *
             * detail::enqueue_cross splits its second argument into the tasks it dispatches and loops the first inside
             * each of them, so a pair of very different sizes would otherwise collapse to a single task holding the
             * whole product. The pairs are the same either way, since a distance does not care which side it is read from.
             */
            template<int pair_factor, detail::Target T>
            void enqueue_balanced_cross(const CompactCoordinates_t& a, const CompactCoordinates_t& b, T target) {
                if (a.size() < b.size()) {
                    detail::enqueue_cross<variable_bin_width, pair_factor>(a, b, target);
                } else {
                    detail::enqueue_cross<variable_bin_width, pair_factor>(b, a, target);
                }
            }

            TargetAA target_aa(int ff1, int ff2) {
                assert(0 <= ff1 && ff1 < n_ff && 0 <= ff2 && ff2 < n_ff && "SimpleFFCPU::target_aa: form factor index out of bounds.");
                return {&aa, ff1, ff2, bin_count};
            }

            TargetAW target_aw(int ff1) {
                assert(0 <= ff1 && ff1 < n_ff && "SimpleFFCPU::target_aw: form factor index out of bounds.");
                return {&aw, ff1, bin_count};
            }

            TargetWW target_ww() {return {&ww, bin_count};}

            /**
             * @brief Split @a source into one unit-weight coordinate set per active form factor, memoized on its address.
             *        Form factors with no atoms get an empty set, which the callers skip.
             */
            const std::vector<CompactCoordinates_t>& partition(const CompactCoordinatesFF_t& source) {
                if (auto it = partitions.find(&source); it != partitions.end()) {return it->second;}

                std::vector<int> counts(n_ff, 0);
                for (int i = 0; i < source.size(); ++i) {++counts[source.get_ff_type(i)];}

                std::vector<CompactCoordinates_t> parts(n_ff);
                for (int ff = 0; ff < n_ff; ++ff) {parts[ff].resize(counts[ff]);}

                std::vector<int> filled(n_ff, 0);
                for (int i = 0; i < source.size(); ++i) {
                    int ff = source.get_ff_type(i);
                    int k = filled[ff]++;
                    parts[ff].set_position(k, source.position(i));
                    parts[ff].get_non_coordinate_value(k) = 1;
                }
                return partitions.emplace(&source, std::move(parts)).first->second;
            }

            /**
             * @brief The whole of @a source as one unit-weight coordinate set, memoized on its address.
             */
            const CompactCoordinates_t& flatten(const CompactCoordinatesFF_t& source) {
                if (auto it = flattened.find(&source); it != flattened.end()) {return it->second;}

                CompactCoordinates_t whole;
                whole.resize(source.size());
                for (int i = 0; i < source.size(); ++i) {
                    whole.set_position(i, source.position(i));
                    whole.get_non_coordinate_value(i) = 1;
                }
                return flattened.emplace(&source, std::move(whole)).first->second;
            }
    };
}
