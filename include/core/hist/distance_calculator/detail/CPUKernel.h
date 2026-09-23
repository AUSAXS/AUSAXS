// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/ThreadLocalWrapper.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/distance_calculator/detail/Evaluators.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/GeneralSettings.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>
#include <utility/observer_ptr.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ausaxs::hist::distance_calculator::detail {
    /**
     * @brief Where a queued calculation accumulates.
     */
    template<typename T>
    concept Target = std::copy_constructible<T> && requires (const T& t) {
        {t.get()} -> std::convertible_to<std::span<typename T::entry_type>>;
    };

    /**
     * @brief A Target accumulating into one row of a thread-local distribution: the whole of a 1D distribution, or
     *        the distance axis at fixed leading indices of a 2D or 3D one.
     *
     * The distance axis varies fastest in all of them, so a row is a contiguous run of bins that the evaluators can
     * write into as if it were a histogram of its own.
     */
    template<typename Distribution, std::size_t rank>
    struct RowTarget {
        using entry_type = typename Distribution::value_type;
        observer_ptr<container::ThreadLocalWrapper<Distribution>> results;
        std::array<int, rank> row;
        int bins;
        std::span<entry_type> get() const {
            auto first = std::apply([this] (auto... i) {return results->get().begin(i...);}, row);
            return {&*first, static_cast<std::size_t>(bins)};
        }
    };

    /**
     * @brief The row of @a results at the leading indices @a row, spanning @a bins bins.
     */
    template<typename Distribution, std::same_as<int>... Row>
    RowTarget<Distribution, sizeof...(Row)> row_target(container::ThreadLocalWrapper<Distribution>& results, int bins, Row... row) {
        return {&results, {row...}, bins};
    }

    /**
     * @brief Invoke @a f with the scaling factor as a compile-time constant, as the evaluators need it.
     */
    template<typename F>
    void dispatch_scaling(int scaling, F&& f) {
        switch (scaling) {
            case 1:  f(std::integral_constant<int, 1>{});  return;
            case 2:  f(std::integral_constant<int, 2>{});  return;
            case 3:  f(std::integral_constant<int, 3>{});  return;
            case 4:  f(std::integral_constant<int, 4>{});  return;
            case 5:  f(std::integral_constant<int, 5>{});  return;
            case 6:  f(std::integral_constant<int, 6>{});  return;
            case 7:  f(std::integral_constant<int, 7>{});  return;
            case 8:  f(std::integral_constant<int, 8>{});  return;
            case 9:  f(std::integral_constant<int, 9>{});  return;
            case 10: f(std::integral_constant<int, 10>{}); return;
            case 11: f(std::integral_constant<int, 11>{}); return;
            case 12: f(std::integral_constant<int, 12>{}); return;
            case 13: f(std::integral_constant<int, 13>{}); return;
            case 14: f(std::integral_constant<int, 14>{}); return;
            case 15: f(std::integral_constant<int, 15>{}); return;
            case 16: f(std::integral_constant<int, 16>{}); return;
            case 17: f(std::integral_constant<int, 17>{}); return;
            case 18: f(std::integral_constant<int, 18>{}); return;
            case 19: f(std::integral_constant<int, 19>{}); return;
            case 20: f(std::integral_constant<int, 20>{}); return;
            case 21: f(std::integral_constant<int, 21>{}); return;
            case 22: f(std::integral_constant<int, 22>{}); return;
            case 23: f(std::integral_constant<int, 23>{}); return;
            case 24: f(std::integral_constant<int, 24>{}); return;
            case 25: f(std::integral_constant<int, 25>{}); return;
            case 26: f(std::integral_constant<int, 26>{}); return;
            case 27: f(std::integral_constant<int, 27>{}); return;
            case 28: f(std::integral_constant<int, 28>{}); return;
            case 29: f(std::integral_constant<int, 29>{}); return;
            case 30: f(std::integral_constant<int, 30>{}); return;
            case 60: f(std::integral_constant<int, 60>{}); return;
            default: throw ausaxs::except::runtime_error(
                "distance_calculator::dispatch_scaling: unsupported scaling factor (" + std::to_string(scaling) + "). "
                "Supported factors are 1-30 and 60."
            );
        }
    }

    /**
     * @brief Queue the self-correlation of @a data into @a target.
     *        This is faster than cross-correlating the data with itself, since only the upper triangle has to be evaluated.
     *
     * @tparam pair_factor What each evaluated pair contributes.
     * @tparam self_factor What the zero distance of each point with itself contributes.
     *
     * The work is dispatched to the thread pool immediately; this does not wait for it. @a data must stay alive until it
     * has been waited for.
     */
    template<bool weighted_bins, bool variable_bin_width, int pair_factor, int self_factor, Target T>
    void enqueue_self(const hist::detail::CompactCoordinates<variable_bin_width>& data, T target) {
        auto* pool = utility::multi_threading::get_global_pool();
        int data_size = data.size();
        int job_size = settings::general::detail::get_job_size(data_size);

        // calculate upper triangle
        for (int i = 0; i < data_size; i+=job_size) {
            pool->detach_task(
                [&data, target, data_size, imin = i, imax = std::min(i+job_size, data_size)] () {
                    auto&& p_aa = target.get();
                    for (int i = imin; i < imax; ++i) { // atom
                        int j = i+1;                    // atom
                        for (; j+15 < data_size; j+=16) {
                            hist::detail::evaluate16<variable_bin_width, pair_factor>(p_aa, data, data, i, j);
                        }

                        for (; j+7 < data_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, pair_factor>(p_aa, data, data, i, j);
                        }

                        for (; j+3 < data_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, pair_factor>(p_aa, data, data, i, j);
                        }

                        for (; j < data_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, pair_factor>(p_aa, data, data, i, j);
                        }
                    }
                }
            );
        }

        // calculate skipped diagonal
        pool->detach_task(
            [&data, target, data_size] () {
                auto&& p_aa = target.get();
                double total_weight = 0;
                for (int i = 0; i < data_size; ++i) {
                    double weight = data.get_non_coordinate_value(i);
                    total_weight += weight*weight;
                }
                total_weight *= self_factor;

                if constexpr (weighted_bins) {
                    p_aa[0] += hist::detail::WeightedEntry(total_weight, static_cast<std::int64_t>(total_weight), 0);
                } else {
                    p_aa[0] += total_weight;
                }
            }
        );
    }

    /**
     * @brief Queue the cross-correlation of @a data_1 and @a data_2 into @a target.
     *
     * @tparam pair_factor What each evaluated pair contributes. Every pair of the two sets is visited exactly once, so a
     *                     convention that counts every unordered pair twice passes twice the scaling factor here.
     *
     * The work is dispatched to the thread pool immediately; this does not wait for it. Both sets must stay alive until
     * it has been waited for.
     */
    template<bool variable_bin_width, int pair_factor, Target T>
    void enqueue_cross(
        const hist::detail::CompactCoordinates<variable_bin_width>& data_1,
        const hist::detail::CompactCoordinates<variable_bin_width>& data_2,
        T target
    ) {
        auto* pool = utility::multi_threading::get_global_pool();
        int data_1_size = data_1.size();
        int data_2_size = data_2.size();
        int job_size = settings::general::detail::get_job_size(data_2_size);

        for (int i = 0; i < data_2_size; i+=job_size) {
            pool->detach_task(
                [&data_1, &data_2, target, data_1_size, imin = i, imax = std::min(i+job_size, data_2_size)] () {
                    auto&& p_ab = target.get();
                    for (int i = imin; i < imax; ++i) { // b
                        int j = 0;                      // a
                        for (; j+15 < data_1_size; j+=16) {
                            hist::detail::evaluate16<variable_bin_width, pair_factor>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+7 < data_1_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, pair_factor>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+3 < data_1_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, pair_factor>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j < data_1_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, pair_factor>(p_ab, data_2, data_1, i, j);
                        }
                    }
                }
            );
        }
    }

    /**
     * @brief Queue the cross-correlation of @a a and @a b into @a target, chunked over the larger of the two.
     *
     * enqueue_cross splits its second argument into the tasks it dispatches and loops the first inside each of them, so
     * a pair of very different sizes would otherwise collapse to a single task holding the whole product. The pairs are
     * the same either way, since a distance does not care which side it is read from. An empty set queues nothing.
     */
    template<bool variable_bin_width, int pair_factor, Target T>
    void enqueue_balanced_cross(
        const hist::detail::CompactCoordinates<variable_bin_width>& a,
        const hist::detail::CompactCoordinates<variable_bin_width>& b,
        T target
    ) {
        if (a.empty() || b.empty()) {return;}
        if (a.size() < b.size()) {
            enqueue_cross<variable_bin_width, pair_factor>(a, b, target);
        } else {
            enqueue_cross<variable_bin_width, pair_factor>(b, a, target);
        }
    }

    /**
     * @brief Split @a source into one unit-weight coordinate set per form factor type, indexed by type.
     *        Types with no atoms get an empty set.
     *
     * The form factor of an atom only decides which histogram its pairs land in, never the distance itself, so a
     * calculation restricted to fixed form factor types is an ordinary weighted one over these subsets. The unit
     * weights make every pair count once; the form factor amplitudes are applied later, by the intensity calculator.
     */
    template<bool variable_bin_width>
    std::vector<hist::detail::CompactCoordinates<variable_bin_width>> partition_by_ff(
        const hist::detail::CompactCoordinatesFF<variable_bin_width>& source, int n_ff
    ) {
        std::vector<int> counts(n_ff, 0);
        for (int i = 0; i < source.size(); ++i) {++counts[source.get_ff_type(i)];}

        std::vector<hist::detail::CompactCoordinates<variable_bin_width>> parts(n_ff);
        for (int ff = 0; ff < n_ff; ++ff) {parts[ff].resize(counts[ff]);}

        std::vector<int> filled(n_ff, 0);
        for (int i = 0; i < source.size(); ++i) {
            int ff = source.get_ff_type(i);
            int k = filled[ff]++;
            parts[ff].set_position(k, source.position(i));
            parts[ff].get_non_coordinate_value(k) = 1;
        }
        return parts;
    }

    /**
     * @brief The whole of @a source as one unit-weight coordinate set, disregarding its form factor types.
     *        See partition_by_ff.
     */
    template<bool variable_bin_width>
    hist::detail::CompactCoordinates<variable_bin_width> flatten(const hist::detail::CompactCoordinatesFF<variable_bin_width>& source) {
        hist::detail::CompactCoordinates<variable_bin_width> whole;
        whole.resize(source.size());
        for (int i = 0; i < source.size(); ++i) {
            whole.set_position(i, source.position(i));
            whole.get_non_coordinate_value(i) = 1;
        }
        return whole;
    }
}
