// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/Evaluators.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/GeneralSettings.h>
#include <utility/Exceptions.h>
#include <utility/MultiThreading.h>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <type_traits>

namespace ausaxs::hist::distance_calculator::detail {
    /**
     * @brief Where a queued calculation accumulates.
     */
    template<typename T>
    concept Target = std::copy_constructible<T> && requires (const T& t) {
        {t.get()} -> std::convertible_to<std::span<typename T::entry_type>>;
    };

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
}
