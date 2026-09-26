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
     * @brief Where a queued calculation accumulates, see HistogramStore::Target.
     *
     * get() returns the calling thread's own bins.
     */
    template<typename T>
    concept Target = std::copy_constructible<T> && requires (const T& t) {
        {t.get()} -> std::convertible_to<std::span<typename T::entry_type>>;
    };

    /**
     * @brief Invoke @a f with the pair factor as a compile-time constant, as the evaluators need it.
     */
    template<typename F>
    void dispatch_scaling(int pair_factor, F&& f) {
        switch (pair_factor) {
            case 1:   f(std::integral_constant<int, 1>{});   return;
            case 2:   f(std::integral_constant<int, 2>{});   return;
            case 4:   f(std::integral_constant<int, 4>{});   return;
            case 6:   f(std::integral_constant<int, 6>{});   return;
            case 8:   f(std::integral_constant<int, 8>{});   return;
            case 10:  f(std::integral_constant<int, 10>{});  return;
            case 12:  f(std::integral_constant<int, 12>{});  return;
            case 14:  f(std::integral_constant<int, 14>{});  return;
            case 16:  f(std::integral_constant<int, 16>{});  return;
            case 18:  f(std::integral_constant<int, 18>{});  return;
            case 20:  f(std::integral_constant<int, 20>{});  return;
            case 22:  f(std::integral_constant<int, 22>{});  return;
            case 24:  f(std::integral_constant<int, 24>{});  return;
            case 26:  f(std::integral_constant<int, 26>{});  return;
            case 28:  f(std::integral_constant<int, 28>{});  return;
            case 30:  f(std::integral_constant<int, 30>{});  return;
            case 32:  f(std::integral_constant<int, 32>{});  return;
            case 34:  f(std::integral_constant<int, 34>{});  return;
            case 36:  f(std::integral_constant<int, 36>{});  return;
            case 38:  f(std::integral_constant<int, 38>{});  return;
            case 40:  f(std::integral_constant<int, 40>{});  return;
            case 42:  f(std::integral_constant<int, 42>{});  return;
            case 44:  f(std::integral_constant<int, 44>{});  return;
            case 46:  f(std::integral_constant<int, 46>{});  return;
            case 48:  f(std::integral_constant<int, 48>{});  return;
            case 50:  f(std::integral_constant<int, 50>{});  return;
            case 52:  f(std::integral_constant<int, 52>{});  return;
            case 54:  f(std::integral_constant<int, 54>{});  return;
            case 56:  f(std::integral_constant<int, 56>{});  return;
            case 58:  f(std::integral_constant<int, 58>{});  return;
            case 60:  f(std::integral_constant<int, 60>{});  return;
            case 120: f(std::integral_constant<int, 120>{}); return;
            default: throw ausaxs::except::runtime_error(
                "distance_calculator::dispatch_scaling: unsupported pair factor (" + std::to_string(pair_factor) + "). "
                "Supported factors are 1, the even numbers 2-60, and 120."
            );
        }
    }

    /**
     * @brief Queue the self-correlation of @a data into @a target.
     *
     * @tparam pair_factor What each evaluated pair contributes.
     * @tparam self_factor What the zero distance of each point with itself contributes.
     * @tparam unit_weights Whether every point weighs 1, in which case the stored weights are never read.
     */
    template<int pair_factor, int self_factor, bool unit_weights, bool variable_bin_width, Target T>
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
                            hist::detail::evaluate16<variable_bin_width, pair_factor, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j+7 < data_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, pair_factor, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j+3 < data_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, pair_factor, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j < data_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, pair_factor, unit_weights>(p_aa, data, data, i, j);
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
                if constexpr (unit_weights) {
                    total_weight = data_size;
                } else {
                    for (int i = 0; i < data_size; ++i) {
                        double weight = data.get_weight(i);
                        total_weight += weight*weight;
                    }
                }
                total_weight *= self_factor;

                if constexpr (std::is_same_v<typename T::entry_type, hist::detail::WeightedEntry>) {
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
     * @tparam pair_factor What each evaluated pair contributes.
     * @tparam unit_weights Whether every point weighs 1, in which case the stored weights are never read.
     */
    template<int pair_factor, bool unit_weights, bool variable_bin_width, Target T>
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
                            hist::detail::evaluate16<variable_bin_width, pair_factor, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+7 < data_1_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, pair_factor, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+3 < data_1_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, pair_factor, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j < data_1_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, pair_factor, unit_weights>(p_ab, data_2, data_1, i, j);
                        }
                    }
                }
            );
        }
    }

    /**
     * @brief Queue the cross-correlation of @a a and @a b into @a target, chunked over the larger of the two.
     *        This leads to a more balanced work distribution for strongly asymmetric sizes. 
     */
    template<int pair_factor, bool unit_weights, bool variable_bin_width, Target T>
    void enqueue_balanced_cross(
        const hist::detail::CompactCoordinates<variable_bin_width>& a,
        const hist::detail::CompactCoordinates<variable_bin_width>& b,
        T target
    ) {
        if (a.size() < b.size()) {
            enqueue_cross<pair_factor, unit_weights>(a, b, target);
        } else {
            enqueue_cross<pair_factor, unit_weights>(b, a, target);
        }
    }
}
