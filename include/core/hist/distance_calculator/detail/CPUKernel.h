// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/detail/Evaluators.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <settings/GeneralSettings.h>
#include <utility/MultiThreading.h>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <span>
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
     * @brief The zero distance of every point of @a data with itself, i.e. the sum of their squared weights.
     *
     * @tparam unit_weights Whether every point weighs 1, in which case the stored weights are never read.
     */
    template<bool unit_weights, bool variable_bin_width>
    double self_weight(const hist::detail::CompactCoordinates<variable_bin_width>& data) {
        if constexpr (unit_weights) {return data.size();}
        double total_weight = 0;
        for (int i = 0; i < data.size(); ++i) {
            double weight = data.get_weight(i);
            total_weight += weight*weight;
        }
        return total_weight;
    }

    /**
     * @brief Queue the self-correlation of @a data into @a target: each distinct pair counted twice, once in either order,
     *        and the zero distance of every point with itself once.
     *
     * @tparam unit_weights Whether every point weighs 1, in which case the stored weights are never read.
     */
    template<bool unit_weights, bool variable_bin_width, Target T>
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
                            hist::detail::evaluate16<variable_bin_width, 2, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j+7 < data_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, 2, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j+3 < data_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, 2, unit_weights>(p_aa, data, data, i, j);
                        }

                        for (; j < data_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, 2, unit_weights>(p_aa, data, data, i, j);
                        }
                    }
                }
            );
        }

        // calculate skipped diagonal
        pool->detach_task(
            [&data, target] () {
                auto&& p_aa = target.get();
                double weight = self_weight<unit_weights>(data);
                if constexpr (std::is_same_v<typename T::entry_type, hist::detail::WeightedEntry>) {
                    p_aa[0] += hist::detail::WeightedEntry(weight, static_cast<std::int64_t>(weight), 0);
                } else {
                    p_aa[0] += weight;
                }
            }
        );
    }

    /**
     * @brief Queue the cross-correlation of @a data_1 and @a data_2 into @a target, each pair counted once.
     *
     * @tparam unit_weights Whether every point weighs 1, in which case the stored weights are never read.
     */
    template<bool unit_weights, bool variable_bin_width, Target T>
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
                            hist::detail::evaluate16<variable_bin_width, 1, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+7 < data_1_size; j+=8) {
                            hist::detail::evaluate8<variable_bin_width, 1, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j+3 < data_1_size; j+=4) {
                            hist::detail::evaluate4<variable_bin_width, 1, unit_weights>(p_ab, data_2, data_1, i, j);
                        }

                        for (; j < data_1_size; ++j) {
                            hist::detail::evaluate1<variable_bin_width, 1, unit_weights>(p_ab, data_2, data_1, i, j);
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
    template<bool unit_weights, bool variable_bin_width, Target T>
    void enqueue_balanced_cross(
        const hist::detail::CompactCoordinates<variable_bin_width>& a,
        const hist::detail::CompactCoordinates<variable_bin_width>& b,
        T target
    ) {
        if (a.size() < b.size()) {
            enqueue_cross<unit_weights>(a, b, target);
        } else {
            enqueue_cross<unit_weights>(b, a, target);
        }
    }
}
