// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container3D.h>
#include <constants/ConstantsAxes.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <utility/TypeTraits.h>
#include <settings/Flags.h>

#include <cmath>

namespace ausaxs::hist {
    class Distribution3D;

    /**
     * @brief This is a small wrapper around the Container3D class. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution3D class.
     */
    class WeightedDistribution3D : public container::Container3D<detail::WeightedEntry> {
        public:
            using Container3D::Container3D;
            WeightedDistribution3D(const Distribution3D& other);

            /**
             * @brief Add twice the value for a given distance.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add(int x, int y, float distance, constants::axes::d_type value) {
                int i = std::round(distance*settings::flags::inv_bin_width);
                index(x, y, i).add<N>(distance, value);
            }

            /**
             * @brief Add a value for a given index.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param i The index to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add_index(int x, int y, int32_t i, const detail::WeightedEntry& value) {
                index(x, y, i) += N*value;
            }

            /**
             * @brief Increment the value for a given distance.
             *
             * @param x The form factor index.
             * @param y The second form factor index.
             * @param distance The distance to increment the value for.
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void increment(int x, int y, float distance) {
                int i = std::round(distance*settings::flags::inv_bin_width);
                index(x, y, i).increment<N>(distance);
            }

            /**
             * @brief Increment the value for a given index.
             *
             * @param x The form factor index.
             * @param y The second form factor index.
             * @param i The index to increment the value for.
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void increment_index(int x, int y, int32_t i, float distance) {
                index(x, y, i).increment<N>(distance);
            }

            template<int N = 1>
            void increment_index(int x, int y, int32_t i) {
                index(x, y, i).increment<N>();
            }
            
            /**
             * @brief Increment the value for a given linear index. 
             * 
             * @param i The index to increment.
             * @tparam N A multiplicative factor for the value. 
             */
            template<int N = 1>
            void increment_linear_index(int32_t i, float distance) {
                linear_index(i).increment<N>(distance);
            }

            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weights() const;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution3D>, "WeightedDistribution3D should support nothrow move semantics.");
}