// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container2D.h>
#include <constants/ConstantsAxes.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <utility/TypeTraits.h>
#include <settings/Flags.h>

#include <cmath>

namespace ausaxs::hist {
    class Distribution2D;

    /**
     * @brief This is a small wrapper around the Container2D class. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution2D class.
     */
    class WeightedDistribution2D : public container::Container2D<detail::WeightedEntry> {
        public:
            using Container2D::Container2D;
            WeightedDistribution2D(const Distribution2D& other);

            /**
             * @brief Add a value for a given index.
             * 
             * @param x The form factor index. 
             * @param i The index to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add_index(int x, int32_t i, const detail::WeightedEntry& value) {
                index(x, i) += N*value;
            }

            /**
             * @brief Increment the value for a given index.
             *
             * @param x The form factor index.
             * @param i The index to increment the value for.
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void increment_bin(int x, int32_t i, float distance) {
                index(x, i).increment<N>(distance);
            }

            template<int N = 1>
            void increment_bin(int x, int32_t i) {
                index(x, i).increment<N>();
            }

            /**
             * @brief Increment the value for a given linear index. 
             * 
             * @param i The index to increment.
             * @tparam N A multiplicative factor for the value. 
             */
            template<int N = 1>
            void increment_linear(int32_t i, float distance) {
                linear_index(i).increment<N>(distance);
            }
            
            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weighted_axis() const;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution2D>, "WeightedDistribution2D should support nothrow move semantics.");
}