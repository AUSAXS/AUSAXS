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

            template<int N = 1>
            void add_index(int x, int y, int32_t i, float distance, float weight) {
                index(x, y, i).add<N>(distance, weight);
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

            /**
             * @brief Increment the value using a linear 2D form factor bin index.
             *
             * @param xy The linear 2D form factor bin index (computed as ff2 + ff1*N_ff).
             * @param i The distance bin index.
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void increment_linear_index(int32_t xy, int32_t i, float distance) {
                linear_index(xy, i).increment<N>(distance);
            }
            
            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weights() const;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution3D>, "WeightedDistribution3D should support nothrow move semantics.");
}