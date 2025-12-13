// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container3D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <constants/ConstantsAxes.h>
#include <utility/TypeTraits.h>

#include <cmath>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container3D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution3D : public container::Container3D<constants::axes::d_type> {
        public:
            using Container3D::Container3D;
            explicit Distribution3D(const WeightedDistribution3D& other);

            /**
             * @brief Add a value for a given bin index.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param i The bin index to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add_index(int x, int y, int32_t i, constants::axes::d_type value) {
                index(x, y, i) += N*value;
            }

            template<int N = 1>
            void increment_index(int x, int y, int32_t i) {
                index(x, y, i) += N;
            }
            
            /**
             * @brief Increment the value for a given linear index. 
             * 
             * @param xy The combined form factor index.
             * @param i The index to increment.
             * @tparam N A multiplicative factor for the value. 
             */
            template<int N = 1>
            void increment_linear_index(int xy, int32_t i) {
                linear_index(xy, i) += N;                
            }
    };
    static_assert(supports_nothrow_move_v<Distribution3D>, "Distribution3D should support nothrow move semantics.");
}