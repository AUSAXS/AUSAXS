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
             * @brief Add a value for a given distance.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value) {
                index(x, y, std::round(distance)) += N*value;
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
            void add_index(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value) {
                index(x, y, i) += N*value;
            }
    };
    static_assert(supports_nothrow_move_v<Distribution3D>, "Distribution3D should support nothrow move semantics.");
}