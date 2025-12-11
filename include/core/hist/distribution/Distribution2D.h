// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container2D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <constants/ConstantsAxes.h>
#include <utility/TypeTraits.h>

#include <cmath>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution2D : public container::Container2D<constants::axes::d_type> {
        public:
            using Container2D::Container2D;
            explicit Distribution2D(const WeightedDistribution2D& other);

            /**
             * @brief Get a bin value from this distribution.
             */
            constants::axes::d_type& get_content(int i, int j);
            const constants::axes::d_type& get_content(int i, int j) const; // @copydoc get_content(int i, int j)

            /**
             * @brief Add a value for a given bin index.
             * 
             * @param x The form factor index. 
             * @param i The bin index to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add_index(int x, int32_t i, constants::axes::d_type value) {
                index(x, i) += N*value;
            }
            
            /**
             * @brief Increment the value for a given bin index.dex.
             * 
             * @param x The form factor index. 
             * @param i The bin index to increment.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void increment(int x, int32_t i) {
                index(x, i) += N;
            }

            template<int N = 1>
            void increment_bin(int x, int32_t i) {
                index(x, i) += N;
            }

            /**
             * @brief Increment the value for a given linear index. 
             * 
             * @param i The index to increment.
             * @tparam N A multiplicative factor for the value. 
             */
            template<int N = 1>
            void increment_linear_bin(int32_t i) {
                linear_index(i) += N;
            }
    };
    static_assert(supports_nothrow_move_v<Distribution2D>, "Distribution2D should support nothrow move semantics.");
}