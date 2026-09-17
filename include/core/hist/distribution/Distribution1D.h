// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <constants/ConstantsAxes.h>
#include <container/Container1D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution1D : public container::Container1D<double> {
        public:
            using Container1D::Container1D;
            explicit Distribution1D(const WeightedDistribution1D& other);

            /**
             * @brief Convert this distribution to a vector format. 
             */
            std::vector<double> as_vector() const;

            /**
             * @brief Get the bin values from this distribution.
             */
            const std::vector<double>& get_content() const;

            /**
             * @brief Get a bin value from this distribution.
             */
            double& get_content(int i);
            const double& get_content(int i) const; // @copydoc get_content(int i)

            /**
             * @brief Set the value of the ith bin.
             */
            void set_content(int i, double value);

            /**
             * @brief Add a value for a given bin index.
             * 
             * @param i The bin index to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add_index(int32_t i, double value) {
                index(i) += N*value;
            }

            template<int N = 1>
            void increment_index(int32_t i) {
                index(i) += N;
            }

            /**
             * @brief Increment the value for a given linear index. 
             * 
             * @param i The index to increment.
             * @tparam N A multiplicative factor for the value. 
             */
            template<int N = 1>
            void increment_linear_index(int32_t i) {
                linear_index(i) += N;
            }

            /**
             * @brief Clear the value for a given distance.
             * 
             * @param distance The index to clear.
             */
            void clear(int32_t i);

            Distribution1D& operator+=(const Distribution1D& rhs);
            Distribution1D& operator-=(const Distribution1D& rhs);
    };
    Distribution1D operator*(double factor, Distribution1D dist);
    static_assert(supports_nothrow_move_v<Distribution1D>, "Distribution1D should support nothrow move semantics.");
}