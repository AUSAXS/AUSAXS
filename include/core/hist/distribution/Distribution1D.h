#pragma once

#include <container/Container1D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <constants/ConstantsAxes.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution1D : public container::Container1D<constants::axes::d_type> {
        public:
            using Container1D::Container1D;
            explicit Distribution1D(const WeightedDistribution1D& other);

            /**
             * @brief Convert this distribution to a vector format. 
             */
            std::vector<constants::axes::d_type> as_vector() const;

            /**
             * @brief Get the bin values from this distribution.
             */
            const std::vector<constants::axes::d_type>& get_content() const;

            /**
             * @brief Get a bin value from this distribution.
             */
            constants::axes::d_type& get_content(int i);
            const constants::axes::d_type& get_content(int i) const; // @copydoc get_content(int i)

            /**
             * @brief Add a value for a given distance.
             * 
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add(float distance, constants::axes::d_type value);

            /**
             * @brief Add twice the value for a given distance.
             * 
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add2(float distance, constants::axes::d_type value);

            /**
             * @brief Add a value for a given index.
             * 
             * @param i The index to add the value to.
             * @param value The value to add.
             */
            void add_index(int32_t i, constants::axes::d_type value);

            /**
             * @brief Clear the value for a given distance.
             * 
             * @param distance The index to clear.
             */
            void clear(int32_t i);

            Distribution1D& operator+=(const Distribution1D& other);
            Distribution1D& operator-=(const Distribution1D& other);
    };
    Distribution1D operator*(double factor, Distribution1D dist);
    static_assert(supports_nothrow_move_v<Distribution1D>, "Distribution1D should support nothrow move semantics.");
}