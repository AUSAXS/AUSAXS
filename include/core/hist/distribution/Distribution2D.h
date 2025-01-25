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
             * @brief Add a value for a given distance.
             * 
             * @param x The form factor index. 
             * @param distance The distance to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add(unsigned int x, float distance, constants::axes::d_type value) {
                index(x, std::round(distance)) += N*value;
            }

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
            void add_index(unsigned int x, int32_t i, constants::axes::d_type value) {
                index(x, i) += N*value;
            }
    };
    static_assert(supports_nothrow_move_v<Distribution2D>, "Distribution2D should support nothrow move semantics.");
}