#pragma once

#include <container/Container2D.h>
#include <hist/distribution/WeightedDistribution2D.h>
#include <constants/Axes.h>

namespace hist {
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
             */
            void add(unsigned int x, float distance, constants::axes::d_type value);
            void add(unsigned int x, int32_t i, constants::axes::d_type value);
    };
}