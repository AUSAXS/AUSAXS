#pragma once

#include <container/Container1D.h>
#include <hist/distribution/WeightedDistribution1D.h>
#include <constants/Axes.h>

#include <cmath>

namespace hist {
    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution1D : public container::Container1D<constants::axes::d_type> {
        public:
            using Container1D::Container1D;
            Distribution1D(const WeightedDistribution1D& other);

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, constants::axes::d_type value);

            /**
             * @brief Add a value for a given distance.
             */
            void add(int32_t i, constants::axes::d_type value);
    };
}