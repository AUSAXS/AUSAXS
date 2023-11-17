#pragma once

#include <container/Container2D.h>
#include <constants/Axes.h>

namespace hist {
    class WeightedDistribution2D;

    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution2D : public container::Container2D<constants::axes::d_type> {
        public:
            using Container2D::Container2D;
            Distribution2D(WeightedDistribution2D&& other);

            /**
             * @brief Add a value for a given distance.
             */
            void add(unsigned int x, float distance, constants::axes::d_type value);
            void add(unsigned int x, int32_t i, constants::axes::d_type value);
    };
}