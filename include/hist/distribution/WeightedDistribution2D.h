#pragma once

#include <container/Container2D.h>
#include <constants/Axes.h>

namespace hist {
    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution2D : public container::Container2D<constants::axes::d_type> {
        public:
            using Container2D::Container2D;

            /**
             * @brief Add a value for a given distance.
             */
            void add(int x, float distance, constants::axes::d_type value);
    };
}