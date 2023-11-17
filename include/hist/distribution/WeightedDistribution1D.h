#pragma once

#include <container/Container1D.h>
#include <constants/Axes.h>

namespace hist {
    class Distribution1D;

    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution1D : public container::Container1D<constants::axes::d_type> {
        public:
            using Container1D::Container1D;
            WeightedDistribution1D(Distribution1D&& other);

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, constants::axes::d_type value);
    };
}