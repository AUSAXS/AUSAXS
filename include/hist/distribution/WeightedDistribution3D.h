#pragma once

#include <container/Container3D.h>
#include <constants/Axes.h>

namespace hist {
    class Distribution3D;

    /**
     * @brief This is a small wrapper around the Container3D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution3D : public container::Container3D<constants::axes::d_type> {
        public:
            using Container3D::Container3D;
            WeightedDistribution3D(Distribution3D&& other);

            /**
             * @brief Add a value for a given distance.
             */
            void add(int x, int y, float distance, constants::axes::d_type value);
    };
}