#pragma once

#include <container/Container2D.h>
#include <constants/Axes.h>
#include <hist/distribution/detail/Entry.h>

namespace hist {
    class Distribution2D;

    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution2D : public container::Container2D<detail::Entry> {
        public:
            using Container2D::Container2D;
            WeightedDistribution2D(const Distribution2D& other);

            /**
             * @brief Add a value for a given distance.
             */
            void add(int x, float distance, constants::axes::d_type value);

            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weighted_axis() const;
    };
}