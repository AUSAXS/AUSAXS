#pragma once

#include <container/Container1D.h>
#include <constants/Axes.h>
#include <hist/distribution/detail/Entry.h>

namespace hist {
    class Distribution1D;

    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class WeightedDistribution1D : public container::Container1D<detail::Entry> {
        public:
            using Container1D::Container1D;
            WeightedDistribution1D(Distribution1D& other);

            std::vector<constants::axes::d_type> get_bins() const;

            std::vector<double> get_weights() const;

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, constants::axes::d_type value);
    };
}