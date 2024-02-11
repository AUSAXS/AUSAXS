#pragma once

#include <container/Container1D.h>
#include <constants/Axes.h>
#include <hist/distribution/detail/WeightedEntry.h>

namespace hist {
    class Distribution1D;

    /**
     * @brief This is a small wrapper around the Container1D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class WeightedDistribution1D : public container::Container1D<detail::WeightedEntry> {
        public:
            using Container1D::Container1D;
            WeightedDistribution1D(const Distribution1D& other);
            WeightedDistribution1D(const std::vector<constants::axes::d_type>& bins);

            /**
             * @brief Convert this distribution to a vector format. 
             *        This is equivalent to get_content() for this class. 
             */
            std::vector<constants::axes::d_type> as_vector() const;

            /**
             * @brief Get the bin values from this distribution.
             */
            std::vector<constants::axes::d_type> get_content() const;

            /**
             * @brief Get a bin value from this distribution.
             */
            constants::axes::d_type& get_content(int i);
            const constants::axes::d_type& get_content(int i) const; // @copydoc get_content(int i)

            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weighted_axis() const;

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, constants::axes::d_type value);

            WeightedDistribution1D& operator+=(const WeightedDistribution1D& other);
            WeightedDistribution1D& operator-=(const WeightedDistribution1D& other);
    };
}