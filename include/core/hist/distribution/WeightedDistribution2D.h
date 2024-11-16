#pragma once

#include <container/Container2D.h>
#include <constants/ConstantsAxes.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <utility/TypeTraits.h>

namespace ausaxs::hist {
    class Distribution2D;

    /**
     * @brief This is a small wrapper around the Container2D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution2D : public container::Container2D<detail::WeightedEntry> {
        public:
            using Container2D::Container2D;
            WeightedDistribution2D(const Distribution2D& other);

            /**
             * @brief Add a value for a given distance.
             * 
             * @param x The first form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add(unsigned int x, float distance, constants::axes::d_type value);

            /**
             * @brief Add twice the value for a given distance.
             * 
             * @param x The first form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add2(unsigned int x, float distance, constants::axes::d_type value);

            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weighted_axis() const;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution2D>, "WeightedDistribution2D should support nothrow move semantics.");
}