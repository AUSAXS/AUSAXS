#pragma once

#include <container/Container3D.h>
#include <constants/ConstantsAxes.h>
#include <hist/distribution/detail/WeightedEntry.h>
#include <utility/TypeTraits.h>

#include <cmath>

namespace ausaxs::hist {
    class Distribution3D;

    /**
     * @brief This is a small wrapper around the Container3D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis. Anything added to this
     *        distribution will be tracked by the WeightedDistribution class, which may add
     *        a significant overhead compared to a pure Distribution1D class.
     */
    class WeightedDistribution3D : public container::Container3D<detail::WeightedEntry> {
        public:
            using Container3D::Container3D;
            WeightedDistribution3D(const Distribution3D& other);

            /**
             * @brief Add twice the value for a given distance.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             *
             * @tparam N A multiplicative factor for the value.
             */
            template<int N = 1>
            void add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value) {
                int i = std::round(distance*constants::axes::d_inv_width);
                index(x, y, i).add<N>(distance, value);
            }

            /**
             * @brief Extract the weights from this distribution.
             */
            std::vector<double> get_weights() const;
    };
    static_assert(supports_nothrow_move_v<WeightedDistribution3D>, "WeightedDistribution3D should support nothrow move semantics.");
}