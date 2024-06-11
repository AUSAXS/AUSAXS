#pragma once

#include <container/Container3D.h>
#include <hist/distribution/WeightedDistribution3D.h>
#include <constants/Axes.h>

namespace hist {
    /**
     * @brief This is a small wrapper around the Container3D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution3D : public container::Container3D<constants::axes::d_type> {
        public:
            using Container3D::Container3D;
            explicit Distribution3D(const WeightedDistribution3D& other);

            /**
             * @brief Add a value for a given distance.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value);

            /**
             * @brief Add twice the value for a given distance.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param distance The distance to add the value to.
             * @param value The value to add.
             */
            void add2(unsigned int x, unsigned int y, float distance, constants::axes::d_type value);

            /**
             * @brief Add a value for a given index.
             * 
             * @param x The first form factor index.
             * @param y The second form factor index.
             * @param i The index to add the value to.
             * @param value The value to add.
             */
            void add(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value);
    };
}