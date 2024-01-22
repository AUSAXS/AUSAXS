#pragma once

#include <container/Container3D.h>
#include <constants/Axes.h>

#include <cmath>

namespace hist {
    /**
     * @brief This is a small wrapper around the Container3D class, indicating that the data
     *        is distributed along the constants::axes::d_vals axis.
     */
    class Distribution3D : public container::Container3D<constants::axes::d_type> {
        public:
            using Container3D::Container3D;

            void add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value);
            void add(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value);
    };
}