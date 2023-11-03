#pragma once

#include <container/Container1D.h>
#include <constants/Constants.h>

#include <cmath>

namespace hist {
    class Distribution1D : public container::Container1D<constants::axes::d_type> {
        public:
            using Container1D::Container1D;

            /**
             * @brief Add a value for a given distance.
             */
            void add(float distance, constants::axes::d_type value);

            /**
             * @brief Add a value for a given distance.
             */
            void add(int32_t i, constants::axes::d_type value);
    };
}