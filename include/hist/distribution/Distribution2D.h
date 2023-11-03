#pragma once

#include <container/Container2D.h>
#include <constants/Constants.h>

namespace hist {
    class Distribution2D : public container::Container2D<constants::axes::d_type> {
        public:
            using Container2D::Container2D;

            /**
             * @brief Add a value for a given distance.
             */
            void add(unsigned int x, float distance, constants::axes::d_type value);
            void add(unsigned int x, int32_t i, constants::axes::d_type value);
    };
}