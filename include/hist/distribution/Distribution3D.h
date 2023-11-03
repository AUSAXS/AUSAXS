#pragma once

#include <container/Container3D.h>
#include <constants/Constants.h>

#include <cmath>

namespace hist {
    class Distribution3D : public container::Container3D<constants::axes::d_type> {
        public:
            using Container3D::Container3D;

            void add(unsigned int x, unsigned int y, float distance, constants::axes::d_type value);
            void add(unsigned int x, unsigned int y, int32_t i, constants::axes::d_type value);
    };
}