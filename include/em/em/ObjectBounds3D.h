// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/EMFwd.h>

#include <vector>

namespace ausaxs::em {
    class ObjectBounds3D {
        public: 
            ObjectBounds3D(int size_x, int size_y, int size_z);

            ~ObjectBounds3D();

            [[nodiscard]] ObjectBounds2D& operator[](int z);

            [[nodiscard]] const ObjectBounds2D& operator[](int z) const;

            [[nodiscard]] int total_volume() const;

            [[nodiscard]] int bounded_volume() const;

            [[nodiscard]] int size_x() const;

            [[nodiscard]] int size_y() const;

            [[nodiscard]] int size_z() const;

        private:
            std::vector<ObjectBounds2D> bounds;
            int _size_x, _size_y, _size_z;
    };
}