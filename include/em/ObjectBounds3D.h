#pragma once

#include <em/EMFwd.h>

#include <vector>

namespace em {
    class ObjectBounds3D {
        public: 
            ObjectBounds3D(unsigned int size_x, unsigned int size_y, unsigned int size_z);

            ~ObjectBounds3D();

            [[nodiscard]] ObjectBounds2D& operator[](unsigned int z);

            [[nodiscard]] const ObjectBounds2D& operator[](unsigned int z) const;

            [[nodiscard]] unsigned int total_volume() const;

            [[nodiscard]] unsigned int bounded_volume() const;

            [[nodiscard]] unsigned int size_x() const;

            [[nodiscard]] unsigned int size_y() const;

            [[nodiscard]] unsigned int size_z() const;

        private:
            std::vector<ObjectBounds2D> bounds;
            unsigned int _size_x, _size_y, _size_z;
    };
}