#pragma once

#include <vector>

namespace em {
    class ObjectBounds2D;
    class ObjectBounds3D {
        public: 
            ObjectBounds3D(unsigned int size_x, unsigned int size_y, unsigned int size_z);

            ~ObjectBounds3D();

            ObjectBounds2D& operator[](unsigned int z);

            const ObjectBounds2D& operator[](unsigned int z) const;

            unsigned int total_volume() const;

            unsigned int bounded_volume() const;

            unsigned int size() const;

        private:
            std::vector<ObjectBounds2D> bounds;
            unsigned int size_x, size_y, size_z;
    };
}