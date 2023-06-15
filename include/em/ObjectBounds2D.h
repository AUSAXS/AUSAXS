#pragma once

#include <vector>

class Limit;
namespace em {
    /**
     * @brief Describes the bounds of some object contained within a 2D matrix. 
     */
    class ObjectBounds2D {
        public:
            ObjectBounds2D(unsigned int size_x, unsigned int size_y);

            ~ObjectBounds2D();

            Limit& operator[](unsigned int x);

            const Limit& operator[](unsigned int x) const;

            unsigned int size() const;

            bool empty() const;

            unsigned int bounded_area() const;

            unsigned int total_area() const;

            bool operator==(const ObjectBounds2D& other) const;

        private:
            std::vector<Limit> bounds;
            unsigned int N, M;
    };
}