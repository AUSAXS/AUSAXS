#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

/**
 * @brief Representation of a 3D container. 
 * 
 * This is just a convenience class supporting only indexing.
 */
template <typename T>
class Container3D {
    public:
        Container3D() : N(0), M(0), L(0), data(0) {}
        Container3D(unsigned int width, unsigned int height, unsigned int depth) : N(width), M(height), L(depth), data(width * height * depth) {}

        T& operator()(unsigned int i, unsigned int j, unsigned int k) {
            #if (SAFE_MATH)
                if (i < 0 || i >= N || j < 0 || j >= M || k < 0 || k >= L) {
                    throw except::out_of_bounds("Container3D::operator: Index out of bounds");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        const T& operator()(unsigned int i, unsigned int y, unsigned int z) const {
            #if (SAFE_MATH)
                if (i < 0 || i >= N || j < 0 || j >= M || k < 0 || k >= L) {
                    throw except::out_of_bounds("Container3D::operator: Index out of bounds");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        unsigned int N, M, L;

    private:
        std::vector<T> data;
};
