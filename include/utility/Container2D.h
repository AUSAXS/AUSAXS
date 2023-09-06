#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

/**
 * @brief Representation of a 2D container. 
 * 
 * This is just a convenience class supporting only indexing.
 */
template <typename T>
class Container2D {
    public:
        Container2D() : N(0), M(0), data(0) {}
        Container2D(unsigned int width, unsigned int height) : N(width), M(height), data(width*height) {}

        T& operator()(unsigned int i, unsigned int j) {
            #if (SAFE_MATH)
                if (i < 0 || i >= N || j < 0 || j >= M) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds");
                }
            #endif
            return data[j + M*i];
        }

        const T& operator()(unsigned int i, unsigned int y) const {
            #if (SAFE_MATH)
                if (i < 0 || i >= N || j < 0 || j >= M) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds");
                }
            #endif
            return data[j + M*i];
        }

        T& index(unsigned int i, unsigned int j) {return operator()(i, j);}
        const T& index(unsigned int i, unsigned int j) const {return operator()(i, j);}

        unsigned int N, M;

    private:
        std::vector<T> data;
};
