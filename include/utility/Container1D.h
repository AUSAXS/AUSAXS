#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

/**
 * @brief Representation of a 1D container. 
 * 
 * This is just a convenience class supporting only indexing.
 */
template <typename T>
class Container1D {
    public:
        Container1D() : N(0), data(0) {}
        Container1D(unsigned int size) : N(size), data(size) {}

        T& operator()(unsigned int i) {
            #if (SAFE_MATH)
                if (i < 0 || i >= N) {
                    throw except::out_of_bounds("Container1D::operator: Index out of bounds");
                }
            #endif
            return data[i];
        }

        const T& operator()(unsigned int i) const {
            #if (SAFE_MATH)
                if (i < 0 || i >= N) {
                    throw except::out_of_bounds("Container1D::operator: Index out of bounds");
                }
            #endif
            return data[i];
        }

        T& index(unsigned int i) {return operator()(i);}
        const T& index(unsigned int i) const {return operator()(i);}

        unsigned int N;

    private:
        std::vector<T> data;
};
