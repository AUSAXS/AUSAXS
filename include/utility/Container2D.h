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
        Container2D(unsigned int width, unsigned int height, const T& value) : N(width), M(height), data(width*height, value) {}

        T& operator()(unsigned int i, unsigned int j) {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data[j + M*i];
        }

        const T& operator()(unsigned int i, unsigned int j) const {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data[j + M*i];
        }

        T& index(unsigned int i, unsigned int j) {return operator()(i, j);}
        const T& index(unsigned int i, unsigned int j) const {return operator()(i, j);}

        const typename std::vector<T>::const_iterator begin() const {return data.begin();}
        const typename std::vector<T>::const_iterator end() const {return data.end();}

        typename std::vector<T>::iterator begin() {return data.begin();}
        typename std::vector<T>::iterator end() {return data.end();}

        unsigned int N, M;

    private:
        std::vector<T> data;
};
