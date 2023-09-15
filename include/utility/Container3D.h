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
        Container3D(unsigned int width, unsigned int height, unsigned int depth, const T& value) : N(width), M(height), L(depth), data(width * height * depth, value) {}

        T& operator()(unsigned int i, unsigned int j, unsigned int k) {
            #if (SAFE_MATH)
                if (i >= N || j >= M || k >= L) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        const T& operator()(unsigned int i, unsigned int j, unsigned int k) const {
            #if (SAFE_MATH)
                if (i >= N || j >= M || k >= L) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        T& index(unsigned int i, unsigned int j, unsigned int k) {return operator()(i, j, k);}
        const T& index(unsigned int i, unsigned int j, unsigned int k) const {return operator()(i, j, k);}

        const typename std::vector<T>::const_iterator begin() const {return data.begin();}
        const typename std::vector<T>::const_iterator end() const {return data.begin();}

        typename std::vector<T>::iterator begin() {return data.begin();}
        typename std::vector<T>::iterator end() {return data.begin();}

        unsigned int N, M, L;

    private:
        std::vector<T> data;
};
