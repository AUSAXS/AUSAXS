#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

/**
 * @brief Representation of a dense 3D container. 
 * 
 * This is just a convenience class supporting only basic indexing.
 */
template <typename T>
class Container3D {
    public:
        Container3D() : N(0), M(0), L(0), data(0) {}
        Container3D(unsigned int width, unsigned int height, unsigned int depth) : N(width), M(height), L(depth), data(width * height * depth) {}
        Container3D(unsigned int width, unsigned int height, unsigned int depth, const T& value) : N(width), M(height), L(depth), data(width * height * depth, value) {}

        /**
         * @brief Get the value at index i, j, k. 
         */
        T& operator()(unsigned int i, unsigned int j, unsigned int k) {
            #if (SAFE_MATH)
                if (i >= N || j >= M || k >= L) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        /**
         * @brief Get the value at index i, j, k. 
         */
        const T& operator()(unsigned int i, unsigned int j, unsigned int k) const {
            #if (SAFE_MATH)
                if (i >= N || j >= M || k >= L) {
                    throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                }
            #endif
            return data[k + L*(j + M*i)];
        }

        /**
         * @brief Get the value at index i, j, k. 
         */
        T& index(unsigned int i, unsigned int j, unsigned int k) {return operator()(i, j, k);}

        /**
         * @brief Get the value at index i, j, k. 
         */
        const T& index(unsigned int i, unsigned int j, unsigned int k) const {return operator()(i, j, k);}

        /**
         * @brief Get the vector sum of all members of the container.
         */
        std::vector<T> vector_sum() const {
            std::vector<T> sum(L, 0);
            for (unsigned int i = 0; i < data.size(); ++i) {sum[i%L] += data[i];}
            return sum;
        }

        /**
         * @brief Get an iterator to the beginning of the vector at index i, j.
         */
        const typename std::vector<T>::const_iterator begin(unsigned int i, unsigned int j) const {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data.begin() + L*(j + M*i);
        }

        /**
         * @brief Get an iterator to the end of the vector at index i, j.
         */
        const typename std::vector<T>::const_iterator end(unsigned int i, unsigned int j) const {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data.begin() + L*(j + M*i) + L;
        }

        /**
         * @brief Get an iterator to the beginning of the vector at index i, j.
         */
        typename std::vector<T>::iterator begin(unsigned int i, unsigned int j) {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data.begin() + L*(j + M*i);
        }

        /**
         * @brief Get an iterator to the end of the vector at index i, j.
         */
        typename std::vector<T>::iterator end(unsigned int i, unsigned int j) {
            #if (SAFE_MATH)
                if (i >= N || j >= M) {
                    throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                }
            #endif
            return data.begin() + L*(j + M*i) + L;
        }

        /**
         * @brief Get an iterator to the beginning of the entire container. 
         */
        const typename std::vector<T>::const_iterator begin() const {return data.begin();}

        /**
         * @brief Get an iterator to the end of the entire container. 
         */
        const typename std::vector<T>::const_iterator end() const {return data.end();}

        /**
         * @brief Get an iterator to the beginning of the entire container. 
         */
        typename std::vector<T>::iterator begin() {return data.begin();}

        /**
         * @brief Get an iterator to the end of the entire container. 
         */
        typename std::vector<T>::iterator end() {return data.end();}

        unsigned int N, M, L;

    protected:
        std::vector<T> data;
};
