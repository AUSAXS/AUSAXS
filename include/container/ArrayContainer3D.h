#pragma once

#include <array>

namespace container {
    template <typename T, int size_x, int size_y, int size_z>
    class ArrayContainer3D {
        public:
            ArrayContainer3D() : N(size_x), M(size_y), L(size_z) {}

            /**
             * @brief Get the value at index i, j, k. 
             */
            T& operator()(unsigned int i, unsigned int j, unsigned int k) {return data[k + L*(j + M*i)];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& operator()(unsigned int i, unsigned int j, unsigned int k) const {return data[k + L*(j + M*i)];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            T& index(unsigned int i, unsigned int j, unsigned int k) {return operator()(i, j, k);}

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& index(unsigned int i, unsigned int j, unsigned int k) const {return operator()(i, j, k);}

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
            std::array<T, size_x*size_y*size_z> data;
    };
}