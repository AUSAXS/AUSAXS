// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <array>
#include <cassert>
#include <vector>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::container {
    template <typename T, int N, int M, int L>
    class ArrayContainer3D {
        public:
            constexpr ArrayContainer3D() noexcept = default;

            /**
             * @brief Get the value at index i, j, k. 
             */
            T& operator()(int i, int j, int k) {return data[k + L*(j + M*i)];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& operator()(int i, int j, int k) const {return data[k + L*(j + M*i)];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            T& index(int i, int j, int k) {return operator()(i, j, k);}

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& index(int i, int j, int k) const {return operator()(i, j, k);}

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            typename std::vector<T>::const_iterator begin(int i, int j) const {
                assert([&]() -> bool {
                    if (i < N && j < M) {return true;}
                    std::cout << "ArrayContainer3D::begin: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "ArrayContainer3D::begin: Index out of bounds.");
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::const_iterator end(int i, int j) const {
                assert([&]() -> bool {
                    if (i < N && j < M) {return true;}
                    std::cout << "ArrayContainer3D::end: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "ArrayContainer3D::end: Index out of bounds.");
                return data.begin() + L*(j + M*i) + L;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            typename std::vector<T>::iterator begin(int i, int j) {
                assert([&]() -> bool {
                    if (i < N && j < M) {return true;}
                    std::cout << "ArrayContainer3D::begin: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "ArrayContainer3D::begin: Index out of bounds.");
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::iterator end(int i, int j) {
                assert([&]() -> bool {
                    if (i < N && j < M) {return true;}
                    std::cout << "ArrayContainer3D::end: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "ArrayContainer3D::end: Index out of bounds.");
                return data.begin() + L*(j + M*i) + L;
            }

            /**
             * @brief Get an iterator to the beginning of the entire container. 
             */
            typename std::vector<T>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the end of the entire container. 
             */
            typename std::vector<T>::const_iterator end() const {return data.end();}

            /**
             * @brief Get an iterator to the beginning of the entire container. 
             */
            typename std::vector<T>::iterator begin() {return data.begin();}

            /**
             * @brief Get an iterator to the end of the entire container. 
             */
            typename std::vector<T>::iterator end() {return data.end();}

        protected:
            std::array<T, N*M*L> data;
    };
}