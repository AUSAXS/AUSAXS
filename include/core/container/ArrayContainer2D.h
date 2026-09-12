// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <array>
#include <cassert>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::container {
    template <typename T, int N, int M>
    class ArrayContainer2D {
        public:
            constexpr ArrayContainer2D() noexcept = default;

            /**
             * @brief Get the value at index i, j, k. 
             */
            constexpr T& index(int i, int j) noexcept {return data[j + M*i];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            constexpr const T& index(int i, int j) const noexcept {return data[j + M*i];}

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::array<T, M>::const_iterator begin(int i) const {
                assert([&]() -> bool {
                    if (i < N) {return true;}
                    std::cout << "ArrayContainer2D::begin: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "ArrayContainer2D::begin: Index out of bounds.");
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::array<T, M>::const_iterator end(int i) const {
                assert([&]() -> bool {
                    if (i < N) {return true;}
                    std::cout << "ArrayContainer2D::end: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "ArrayContainer2D::end: Index out of bounds.");
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::array<T, M>::iterator begin(int i) {
                assert([&]() -> bool {
                    if (i < N) {return true;}
                    std::cout << "ArrayContainer2D::begin: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "ArrayContainer2D::begin: Index out of bounds.");
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::array<T, M>::iterator end(int i) {
                assert([&]() -> bool {
                    if (i < N) {return true;}
                    std::cout << "ArrayContainer2D::end: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "ArrayContainer2D::end: Index out of bounds.");
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, M>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, M>::const_iterator end() const {return data.end();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, M>::iterator begin() {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, M>::iterator end() {return data.end();}

        protected:
            std::array<T, N*M> data;
    };
}