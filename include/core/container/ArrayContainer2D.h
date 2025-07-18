// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <array>

#if SAFE_MATH
    #include <utility/Exceptions.h>
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
            const typename std::array<T, M>::const_iterator begin(int i) const {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N)) {
                        throw except::out_of_bounds(
                            "ArrayContainer2D::begin: Index out of bounds "
                            "(" + std::to_string(N) + ") <= (" + std::to_string(i) + ")"
                        );
                    }
                #endif
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            const typename std::array<T, M>::const_iterator end(unsigned int i) const {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N)) {
                        throw except::out_of_bounds(
                            "ArrayContainer2D::end: Index out of bounds "
                            "(" + std::to_string(N) + ") <= (" + std::to_string(i) + ")"
                        );
                    }
                #endif
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::array<T, M>::iterator begin(unsigned int i) {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N)) {
                        throw except::out_of_bounds(
                            "ArrayContainer2D::begin: Index out of bounds "
                            "(" + std::to_string(N) + ") <= (" + std::to_string(i) + ")"
                        );
                    }
                #endif            
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::array<T, M>::iterator end(unsigned int i) {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N)) {
                        throw except::out_of_bounds(
                            "ArrayContainer2D::end: Index out of bounds "
                            "(" + std::to_string(N) + ") <= (" + std::to_string(i) + ")"
                        );
                    }
                #endif            
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::array<T, M>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::array<T, M>::const_iterator end() const {return data.end();}

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