#pragma once

#include <array>

namespace container {
    template <typename T, int size_x, int size_y>
    class ArrayContainer2D {
        public:
            constexpr ArrayContainer2D() noexcept : N(size_x), M(size_y) {}

            /**
             * @brief Get the value at index i, j, k. 
             */
            constexpr T& index(unsigned int i, unsigned int j, unsigned int k) noexcept {return data[j + M*i];}

            /**
             * @brief Get the value at index i, j, k. 
             */
            constexpr const T& index(unsigned int i, unsigned int j, unsigned int k) const noexcept {return data[j + M*i];}

             /**
             * @brief Get the value at index i, j. 
             */
            constexpr T& index(unsigned int i, unsigned int j) noexcept {return data[j + M*i];}

            /**
             * @brief Get the value at index i, j. 
             */
            constexpr const T& index(unsigned int i, unsigned int j) const noexcept {return data[j + M*i];}

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            const typename std::array<T, size_y>::const_iterator begin(unsigned int i) const {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            const typename std::array<T, size_y>::const_iterator end(unsigned int i) const {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::array<T, size_y>::iterator begin(unsigned int i) {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif            
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::array<T, size_y>::iterator end(unsigned int i) {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif            
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::array<T, size_y>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::array<T, size_y>::const_iterator end() const {return data.end();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, size_y>::iterator begin() {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::array<T, size_y>::iterator end() {return data.end();}

            unsigned int N, M;

        protected:
            std::array<T, size_x*size_y> data;
    };
}