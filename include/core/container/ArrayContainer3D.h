#pragma once

#include <array>
#include <vector>

#if SAFE_MATH
    #include <utility/Exceptions.h>
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
            const typename std::vector<T>::const_iterator begin(int i, int j) const {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M)) {
                        throw except::out_of_bounds(
                            "ArrayContainer3D::begin: Index out of bounds "
                            "(" + std::to_string(N) + ", " + std::to_string(M) + ") <= "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
                    }
                #endif
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            const typename std::vector<T>::const_iterator end(int i, int j) const {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M)) {
                        throw except::out_of_bounds(
                            "ArrayContainer3D::end: Index out of bounds "
                            "(" + std::to_string(N) + ", " + std::to_string(M) + ") <= "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
                    }
                #endif
                return data.begin() + L*(j + M*i) + L;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            typename std::vector<T>::iterator begin(int i, int j) {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M)) {
                        throw except::out_of_bounds(
                            "ArrayContainer3D::begin: Index out of bounds "
                            "(" + std::to_string(N) + ", " + std::to_string(M) + ") <= "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
                    }
                #endif
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::iterator end(int i, int j) {
                #if SAFE_MATH
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M)) {
                        throw except::out_of_bounds(
                            "ArrayContainer3D::end: Index out of bounds "
                            "(" + std::to_string(N) + ", " + std::to_string(M) + ") <= "
                            "(" + std::to_string(i) + ", " + std::to_string(j) + ")"
                        );
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

        protected:
            std::array<T, N*M*L> data;
    };
}