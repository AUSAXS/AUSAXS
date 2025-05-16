#pragma once

#include <utility/Exceptions.h>

#include <vector>

namespace ausaxs::container {
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
            T& operator()(int i, int j, int k) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || k >= static_cast<int>(L) || i < 0 || j < 0 || k < 0) {
                        throw except::out_of_bounds("Container3D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                    }
                #endif
                return data[k + L*(j + M*i)];
            }

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& operator()(int i, int j, int k) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || k >= static_cast<int>(L) || i < 0 || j < 0 || k < 0) {
                        throw except::out_of_bounds("Container3D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ", " + std::to_string(L) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ")");
                    }
                #endif
                return data[k + L*(j + M*i)];
            }

            /**
             * @brief Get the value at index i, j, k. 
             */
            T& index(int i, int j, int k) {return operator()(i, j, k);}

            /**
             * @brief Get the value at index i, j, k. 
             */
            const T& index(int i, int j, int k) const {return operator()(i, j, k);}

            /**
             * @brief Get the vector sum of all members of the container.
             */
            std::vector<T> vector_sum() const {
                std::vector<T> sum(L, 0);
                for (int i = 0; i < data.size(); ++i) {sum[i%L] += data[i];}
                return sum;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            const typename std::vector<T>::const_iterator begin(int i, int j) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || i < 0 || j < 0) {
                        throw except::out_of_bounds("Container3D::begin: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                    }
                #endif
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            const typename std::vector<T>::const_iterator end(int i, int j) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || i < 0 || j < 0) {
                        throw except::out_of_bounds("Container3D::end: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                    }
                #endif
                return data.begin() + L*(j + M*i) + L;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            typename std::vector<T>::iterator begin(int i, int j) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || i < 0 || j < 0) {
                        throw except::out_of_bounds("Container3D::begin: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                    }
                #endif
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::iterator end(int i, int j) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || j >= static_cast<int>(M) || i < 0 || j < 0) {
                        throw except::out_of_bounds("Container3D::end: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
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

            /**
             * @brief Get the number of elements in the x direction.
             */
            std::size_t size_x() const {return N;}

            /**
             * @brief Get the number of elements in the y direction.
             */
            std::size_t size_y() const {return M;}

            /**
             * @brief Get the length of each (x, y) element.
             */
            std::size_t size_z() const {return L;}

            /**
             * @brief Resize the container to contain @a size elements for each (x, y) index.
             */
            void resize(int size) {
                Container3D tmp(N, M, size);
                for (int i = 0; i < static_cast<int>(N); i++) {
                    for (int j = 0; j < static_cast<int>(M); j++) {
                        std::move(begin(i, j), begin(i, j)+std::min<int>(size, L), tmp.begin(i, j));
                    }
                }
                L = size;
                data = std::move(tmp.data);
            }

            /**
             * @brief Check if the container is empty.
             */
            bool empty() const {return data.empty();}

        protected:
            std::size_t N, M, L;
            std::vector<T> data;
    };
}