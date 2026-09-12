// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/indexers/Indexer3D.h>

#include <cassert>
#include <vector>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::container {
    /**
     * @brief Representation of a dense 3D container. 
     * 
     * This is just a convenience class supporting only basic indexing.
     */
    template <typename T>
    class Container3D : utility::indexer::Indexer3D<Container3D<T>> {
        friend class utility::indexer::Indexer3D<Container3D<T>>;
        public:
            Container3D() : N(0), M(0), L(0), data(0) {}
            Container3D(int width, int height, int depth) : N(width), M(height), L(depth), data(width * height * depth) {}
            Container3D(int width, int height, int depth, const T& value) : N(width), M(height), L(depth), data(width * height * depth, value) {}

            using utility::indexer::Indexer3D<Container3D<T>>::index;
            using utility::indexer::Indexer3D<Container3D<T>>::linear_index;
            T& operator()(int i, int j, int k) {return this->index(i, j, k);}
            const T& operator()(int i, int j, int k) const {return this->index(i, j, k);}

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
            typename std::vector<T>::const_iterator begin(int i, int j) const {
                assert([&]() -> bool {
                    if (0 <= i && i < N && 0 <= j && j < M) {return true;}
                    std::cout << "Container3D::begin: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "Container3D::begin: Index out of bounds.");
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::const_iterator end(int i, int j) const {
                assert([&]() -> bool {
                    if (0 <= i && i < N && 0 <= j && j < M) {return true;}
                    std::cout << "Container3D::end: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "Container3D::end: Index out of bounds.");
                return data.begin() + L*(j + M*i) + L;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            typename std::vector<T>::iterator begin(int i, int j) {
                assert([&]() -> bool {
                    if (0 <= i && i < N && 0 <= j && j < M) {return true;}
                    std::cout << "Container3D::begin: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "Container3D::begin: Index out of bounds.");
                return data.begin() + L*(j + M*i);
            }

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            typename std::vector<T>::iterator end(int i, int j) {
                assert([&]() -> bool {
                    if (0 <= i && i < N && 0 <= j && j < M) {return true;}
                    std::cout << "Container3D::end: Index out of bounds (" << N << ", " << M << ") <= (" << i << ", " << j << ")" << std::endl;
                    return false;
                }() && "Container3D::end: Index out of bounds.");
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

            /**
             * @brief Get the number of elements in the x direction.
             */
            int size_x() const {return N;}

            /**
             * @brief Get the number of elements in the y direction.
             */
            int size_y() const {return M;}

            /**
             * @brief Get the length of each (x, y) element.
             */
            int size_z() const {return L;}

            /**
             * @brief Resize the container to contain @a size elements for each (x, y) index.
             */
            void resize(int size) {
                Container3D tmp(N, M, size);
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < M; j++) {
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
            int N, M, L;
            std::vector<T> data;
    };
}