// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/indexers/Indexer2D.h>

#include <cassert>
#include <span>
#include <vector>

#ifndef NDEBUG
    #include <iostream>  // only the asserts below print
#endif

namespace ausaxs::container {
    using utility::indexer::Shape;

    /**
     * @brief Representation of a dense 2D container. 
     * 
     * This is just a convenience class supporting only indexing. With Shape::Triangular, both dimensions must be equal, and only
     * one element per unordered pair (i, j) is stored; (i, j) and (j, i) then name the same element. Its rows have no fixed
     * length, so only the square layout has the row accessors.
     */
    template <typename T, Shape S = Shape::Square>
    class Container2D : utility::indexer::Indexer2D<Container2D<T, S>, S> {
        using Indexer = utility::indexer::Indexer2D<Container2D<T, S>, S>;
        friend Indexer;
        public:
            using value_type = T;

            Container2D() : N(0), M(0), data(0) {}
            Container2D(int width, int height) : N(width), M(height), data(pair_count()) {}
            Container2D(int width, int height, const T& value) : N(width), M(height), data(pair_count(), value) {}

            using Indexer::index;
            using Indexer::linear_index;
            T& operator()(int i, int j) {return this->index(i, j);}
            const T& operator()(int i, int j) const {return this->index(i, j);}

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::vector<T>::const_iterator begin(int i) const requires (S == Shape::Square) {
                assert([&]() -> bool {
                    if (0 <= i && i < N) {return true;}
                    std::cout << "Container2D::begin: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "Container2D::begin: Index out of bounds.");
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::vector<T>::const_iterator end(int i) const requires (S == Shape::Square) {
                assert([&]() -> bool {
                    if (0 <= i && i < N) {return true;}
                    std::cout << "Container2D::end: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "Container2D::end: Index out of bounds.");
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::vector<T>::iterator begin(int i) requires (S == Shape::Square) {
                assert([&]() -> bool {
                    if (0 <= i && i < N) {return true;}
                    std::cout << "Container2D::begin: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "Container2D::begin: Index out of bounds.");
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::vector<T>::iterator end(int i) requires (S == Shape::Square) {
                assert([&]() -> bool {
                    if (0 <= i && i < N) {return true;}
                    std::cout << "Container2D::end: Index out of bounds (" << N << ") <= (" << i << ")" << std::endl;
                    return false;
                }() && "Container2D::end: Index out of bounds.");
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get the vector at index i.
             */
            std::span<T> row(int i) requires (S == Shape::Square) {return {begin(i), static_cast<std::size_t>(M)};}

            /**
             * @brief Get the vector at index i.
             */
            std::span<const T> row(int i) const requires (S == Shape::Square) {return {begin(i), static_cast<std::size_t>(M)};}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::vector<T>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::vector<T>::const_iterator end() const {return data.end();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::vector<T>::iterator begin() {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            typename std::vector<T>::iterator end() {return data.end();}

            /**
             * @brief Get the number of contained x-elements.
             */
            int size_x() const {return N;}

            /**
             * @brief Get the length of each x-element.
             */
            int size_y() const {return M;}

            /**
             * @brief Resize the container to contain @a size elements for each x index.
             */
            void resize(int size) requires (S == Shape::Square) {
                Container2D tmp(N, size);
                for (int i = 0; i < N; i++) {
                    std::move(begin(i), begin(i)+std::min<int>(size, M), tmp.begin(i));
                }
                M = size;
                data = std::move(tmp.data);                
            }

            /**
             * @brief Check if the container is empty.
             */
            bool empty() const {return data.empty();}

        protected:
            using Indexer::pair_count;
            int N, M;
            std::vector<T> data;
    };
}