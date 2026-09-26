// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <math/indexers/Indexer3D.h>

#include <algorithm>
#include <vector>

namespace ausaxs::container {
    using utility::indexer::Shape;

    /**
     * @brief Representation of a dense 3D container. 
     * 
     * This is just a convenience class supporting only basic indexing. With Shape::Triangular, the first two dimensions must be
     * equal, and only one row per unordered pair (i, j) is stored; (i, j) and (j, i) then name the same row.
     */
    template <typename T, Shape S = Shape::Square>
    class Container3D : utility::indexer::Indexer3D<Container3D<T, S>, S> {
        using Indexer = utility::indexer::Indexer3D<Container3D<T, S>, S>;
        friend Indexer;
        public:
            using value_type = T;
            using Indexer::row;
            using Indexer::rows;
            using Indexer::index;
            using Indexer::linear_index;

            Container3D() : N(0), M(0), L(0), data(0) {}
            Container3D(int width, int height, int depth) : N(width), M(height), L(depth), data(pair_count() * depth) {}
            Container3D(int width, int height, int depth, const T& value) : N(width), M(height), L(depth), data(pair_count() * depth, value) {}

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
            auto begin(int i, int j) const {return row(i, j).begin();}

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            auto end(int i, int j) const {return row(i, j).end();}

            /**
             * @brief Get an iterator to the beginning of the vector at index i, j.
             */
            auto begin(int i, int j) {return row(i, j).begin();}

            /**
             * @brief Get an iterator to the end of the vector at index i, j.
             */
            auto end(int i, int j) {return row(i, j).end();}

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
                auto to = tmp.rows().begin();
                for (auto from : rows()) {
                    std::move(from.begin(), from.begin() + std::min(size, L), (*to++).begin());
                }
                L = size;
                data = std::move(tmp.data);
            }

            /**
             * @brief Check if the container is empty.
             */
            bool empty() const {return data.empty();}

        protected:
            using Indexer::pair_count;
            int N, M, L;
            std::vector<T> data;
    };
}