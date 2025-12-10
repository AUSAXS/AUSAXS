// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Exceptions.h>
#include <math/indexers/Indexer2D.h>

#include <vector>

namespace ausaxs::container {
    /**
     * @brief Representation of a dense 2D container. 
     * 
     * This is just a convenience class supporting only indexing.
     */
    template <typename T>
    class Container2D : utility::indexer::Indexer2D<Container2D<T>> {
        friend class utility::indexer::Indexer2D<Container2D<T>>;
        public:
            Container2D() : N(0), M(0), data(0) {}
            Container2D(int width, int height) : N(width), M(height), data(width*height) {}
            Container2D(int width, int height, const T& value) : N(width), M(height), data(width*height, value) {}

            using utility::indexer::Indexer2D<Container2D<T>>::index;
            using utility::indexer::Indexer2D<Container2D<T>>::linear_index;
            T& operator()(int i, int j) {return this->index(i, j);}
            const T& operator()(int i, int j) const {return this->index(i, j);}

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            const typename std::vector<T>::const_iterator begin(int i) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            const typename std::vector<T>::const_iterator end(int i) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            typename std::vector<T>::iterator begin(int i) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container2D::begin: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif            
                return data.begin() + i*M;
            }

            /**
             * @brief Get an iterator to the end of the vector at index i.
             */
            typename std::vector<T>::iterator end(int i) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container2D::end: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif            
                return data.begin() + i*M + M;
            }

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::vector<T>::const_iterator begin() const {return data.begin();}

            /**
             * @brief Get an iterator to the beginning of the entire container.
             */
            const typename std::vector<T>::const_iterator end() const {return data.end();}

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
            std::size_t size_x() const {return N;}

            /**
             * @brief Get the length of each x-element.
             */
            std::size_t size_y() const {return M;}

            /**
             * @brief Resize the container to contain @a size elements for each x index.
             */
            void resize(int size) {
                Container2D tmp(N, size);
                for (int i = 0; i < static_cast<int>(N); i++) {
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
            std::size_t N, M;
            std::vector<T> data;
    };
}