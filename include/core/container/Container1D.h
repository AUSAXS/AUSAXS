// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Exceptions.h>

#include <vector>

namespace ausaxs::container {
    /**
     * @brief Representation of a 1D container. 
     * 
     * This is just a convenience class supporting only indexing. 
     * Technically this is redundant since it is just a simpler std::vector, but it is used for consistency with the other two containers.
     */
    template <typename T>
    class Container1D {
        public:
            Container1D() : N(0), data(0) {}
            Container1D(int size) : N(size), data(size) {}
            Container1D(int size, const T& value) : N(size), data(size, value) {}
            Container1D(const std::vector<T>& data) : N(data.size()), data(data) {}
            Container1D(std::vector<T>&& data) : N(data.size()), data(std::move(data)) {}

            T& operator()(int i) {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container1D::operator: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data[i];
            }

            const T& operator()(int i) const {
                #if (SAFE_MATH)
                    if (i >= static_cast<int>(N) || i < 0) {
                        throw except::out_of_bounds("Container1D::operator: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data[i];
            }

            T& operator[](int i) {return operator()(i);}
            const T& operator[](int i) const {return operator()(i);}

            T& index(int i) {return operator()(i);}
            const T& index(int i) const {return operator()(i);}

            const typename std::vector<T>::const_iterator begin() const {return data.begin();}
            const typename std::vector<T>::const_iterator end() const {return data.end();}

            typename std::vector<T>::iterator begin() {return data.begin();}
            typename std::vector<T>::iterator end() {return data.end();}

            /**
             * @brief Get the size of the container.
             */
            std::size_t size() const {return N;}

            /**
             * @brief Resize the container to contain @a size elements.
             */
            void resize(int size) {
                N = size;
                data.resize(size);
            }

            operator std::vector<T>&() {return data;}

            std::vector<T>& get_data() {return data;}

            /**
             * @brief Check if the container is empty.
             */
            bool empty() const {return data.empty();}

        protected:
            int N;
            std::vector<T> data;
    };
}