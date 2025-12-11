// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Exceptions.h>
#include <math/indexers/Indexer1D.h>

#include <vector>

namespace ausaxs::container {
    /**
     * @brief Representation of a 1D container. 
     * 
     * This is just a convenience class supporting only indexing. 
     * Technically this is redundant since it is just a simpler std::vector, but it is used for consistency with the other two containers.
     */
    template <typename T>
    class Container1D : public utility::indexer::Indexer1D<Container1D<T>> {
        friend class utility::indexer::Indexer1D<Container1D<T>>;
        public:
            Container1D() : N(0), data(0) {}
            Container1D(int size) : N(size), data(size) {}
            Container1D(int size, const T& value) : N(size), data(size, value) {}
            Container1D(const std::vector<T>& data) : N(data.size()), data(data) {}
            Container1D(std::vector<T>&& data) : N(data.size()), data(std::move(data)) {}

            using utility::indexer::Indexer1D<Container1D<T>>::index;
            using utility::indexer::Indexer1D<Container1D<T>>::linear_index;
            T& operator()(int i) {return this->index(i);}
            const T& operator()(int i) const {return this->index(i);}
            T& operator[](int i) {return this->index(i);}
            const T& operator[](int i) const {return this->index(i);}

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