#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

namespace container {
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
            Container1D(unsigned int size) : N(size), data(size) {}
            Container1D(unsigned int size, const T& value) : N(size), data(size, value) {}

            T& operator()(unsigned int i) {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data[i];
            }

            const T& operator()(unsigned int i) const {
                #if (SAFE_MATH)
                    if (i >= N) {
                        throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ") <= (" + std::to_string(i) + ")");
                    }
                #endif
                return data[i];
            }

            T& index(unsigned int i) {return operator()(i);}
            const T& index(unsigned int i) const {return operator()(i);}

            const typename std::vector<T>::const_iterator begin() const {return data.begin();}
            const typename std::vector<T>::const_iterator end() const {return data.end();}

            typename std::vector<T>::iterator begin() {return data.begin();}
            typename std::vector<T>::iterator end() {return data.end();}

            /**
             * @brief Get the size of the container.
             */
            unsigned int size() const {return N;}

            /**
             * @brief Resize the container to contain @a size elements.
             */
            void resize(unsigned int size) {
                N = size;
                data.resize(size);
            }

        protected:
            unsigned int N;
            std::vector<T> data;
    };
}