#pragma once

#include <utility/Exceptions.h>
#include <Symbols.h>

#include <vector>

namespace container {
    /**
     * @brief Representation of a dense 2D container. 
     * 
     * This is just a convenience class supporting only indexing.
     */
    template <typename T>
    class Container2D {
        public:
            Container2D() : N(0), M(0), data(0) {}
            Container2D(unsigned int width, unsigned int height) : N(width), M(height), data(width*height) {}
            Container2D(unsigned int width, unsigned int height, const T& value) : N(width), M(height), data(width*height, value) {}

            /**
             * @brief Get the value at index i, j. 
             */
            T& operator()(unsigned int i, unsigned int j) {
                #if (SAFE_MATH)
                    if (i >= N || j >= M) {
                        throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                    }
                #endif
                return data[j + M*i];
            }

            /**
             * @brief Get the value at index i, j. 
             */
            const T& operator()(unsigned int i, unsigned int j) const {
                #if (SAFE_MATH)
                    if (i >= N || j >= M) {
                        throw except::out_of_bounds("Container2D::operator: Index out of bounds (" + std::to_string(N) + ", " + std::to_string(M) + ") <= (" + std::to_string(i) + ", " + std::to_string(j) + ")");
                    }
                #endif
                return data[j + M*i];
            }

            /**
             * @brief Get the value at index i, j. 
             */
            T& index(unsigned int i, unsigned int j) {return operator()(i, j);}

            /**
             * @brief Get the value at index i, j. 
             */
            const T& index(unsigned int i, unsigned int j) const {return operator()(i, j);}

            /**
             * @brief Get an iterator to the beginning of the vector at index i.
             */
            const typename std::vector<T>::const_iterator begin(unsigned int i) const {
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
            const typename std::vector<T>::const_iterator end(unsigned int i) const {
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
            typename std::vector<T>::iterator begin(unsigned int i) {
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
            typename std::vector<T>::iterator end(unsigned int i) {
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
            unsigned int size_x() const {return N;}

            /**
             * @brief Get the length of each x-element.
             */
            unsigned int size_y() const {return M;}

            /**
             * @brief Resize the container to contain @a size elements for each x index.
             */
            void resize(unsigned int size) {
                Container2D tmp(N, size);
                for (unsigned int i = 0; i < N; i++) {
                    std::move(begin(i), begin(i)+size, tmp.begin(i));
                }
                M = size;
                data = std::move(tmp.data);                
            }

        protected:
            unsigned int N, M;
            std::vector<T> data;
    };
}