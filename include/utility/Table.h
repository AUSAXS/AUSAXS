#pragma once

#include <vector>
#include <stdexcept>

namespace table {
    /**
     * @brief A representation of a table. It uses only a single contiguous vector as storage for better locality. 
     */
    class Table {
        public: 
            /**
             * @brief Default constructor.
             */
            Table() : N(0), M(0), data(0) {}

            /**
             * @brief Construct a new Table of a given size.
             * 
             * @param N The number of rows to reserve space for. 
             * @param M The number of columns to reserve space for. 
             */
            Table(unsigned int N, unsigned int M) : N(N), M(M), data(N*M) {}

            /**
             * @brief Read-only access to the given table index.  
             * 
             * @param i The row index.
             * @param j The column index. 
             */
            const double& index(unsigned int i, unsigned int j) const {
                #ifdef DEBUG
                    if (i >= N || j >= M) {
                        throw std::out_of_range("Table::index: Table index out of bounds. \nAccessed (" + std::to_string(i) + ", " + std::to_string(j) + ") in table of size (" + std::to_string(N) + ", " + std::to_string(M) + ")");
                    }
                #endif
                return data[M*i + j];
            }

            /**
             * @brief Read-write access to the given table index.  
             * 
             * @param i The row index.
             * @param j The column index. 
             */
            double& index(unsigned int i, unsigned int j) {
                #ifdef DEBUG
                    if (i >= N || j >= M) {
                        throw std::out_of_range("Table::index: Table index out of bounds. \nAccessed (" + std::to_string(i) + ", " + std::to_string(j) + ") in table of size (" + std::to_string(N) + ", " + std::to_string(M) + ")");
                    }
                #endif
                return data[M*i + j];
            }

            unsigned int N, M;
            std::vector<double> data;
    };
}