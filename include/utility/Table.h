#pragma once

#include <vector>
#include <unordered_map>
#include <iostream>

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
             * @brief Constructor.
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
            const double& index(unsigned int i, unsigned int j) const {return data.at(M*i + j);}

            /**
             * @brief Read-write access to the given table index.  
             * 
             * @param i The row index.
             * @param j The column index. 
             */
            double& index(unsigned int i, unsigned int j) {return data.at(M*i + j);}

            unsigned int N, M;
            std::vector<double> data;
    };
}