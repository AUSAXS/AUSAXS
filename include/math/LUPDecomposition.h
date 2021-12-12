#pragma once

// forwards declaration
class Matrix;

#include "Decomposition.h"
#include "Matrix.h"
#include "Vector.h"

#include <stdexcept>

class LUPDecomposition : public Decomposition {
public: 
    LUPDecomposition(const Matrix& A) : A(A.copy()) {decompose();}

    // follows the C implementation from Wikipedia: https://en.wikipedia.org/wiki/LU_decomposition
    void decompose() {
        Vector row;
        P = Vector(A.N); // enumerate our matrix rows so we can keep track of permutations.
        for (size_t i = 0; i < A.N; ++i) {P[i] = i;}
        permutations = 0;

        double A_max; size_t i_max;
        size_t j, k;
        for (size_t i = 0; i < A.N; ++i) {
            // find pivot point
            A_max = 0; i_max = 0;
            for (k = i; k < A.N; ++k) {
                if (abs(A[k][i]) > A_max) {
                    A_max = abs(A[k][i]);
                    i_max = k;
                }
            }
            if (A_max < precision) {throw std::invalid_argument("Error in LUDecomposition::decompose: Matrix is degenerate.");}
            if (i_max != i) {
                // pivot P
                j = P[i];
                P[i] = P[i_max];
                P[i_max] = j;

                // exchange rows of A
                row = A[i];
                A[i] = A[i_max];
                A[i_max] = row;

                permutations++;
            }

            for (j = i+1; j < A.N; ++j) {
                A[j][i] /= A[i][i];
                for (k = i+1; k < A.N; ++k) {
                    A[j][k] -= A[j][i]*A[i][k];
                }
            }
        }
    }

    double determinant() const {
        double det = A[0][0];
        for (size_t i = 1; i < A.N; ++i) {
            det *= A[i][i];
        }
        return permutations % 2 == 0 ? det : -det;
    }

    int permutations;
private: 
    Vector P;
    Matrix A;
};