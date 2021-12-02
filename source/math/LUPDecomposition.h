#pragma once

// forwards declaration
class Matrix;

#include "Decomposition.h"
#include "Matrix.h"
#include "Vector.h"

#include <stdexcept>

class LUPDecomposition : public Decomposition {
public: 
    LUPDecomposition(const Matrix& A) : A(A.copy()) {}

    // follows the C implementation from Wikipedia: https://en.wikipedia.org/wiki/LU_decomposition
    void decompose() {
        int N = A.N;

        std::vector<double> row;
        P = Vector(A.N); // enumerate our matrix rows so we can keep track of permutations.
        for (int i = 0; i < N; ++i) {P[i] = i;}
        permutations = 0;

        double A_max; int i_max;
        int j, k;
        for (int i = 0; i < N; ++i) {
            // find pivot point
            A_max = 0; i_max = 0;
            for (k = i; k < N; ++k) {
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

            for (j = i+1; j < N; ++j) {
                A[j][i] /= A[i][i];
                for (k = i+1; k < N; ++k) {
                    A[j][k] -= A[j][i]*A[i][k];
                }
            }
        }
    }

    double determinant() {
        if (P.size() == 0) {decompose();}
        double det = A[0][0];
        for (int i = 1; i < A.N; ++i) {
            det *= A[i][i];
        }
        return permutations % 2 == 0 ? det : -det;
    }

    int permutations;
private: 
    Vector P;
    Matrix A;
};