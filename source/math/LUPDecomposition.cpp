#include <math/LUPDecomposition.h>
#include <math/Matrix.h>
#include <math/Vector.h>

LUPDecomposition::LUPDecomposition(const Matrix<double>& A) : Ap(std::make_unique<Matrix<double>>(A.copy())) {decompose();}

void LUPDecomposition::decompose() {
    Matrix<double>& A = *Ap;

    Vector<double> row;
    P = Vector<double>(A.N); // enumerate our matrix rows so we can keep track of permutations.
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

double LUPDecomposition::determinant() const {
    Matrix<double>& A = *Ap;

    double det = A[0][0];
    for (size_t i = 1; i < A.N; ++i) {
        det *= A[i][i];
    }
    return permutations % 2 == 0 ? det : -det;
}