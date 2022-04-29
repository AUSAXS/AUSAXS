#pragma once

#include <math/Decomposition.h>
#include <math/Vector.h>

template<typename T> class Matrix;

class LUPDecomposition : public Decomposition {
  public: 
    LUPDecomposition(const Matrix<double>& A);

    // follows the C implementation from Wikipedia: https://en.wikipedia.org/wiki/LU_decomposition
    void decompose();

    double determinant() const;

    int permutations;
  private: 
    Vector<double> P;
    std::unique_ptr<Matrix<double>> Ap;
};