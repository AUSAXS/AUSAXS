#include "Matrix.h"
#include "Vector.h"
#include "LUPDecomposition.h"

Matrix::Matrix(const Vector& v) : _N(v.N), _M(1), _data(v.data) {} // vector --> matrix constructor

double Matrix::det() const {
    if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

// Read-only indexing, A[i]
const ConstRowSlice Matrix::operator[](const int& i) const {return ConstRowSlice(this, i);}

// Read/write indexing, A[i] = ...
RowSlice Matrix::operator[](const int& i) {return RowSlice(this, i);}