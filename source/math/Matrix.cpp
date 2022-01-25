#include "math/Matrix.h"
#include "math/Vector.h"
#include "math/LUPDecomposition.h"

Matrix::Matrix(const Vector& v) : _N(v.N), _M(1), _data(v.data) {} // vector --> matrix constructor

double Matrix::det() const {
    if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

Row Matrix::row(const int i) {return Row(_data, N, M, i);}
Row Matrix::operator[](const int i) {return Row(_data, N, M, i);}
Column Matrix::col(const int i) {return Column(_data, N, M, i);}

const ConstRow Matrix::row(const int i) const {return ConstRow(data, N, M, i);}
const ConstRow Matrix::operator[](const int i) const {return ConstRow(data, N, M, i);}
const ConstColumn Matrix::col(const int i) const {return ConstColumn(data, N, M, i);}