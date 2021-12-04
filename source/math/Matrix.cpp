#include "Matrix.h"
#include "LUPDecomposition.h"

double Matrix::det() const {
    if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

const Matrix::RowIterator Matrix::row_begin() const {return Matrix::RowIterator(this, 0);}
Matrix::RowIterator Matrix::row_begin() {return Matrix::RowIterator(this, 0);}
const size_t Matrix::row_end() const {return N+1;}
// Matrix::RowIterator Matrix::row_end() {return Matrix::RowIterator(this, N+1);}

const Matrix::ColumnIterator Matrix::col_begin() const {return Matrix::ColumnIterator(this, 0);}
Matrix::ColumnIterator Matrix::col_begin() {return Matrix::ColumnIterator(this, 0);}
const size_t Matrix::col_end() const {return M+1;}
// Matrix::ColumnIterator Matrix::col_end() {return Matrix::ColumnIterator(this, N+1);}

// const Matrix::RowIterator Matrix::row_begin() const {return Matrix::RowIterator(this, 0);}
// Matrix::RowIterator Matrix::row_begin() {return Matrix::RowIterator(this, 0);}
// const Matrix::RowIterator Matrix::row_end() const {return Matrix::RowIterator(this, N+1);}
// Matrix::RowIterator Matrix::row_end() {return Matrix::RowIterator(this, N+1);}

// const Matrix::ColumnIterator Matrix::col_begin() const {return Matrix::ColumnIterator(this, 0);}
// Matrix::ColumnIterator Matrix::col_begin() {return Matrix::ColumnIterator(this, 0);}
// const Matrix::ColumnIterator Matrix::col_end() const {return Matrix::ColumnIterator(this, M+1);}
// Matrix::ColumnIterator Matrix::col_end() {return Matrix::ColumnIterator(this, M+1);}
