#pragma once

#include <math/Matrix.h>
#include <math/LUPDecomposition.h>
#include <utility/Exceptions.h>

#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>

template<numeric Q>
Matrix<Q>::Matrix(Matrix<Q>&& A) noexcept : N(A.N), M(A.M), data(std::move(A.data)) {}

template<numeric Q>
Matrix<Q>::Matrix(const Matrix<Q>& A) : N(A.N), M(A.M), data(A.data) {}

template<numeric Q>
Matrix<Q>::Matrix(std::initializer_list<std::initializer_list<Q>> l) : N(l.size()), M(l.begin()->size()) {
    for (const auto& row : l) {
        if (row.size() != M) [[unlikely]] {throw except::invalid_argument("Matrix::Matrix: columns must be of equal size!");}
        for (const auto& e : row) {
            data.push_back(e);
        }
    }
}

template<numeric Q>
Matrix<Q>::Matrix(const std::vector<std::vector<Q>>& cols) : N(cols[0].size()), M(cols.size()), data(N*M) {
    for (unsigned int col = 0; col < M; col++) {
        if (cols[col].size() != N) [[unlikely]] {throw except::invalid_argument("Matrix::Matrix: columns must be of equal size!");}
        for (unsigned int row = 0; row < N; row++) {
            index(row, col) = cols[col][row];
        }
    }
}

template<numeric Q>
Matrix<Q>::Matrix(const Vector<Q>& v) : N(v.N), M(1), data(v.data) {}

template<numeric Q>
Matrix<Q>::Matrix(unsigned int n, unsigned int m) : N(n), M(m), data(N*M) {} 

template<numeric Q>
Matrix<Q>::Matrix() : N(0), M(0), data(0) {}

template<numeric Q>
Matrix<Q>::~Matrix() = default;

template<numeric Q>
void Matrix<Q>::push_back(const std::vector<double>& r) {
    compatibility_check_M(r.size());
    extend(1);
    row(N-1) = r;
}

template<numeric Q>
Matrix<Q>& Matrix<Q>::operator*=(double a) {
    std::transform(begin(), end(), begin(), [a] (const Q& x) {return x*a;});
    return *this;
}

template<numeric Q>
Matrix<Q>& Matrix<Q>::operator/=(double a) {
    std::transform(begin(), end(), begin(), [a] (const Q& x) {return x/a;});
    return *this;
}

template<numeric Q>
Matrix<Q>& Matrix<Q>::operator=(const Matrix<Q>& A) {
    N = A.N; M = A.M;
    data = A.data;
    return *this;
}

template<numeric Q>
Matrix<Q>& Matrix<Q>::operator=(Matrix<Q>&& A) {
    N = A.N; M = A.M;
    data = std::move(A.data);
    return *this;
}

template<numeric Q>
Matrix<Q> Matrix<Q>::operator-() const {
    Matrix<Q> A(N, M);
    std::transform(begin(), end(), A.begin(), std::negate<Q>());
    return A;
}

template<numeric Q> template<numeric R>
Matrix<Q>& Matrix<Q>::operator+=(const Matrix<R>& A) {
    compatibility_check(A);
    std::transform(begin(), end(), A.begin(), begin(), std::plus<Q>());
    return *this;
}

template<numeric Q> template<numeric R>
Matrix<Q>& Matrix<Q>::operator-=(const Matrix<R>& A) {
    compatibility_check(A);
    std::transform(begin(), end(), A.begin(), begin(), std::minus<Q>());
    return *this;
}

template<numeric Q>
void Matrix<Q>::extend(int n) {
    N += n;
    data.resize(N*M);
}

template<numeric Q>
void Matrix<Q>::resize(int n, int m) {
    N = n; M = m;
    data.resize(N*M);
}

template<numeric Q> const ConstRow<Q> Matrix<Q>::operator[](unsigned int i) const {return row(i);}
template<numeric Q> MutableRow<Q> Matrix<Q>::operator[](unsigned int i) {return row(i);}

template<numeric Q> const ConstColumn<Q> Matrix<Q>::col(unsigned int j) const {return ConstColumn<Q>(data, N, M, j);}
template<numeric Q> MutableColumn<Q> Matrix<Q>::col(unsigned int j) {return MutableColumn<Q>(data, N, M, j);}

template<numeric Q> const ConstRow<Q> Matrix<Q>::row(unsigned int i) const {return ConstRow<Q>(data, N, M, i);}
template<numeric Q> MutableRow<Q> Matrix<Q>::row(unsigned int i) {return MutableRow<Q>(data, N, M, i);}

template<numeric Q> template<numeric R>
bool Matrix<Q>::operator==(const Matrix<R>& A) const {
    compatibility_check(A);
    Matrix<Q> diff = *this - A; // difference matrix
    return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, Q x) {return sum + std::abs(x);}) < precision;
}

template<numeric Q>
double Matrix<Q>::det() const {
    #if (SAFE_MATH)
        if (N != M) [[unlikely]] {throw except::invalid_argument("Matrix::det: Matrix is not square.");}
    #endif

    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

template<numeric Q>
Matrix<Q> Matrix<Q>::copy() const {
    Matrix A(N, M);
    A.data.assign(begin(), end());
    return A;
}

template<numeric Q>
Matrix<Q> Matrix<Q>::T() const {
    Matrix A(M, N);
    for (unsigned int row = 0; row < A.N; ++row) {
        for (unsigned int col = 0; col < A.M; ++col) {
            A[row][col] = index(col, row);
        }
    }
    return A;
}

template<numeric Q>
Matrix<Q> Matrix<Q>::transpose() const {
    return T();
}

template<numeric Q>
const Q& Matrix<Q>::operator()(unsigned int i, unsigned int j) const {
    return index(i, j);
}

template<numeric Q>
Q& Matrix<Q>::operator()(unsigned int i, unsigned int j) {
    return index(i, j);
}

template<numeric Q>
const Q& Matrix<Q>::index(unsigned int i, unsigned int j) const {
    #if (SAFE_MATH)
        if (i >= N) [[unlikely]] {throw except::out_of_bounds("Matrix::index: Row index out of range. (" + std::to_string(i) + " >= " + std::to_string(N) + ")");}
        if (j >= M) [[unlikely]] {throw except::out_of_bounds("Matrix::index: Column index out of range. (" + std::to_string(j) + " >= " + std::to_string(M) + ")");}
    #endif
    return data[M*i + j];
}

template<numeric Q>
Q& Matrix<Q>::index(unsigned int i, unsigned int j) {
    #if (SAFE_MATH)
        if (i >= N) [[unlikely]] {throw except::out_of_bounds("Matrix::index: Row index out of range. (" + std::to_string(i) + " >= " + std::to_string(N) + ")");}
        if (j >= M) [[unlikely]] {throw except::out_of_bounds("Matrix::index: Column index out of range. (" + std::to_string(j) + " >= " + std::to_string(M) + ")");}
    #endif
    return data[M*i + j];
}

template<numeric Q>
const typename std::vector<Q>::const_iterator Matrix<Q>::begin() const {return data.cbegin();}

template<numeric Q>
const typename std::vector<Q>::const_iterator Matrix<Q>::end() const {return data.cend();}

template<numeric Q>
typename std::vector<Q>::iterator Matrix<Q>::begin() {return data.begin();}

template<numeric Q>
typename std::vector<Q>::iterator Matrix<Q>::end() {return data.end();}

template<numeric Q>
std::string Matrix<Q>::to_string() const {
    std::stringstream ss;
    for (unsigned int i = 0; i < N; i++) {
        ss << "\t" << std::setprecision(3);
        for (unsigned int j = 0; j < M; j++) {
            ss << std::setw(8) << index(i, j);
        }
        ss << std::endl;
    }
    return ss.str();
}

template<numeric Q> template<numeric R>
void Matrix<Q>::compatibility_check(const Matrix<R>& A) const {
    #if (SAFE_MATH)
        if (N != A.N || M != A.M) [[unlikely]] {
            throw except::invalid_argument("Matrix::compatibility_check: Matrix dimensions do not match (got: [" + std::to_string(N) + ", " + std::to_string(M) + "] and [" + 
                std::to_string(A.N) + ", " + std::to_string(A.M) + "]).");
        }
    #endif
}

template<numeric Q>
void Matrix<Q>::compatibility_check_N(unsigned int N) const {
    #if (SAFE_MATH)
        if (this->N != N) [[unlikely]] {
            throw except::invalid_argument("Matrix::compatibility_check: Matrix dimensions do not match (got: N = " + std::to_string(N) + ", expected " + std::to_string(this->N) + ")");
        }
    #endif
}

template<numeric Q>
void Matrix<Q>::compatibility_check_M(unsigned int M) const {
    #if (SAFE_MATH)
        if (this->M != M) [[unlikely]] {
            throw except::invalid_argument("Matrix::compatibility_check: Matrix dimensions do not match (got: M = " + std::to_string(N) + ", expected " + std::to_string(this->N) + ")");
        }
    #endif
}

template<numeric Q>
Matrix<Q> Matrix<Q>::identity(unsigned int dim) {
    Matrix A(dim, dim);
    for (unsigned int i = 0; i < dim; ++i) {
        A[i][i] = 1;
    }
    return A;
}
