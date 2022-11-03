#pragma once

#include <math/Matrix.h>

template<typename Q>
Matrix<Q>::Matrix(Matrix<Q>&& A) noexcept : N(A.N), M(A.M), data(std::move(A.data)) {}

template<typename Q>
Matrix<Q>::Matrix(const Matrix<Q>& A) : N(A.N), M(A.M), data(A.data) {}

template<typename Q>
Matrix<Q>::Matrix(std::initializer_list<std::initializer_list<Q>> l) : N(l.size()), M(l.begin()->size()) {
    for (const auto& row : l) {
        if (__builtin_expect(row.size() != M, false)) {throw std::invalid_argument("Malformed matrix: columns must be of equal size!");}
        for (const auto& e : row) {
            data.push_back(e);
        }
    }
}

template<typename Q>
Matrix<Q>::Matrix(std::vector<std::vector<Q>> v) : N(v[0].size()), M(v.size()), data(N*M) {
    for (unsigned int col = 0; col < M; col++) {
        if (__builtin_expect(v[col].size() != N, false)) {throw std::invalid_argument("Malformed matrix: columns must be of equal size!");}
        for (unsigned int row = 0; row < N; row++) {
            index(row, col) = v[col][row];
        }
    }
}

template<typename Q>
Matrix<Q>::Matrix(const Vector<Q>& v) : N(v.N), M(1), data(v.data) {}

template<typename Q>
Matrix<Q>::Matrix(int n, int m) : N(n), M(m), data(N*M) {} 

template<typename Q>
Matrix<Q>::Matrix() : N(0), M(0), data(0) {}

template<typename Q>
Matrix<Q>::~Matrix() = default;

template<typename Q>
void Matrix<Q>::push_back(std::vector<double> r) {
    compatibility_check_M(r.size());
    extend(1);
    row(N-1) = r;
}

template<typename Q>
Matrix<Q>& Matrix<Q>::operator*=(double a) {
    std::transform(begin(), end(), begin(), [a] (const Q& x) {return x*a;});
    return *this;
}

template<typename Q>
Matrix<Q>& Matrix<Q>::operator/=(double a) {
    std::transform(begin(), end(), begin(), [a] (const Q& x) {return x/a;});
    return *this;
}

template<typename Q>
Matrix<Q>& Matrix<Q>::operator=(const Matrix<Q>& A) {
    N = A.N; M = A.M;
    data = A.data;
    return *this;
}

template<typename Q>
Matrix<Q> Matrix<Q>::operator-() const {
    Matrix<Q> A(N, M);
    std::transform(begin(), end(), A.begin(), std::negate<Q>());
    return A;
}

template<typename Q> template<typename R>
Matrix<Q>& Matrix<Q>::operator+=(const Matrix<R>& A) {
    compatibility_check(A);
    std::transform(begin(), end(), A.begin(), begin(), std::plus<Q>());
    return *this;
}

template<typename Q> template<typename R>
Matrix<Q>& Matrix<Q>::operator-=(const Matrix<R>& A) {
    compatibility_check(A);
    std::transform(begin(), end(), A.begin(), begin(), std::minus<Q>());
    return *this;
}

template<typename Q>
void Matrix<Q>::extend(int n) {
    N += n;
    data.resize(N*M);
}

template<typename Q>
void Matrix<Q>::resize(int n, int m) {
    N = n; M = m;
    data.resize(N*M);
}

template<typename Q> const ConstRow<Q> Matrix<Q>::operator[](unsigned int i) const {return row(i);}
template<typename Q> Row<Q> Matrix<Q>::operator[](unsigned int i) {return row(i);}

template<typename Q> const ConstColumn<Q> Matrix<Q>::col(unsigned int j) const {return ConstColumn<Q>(data, N, M, j);}
template<typename Q> Column<Q> Matrix<Q>::col(unsigned int j) {return Column<Q>(data, N, M, j);}

template<typename Q> const ConstRow<Q> Matrix<Q>::row(unsigned int i) const {return ConstRow<Q>(data, N, M, i);}
template<typename Q> Row<Q> Matrix<Q>::row(unsigned int i) {return Row<Q>(data, N, M, i);}

template<typename Q> template<typename R>
bool Matrix<Q>::operator==(const Matrix<R>& A) const {
    compatibility_check(A);
    Matrix<Q> diff = *this - A; // difference matrix
    return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, Q x) {return sum + abs(x);}) < precision;
}

template<typename Q> template<typename R>
bool Matrix<Q>::operator!=(const Matrix<R>& A) const {return !operator==(A);}

template<typename Q>
double Matrix<Q>::det() const {
    #if (SAFE_MATH)
        if (__builtin_expect(N != M, false)) {throw std::invalid_argument("Error in matrix determinant: Matrix is not square.");}
    #endif

    LUPDecomposition decomp(*this);
    return decomp.determinant();
}

template<typename Q>
Matrix<Q> Matrix<Q>::copy() const {
    Matrix A(N, M);
    A.data.assign(data.begin(), data.end());
    return A;
}

template<typename Q>
Matrix<Q> Matrix<Q>::T() const {
    Matrix A(M, N);
    for (size_t row = 0; row < A.N; ++row) {
        for (size_t col = 0; col < A.M; ++col) {
            A[row][col] = index(col, row);
        }
    }
    return A;
}

template<typename Q>
Matrix<Q> Matrix<Q>::transpose() const {
    return T();
}

template<typename Q>
const Q& Matrix<Q>::index(unsigned int i, unsigned int j) const {
    #if (SAFE_MATH)
        if (__builtin_expect(i >= N, false)) {raise(SIGSEGV);}
        if (__builtin_expect(j >= M, false)) {raise(SIGSEGV);}
        // if (__builtin_expect(i >= N, false)) {throw std::out_of_range("Error in Matrix::index: Row index out of range.");}
        // if (__builtin_expect(j >= M, false)) {throw std::out_of_range("Error in Matrix::index: Column index out of range.");}
    #endif
    return data[M*i + j];
}

template<typename Q>
Q& Matrix<Q>::index(unsigned int i, unsigned int j) {
    #if (SAFE_MATH)
        if (__builtin_expect(i >= N, false)) {raise(SIGSEGV);}
        if (__builtin_expect(j >= M, false)) {raise(SIGSEGV);}
        // if (__builtin_expect(i >= N, false)) {throw std::out_of_range("Error in Matrix::index: Row index out of range.");}
        // if (__builtin_expect(j >= M, false)) {throw std::out_of_range("Error in Matrix::index: Column index out of range.");}
    #endif
    return data[M*i + j];
}

template<typename Q>
const typename std::vector<Q>::const_iterator Matrix<Q>::begin() const {return data.begin();}

template<typename Q>
const typename std::vector<Q>::const_iterator Matrix<Q>::end() const {return data.end();}

template<typename Q>
typename std::vector<Q>::iterator Matrix<Q>::begin() {return data.begin();}

template<typename Q>
typename std::vector<Q>::iterator Matrix<Q>::end() {return data.end();}

template<typename Q>
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

template<typename Q> template<typename R>
void Matrix<Q>::compatibility_check(const Matrix<R>& A) const {
    #if (SAFE_MATH)
        if (__builtin_expect(N != A.N || M != A.M, false)) {
            throw std::invalid_argument("Matrix dimensions do not match (got: [" + std::to_string(N) + ", " + std::to_string(M) + "] and [" + 
                std::to_string(A.N) + ", " + std::to_string(A.M) + "]).");
        }
    #endif
}

template<typename Q>
void Matrix<Q>::compatibility_check_N(unsigned int N) const {
    #if (SAFE_MATH)
        if (__builtin_expect(this->N != N, false)) {
            throw std::invalid_argument("Matrix dimensions do not match (got: N = " + std::to_string(N) + ", expected " + std::to_string(this->N) + ")");
        }
    #endif
}

template<typename Q>
void Matrix<Q>::compatibility_check_M(unsigned int M) const {
    #if (SAFE_MATH)
        if (__builtin_expect(this->M != M, false)) {
            throw std::invalid_argument("Matrix dimensions do not match (got: M = " + std::to_string(N) + ", expected " + std::to_string(this->N) + ")");
        }
    #endif
}