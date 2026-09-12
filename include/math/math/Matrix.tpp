#pragma once

#include <math/LUPDecomposition.h>
#include <math/Matrix.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

namespace ausaxs {
    template<numeric Q>
    Matrix<Q>::Matrix(std::initializer_list<std::initializer_list<Q>> l) : N(l.size()), M(l.begin()->size()) {
        for (const auto& row : l) {
            assert(static_cast<int>(row.size()) == M && "Matrix::Matrix: columns must be of equal size!");
            for (const auto& e : row) {
                data.push_back(e);
            }
        }
    }

    template<numeric Q>
    Matrix<Q>::Matrix(const std::vector<std::vector<Q>>& cols) : N(cols[0].size()), M(cols.size()), data(N*M) {
        for (int col = 0; col < M; col++) {
            assert(static_cast<int>(cols[col].size()) == N && "Matrix::Matrix: columns must be of equal size!");
            for (int row = 0; row < N; row++) {
                index(row, col) = cols[col][row];
            }
        }
    }

    template<numeric Q>
    Matrix<Q>::Matrix(const Vector<Q>& v) : N(v.size()), M(1), data(v.data) {}

    template<numeric Q>
    Matrix<Q>::Matrix(int n, int m) : N(n), M(m), data(N*M) {} 

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

    template<numeric Q> ConstRow<Q> Matrix<Q>::operator[](int i) const {return row(i);}
    template<numeric Q> MutableRow<Q> Matrix<Q>::operator[](int i) {return row(i);}

    template<numeric Q> ConstColumn<Q> Matrix<Q>::col(int j) const {return ConstColumn<Q>(data, N, M, j);}
    template<numeric Q> MutableColumn<Q> Matrix<Q>::col(int j) {return MutableColumn<Q>(data, N, M, j);}

    template<numeric Q> ConstRow<Q> Matrix<Q>::row(int i) const {return ConstRow<Q>(data, N, M, i);}
    template<numeric Q> MutableRow<Q> Matrix<Q>::row(int i) {return MutableRow<Q>(data, N, M, i);}

    template<numeric Q> template<numeric R>
    bool Matrix<Q>::operator==(const Matrix<R>& A) const {
        compatibility_check(A);
        Matrix<Q> diff = *this - A; // difference matrix
        return std::accumulate(diff.begin(), diff.end(), 0.0, [] (double sum, Q x) {return sum + std::abs(x);}) < precision;
    }

    template<numeric Q>
    double Matrix<Q>::det() const {
        assert(N == M && "Matrix::det: Matrix is not square.");

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
        for (int row = 0; row < A.N; ++row) {
            for (int col = 0; col < A.M; ++col) {
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
    const Q& Matrix<Q>::operator()(int i, int j) const {
        return index(i, j);
    }

    template<numeric Q>
    Q& Matrix<Q>::operator()(int i, int j) {
        return index(i, j);
    }

    template<numeric Q>
    typename std::vector<Q>::const_iterator Matrix<Q>::begin() const {return data.cbegin();}

    template<numeric Q>
    typename std::vector<Q>::const_iterator Matrix<Q>::end() const {return data.cend();}

    template<numeric Q>
    typename std::vector<Q>::iterator Matrix<Q>::begin() {return data.begin();}

    template<numeric Q>
    typename std::vector<Q>::iterator Matrix<Q>::end() {return data.end();}

    template<numeric Q>
    std::string Matrix<Q>::to_string() const {
        std::stringstream ss;
        for (int i = 0; i < N; i++) {
            ss << "\t" << std::setprecision(3);
            for (int j = 0; j < M; j++) {
                ss << std::setw(8) << index(i, j);
            }
            ss << std::endl;
        }
        return ss.str();
    }

    template<numeric Q> template<numeric R>
    void Matrix<Q>::compatibility_check([[maybe_unused]] const Matrix<R>& A) const {
        assert([&]() -> bool {
            if (N == A.N && M == A.M) {return true;}
            std::cout << "Matrix::compatibility_check: Matrix dimensions do not match (got: [" << N << ", " << M << "] and [" << A.N << ", " << A.M << "])." << std::endl;
            return false;
        }() && "Matrix::compatibility_check: Matrix dimensions do not match.");
    }

    template<numeric Q>
    void Matrix<Q>::compatibility_check_N([[maybe_unused]] int N) const {
        assert([&]() -> bool {
            if (this->N == N) {return true;}
            std::cout << "Matrix::compatibility_check: Matrix dimensions do not match (got: N = " << N << ", expected " << this->N << ")" << std::endl;
            return false;
        }() && "Matrix::compatibility_check: Matrix dimensions do not match.");
    }

    template<numeric Q>
    void Matrix<Q>::compatibility_check_M([[maybe_unused]] int M) const {
        assert([&]() -> bool {
            if (this->M == M) {return true;}
            std::cout << "Matrix::compatibility_check: Matrix dimensions do not match (got: M = " << M << ", expected " << this->M << ")" << std::endl;
            return false;
        }() && "Matrix::compatibility_check: Matrix dimensions do not match.");
    }

    template<numeric Q>
    Matrix<Q> Matrix<Q>::identity(int dim) {
        Matrix A(dim, dim);
        for (int i = 0; i < dim; ++i) {
            A[i][i] = 1;
        }
        return A;
    }
}