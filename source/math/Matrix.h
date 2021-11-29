#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "math/Vector.h"

class Matrix {
public: 
    Matrix(std::initializer_list<std::initializer_list<double>> l) : N(N), M(M) {
        data.assign(l.begin(), l.end());
        for (const auto& e : data) {
            if (__builtin_expect(e.size() != data[0].size(), false)) {throw std::invalid_argument("Matrix columns must be of equal size!");}
        }
    }
    Matrix(const int& N, const int& M) : N(N), M(M), data(N, std::vector<double>(M, 0)) {}
    ~Matrix() {}

    // Assignment operator, B = A
    Matrix& operator=(const Matrix& A) {
        compatibility_check(A);
        data.assign(A.begin(), A.end());
        return *this;
    }

    // Plus operator, B + A
    Matrix operator+(const Matrix& A) const {
        compatibility_check(A);
        Matrix B(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B[i][j] = data[i][j] + A[i][j];
            }
        }
        return B;
    }

    // Minus operator, B - A
    Matrix operator-(const Matrix& A) const {
        compatibility_check(A);
        Matrix B(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                B[i][j] = data[i][j] - A[i][j];
            }
        }
        return B;
    }

    // Negation operator, -A
    Matrix operator-() const {
        Matrix A(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                A[i][j] = -data[i][j];
            }
        }
        return A;
    }

    // Scalar multiplication, B*a
    Matrix operator*(const double& a) const {
        Matrix A(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                A[i][j] = a*data[i][j];
            }
        }
        return A;
    }

    // // Scalar division, w/a
    // Vector operator/(const double& a) const {
    //     Vector w(N);
    //     std::transform(begin(), end(), w.begin(), [&a](double x) {return x/a;});
    //     return w;
    // }

    // // Plus-assignment, w += v
    // Vector& operator+=(const Vector& v) {
    //     compatibility_check(v);
    //     std::transform(begin(), end(), v.begin(), begin(), std::plus<double>()); 
    //     return *this;
    // }

    // // Minus-assignment, w -= v
    // Vector& operator-=(const Vector& v) {
    //     compatibility_check(v);
    //     std::transform(begin(), end(), v.begin(), begin(), std::minus<double>()); 
    //     return *this;
    // }

    // Read-only indexing, w[i]
    const std::vector<double>& operator[](const int& i) const {return data[i];}
    
    // Read/write indexing, w[i] = ...
    std::vector<double>& operator[](const int& i) {return data[i];}

    // Vector multiplication, A*v
    friend double operator*(const Matrix& A, const Vector& v) {
        if (__builtin_expect(A.M != v.N, false)) {
            throw std::invalid_argument("Invalid vector dimension (got: " + std::to_string(v.N) + ", expected: " + std::to_string(A.M) + " ).");
        }
        double a = 0;
        for (int i = 0; i < A.N; ++i) {
            for (int j = 0; j < A.M; ++j) {
                a += v[j]*A[i][j];
            }
        }
        return a;
    }

    // Matrix multiplication, A*B
    friend Matrix operator*(const Matrix& A, const Matrix& B) {
        if (__builtin_expect(A.M != B.N || A.N != B.M, false)) {
            throw std::invalid_argument("Invalid vector dimension (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + 
                ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + " ).");
        }
        double a = 0;
        for (int i = 0; i < A.N; ++i) {
            for (int j = 0; j < A.M; ++j) {
                a += v[j]*A[i][j];
            }
        }
        return a;
    }

    // // Approximate equality, w ~ v
    // bool operator==(const Vector& v) const {
    //     compatibility_check(v);
    //     Vector a = operator-(v); // difference vector
    //     return std::accumulate(a.begin(), a.end(), 0.0, [] (double sum, double x) {return sum + abs(x);}) < precision;
    // }

    // // Approximate inequality operator, w != v
    // bool operator!=(const Vector& v) const {return !operator==(v);}

    // check if the matrix is compatible with ours
    virtual void compatibility_check(const Matrix& A) const {
        if (__builtin_expect(N != A.N || M != A.M, false)) {throw std::invalid_argument("Vector dimensions do not match.");}
    }

    // read-only iterators
    const std::vector<std::vector<double>>::const_iterator begin() const {return data.begin();}
    const std::vector<std::vector<double>>::const_iterator end() const {return data.end();}

    // read-write iterators
    std::vector<std::vector<double>>::iterator begin() {return data.begin();}
    std::vector<std::vector<double>>::iterator end() {return data.end();}

    const int N, M;

private: 
    std::vector<std::vector<double>> data;
    static constexpr double precision = 1e-9;
};