#pragma once

#include <initializer_list>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iomanip>

#include "math/Vector.h"

class Matrix {
public: 
    Matrix(std::initializer_list<std::initializer_list<double>> l) : _N(l.size()), _M(l.begin()->size()) {
        _data.assign(l.begin(), l.end());
        for (const auto& e : data) {
            if (__builtin_expect(e.size() != M, false)) {throw std::invalid_argument("Malformed matrix: columns must be of equal size!");}
        }
    }
    Matrix(const Vector& v) : _N(v.N), _M(1), _data(1, std::vector<double>(v.N, 0)) {_data[0] = v.data;}
    Matrix(const int& N, const int& M) : _N(N), _M(M), _data(N, std::vector<double>(M, 0)) {}
    Matrix() : _N(0), _M(0), _data(0, std::vector<double>(0)) {}
    ~Matrix() {}

    // Assignment operator, B = A
    Matrix& operator=(const Matrix& A) {
        _N = A.N; _M = A.M;
        _data = std::vector<std::vector<double>>(A.data);
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
        for (int row = 0; row < N; ++row) {
            for (int col = 0; col < M; ++col) {
                A[row][col] = data[row][col]*a;
            }
        }
        return A;
    }

    // Scalar division, B/a
    Matrix operator/(const double& a) const {
        Matrix A(N, M);
        for (int row = 0; row < N; ++row) {
            for (int col = 0; col < M; ++col) {
                A[row][col] = data[row][col]/a;
            }
        }
        return A;
    }

    // Plus-assignment, B += A
    Matrix& operator+=(const Matrix& A) {
        compatibility_check(A);
        for (int row = 0; row < N; ++row) {
            for (int col = 0; col < M; ++col) {
                _data[row][col] += A[row][col];
            }
        }
        return *this;
    }

    // Minus-assignment, B -= A
    Matrix& operator-=(const Matrix& A) {
        compatibility_check(A);
        for (int row = 0; row < N; ++row) {
            for (int col = 0; col < M; ++col) {
                _data[row][col] -= A[row][col];
            }
        }
        return *this;
    }

    // Read-only indexing, A[i]
    const std::vector<double>& operator[](const int& i) const {return data[i];}
    
    // Read/write indexing, A[i] = ...
    std::vector<double>& operator[](const int& i) {return _data[i];}

    // Vector multiplication, A*v
    friend Vector operator*(const Matrix& A, const Vector& v) {
        if (__builtin_expect(A.M != v.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(v.N) + ", expected: " + std::to_string(A.M) + "]).");
        }
        Vector w(A.N);
        for (int row = 0; row < A.N; ++row) {
            for (int col = 0; col < A.M; ++col) {
                w[row] += v[col]*A[row][col];
            }
        }
        return w;
    }

    // Matrix multiplication, A*B
    friend Matrix operator*(const Matrix& A, const Matrix& B) {
        if (__builtin_expect(A.M != B.N, false)) {
            throw std::invalid_argument("Invalid matrix dimensions (got: " + std::to_string(A.M) + ", " + std::to_string(A.N) + 
                ", expected: " + std::to_string(B.N) + ", " + std::to_string(B.M) + "]).");
        }
        Matrix C(A.N, B.M);
        for (int row = 0; row < C.N; ++row) {
            for (int col = 0; col < C.M; ++col) {
                for (int inner = 0; inner < A.M; ++inner) {
                    C[row][col] += A[row][inner]*B[inner][col];
                }
            }
        }
        return C;
    }

    // Approximate equality, B ~ A
    bool operator==(const Matrix& A) const {
        compatibility_check(A);
        double sum = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                sum += abs(data[i][j] - A[i][j]);
            }
        }
        return sum < precision;
    }

    // Approximate inequality operator, w != v
    bool operator!=(const Matrix& A) const {return !operator==(A);}

    double det() const;

    // Returns a copy of this matrix
    Matrix copy() const {
        Matrix A(N, M);
        A._data = std::vector<std::vector<double>>(data);
        return A;
    }

    // Transpose
    Matrix T() const {
        Matrix A(M, N);
        for (int row = 0; row < A.N; ++row) {
            for (int col = 0; col < A.M; ++col) {
                A[row][col] = data[col][row];
            }
        }
        return A;
    }

    // read-only iterators
    const std::vector<std::vector<double>>::const_iterator begin() const {return data.begin();}
    const std::vector<std::vector<double>>::const_iterator end() const {return data.end();}

    // read-write iterators
    std::vector<std::vector<double>>::iterator begin() {return _data.begin();}
    std::vector<std::vector<double>>::iterator end() {return _data.end();}

    // Print this matrix to the terminal.
    void print() const {
        std::cout << "Printing a (" + std::to_string(N) + ", " + std::to_string(M) + ") matrix: " << std::endl;
        for (const auto& row : data) {
            std::cout << "\t" << std::setprecision(3);
            for (const auto& col : row) {
                std::cout << std::setw(8) << col;
            }
            std::cout << std::endl;
        }
    }

    const int &N = _N, &M = _M; // read-only access to the dimensions
    const std::vector<std::vector<double>>& data = _data; // read-only access to the data container

private: 
    int _N, _M;
    std::vector<std::vector<double>> _data;
    static constexpr double precision = 1e-9;

    // check if the matrix is compatible with ours
    virtual void compatibility_check(const Matrix& A) const {
        if (__builtin_expect(N != A.N || M != A.M, false)) {
            throw std::invalid_argument("Matrix dimensions do not match (got: [" + std::to_string(N) + ", " + std::to_string(M) + "] and [" + 
                std::to_string(A.N) + ", " + std::to_string(A.M) + "]).");
        }
    }
};